function fig6_oddball()
% Figure 6 — Empirical auditory oddball ERP.
%
% Loads pre-extracted deviant trials from cfg.data_dir:
%   oddball_trials_80x351.csv   (N_trials x N_samples, baseline-corrected, uV)
%   oddball_t_ms.csv            (N_samples, time axis in ms)

cfg = n1f.config();

trials_file = fullfile(cfg.data_dir, 'oddball_trials_80x351.csv');
t_ms_file   = fullfile(cfg.data_dir, 'oddball_t_ms.csv');
if ~exist(trials_file,'file') || ~exist(t_ms_file,'file')
    fprintf('  Oddball CSVs not found in %s, skipping fig6.\n', cfg.data_dir);
    return;
end

fprintf('  Loading pre-extracted oddball trials...\n');
trials = readmatrix(trials_file);   % N_trials x N_samples
t_ms   = readmatrix(t_ms_file);
t_ms   = t_ms(:);                   % column
t_s    = t_ms / 1000;

fs_eeg = 1000;
N      = size(trials, 1);
avg    = mean(trials, 1)';
sem    = std(trials, 0, 1)' / sqrt(N);
fprintf('  N = %d trials\n', N);

[fw, psd_avg] = n1f.compute_psd(avg, fs_eeg, numel(avg));
sp = n1f.run_specparam(fw, psd_avg, [2 40]);

fprintf('  Running ERP-SPRiNT...\n');
freqs_an = (2:40)';
res = n1f.run_erp_sprint(trials, t_s', fs_eeg, 'n_cycles',4, 'time_step',0.010, ...
                          'freq_range',[2 40], 'freqs',freqs_an);

% ── Baseline normalise the raw trial-averaged TFR (dB re pre-stim mean) ──
% Frequencies whose wavelet doesn't fully fit inside the pre-stim window
% don't have a usable baseline → row stays NaN and renders white.
fwhm_s = @(f) res.params.n_cycles * sqrt(2*log(2)) / (pi * f);
bl_t   = res.times < 0;
t_pre  = 0 - t_s(1);                       % pre-stim length (s)

tfr_db = nan(size(res.tfr_avg));
for fi = 1:numel(res.freqs)
    if fwhm_s(res.freqs(fi))/2 > t_pre, continue; end
    row  = res.tfr_avg(fi,:);
    bl_v = row(bl_t & ~isnan(row));
    if numel(bl_v) < 3, continue; end
    tfr_db(fi,:) = 10*log10(row ./ mean(bl_v));
end

% Edge-mask each row by its own wavelet FWHM/2
for fi = 1:numel(res.freqs)
    hw = fwhm_s(res.freqs(fi))/2;
    tfr_db(fi, res.times < (t_s(1)+hw) | res.times > (t_s(end)-hw)) = NaN;
end

col = [.77 .20 .21];
fig = figure('Units','centimeters','Position',[0 0 18 16],'Color','w');

% A: ERP
ax = subplot(2,2,1); hold(ax,'on');
fill(ax,[t_ms;flipud(t_ms)],[avg-sem;flipud(avg+sem)],col,'FaceAlpha',.22,'EdgeColor','none');
plot(ax,t_ms,avg,'Color',col,'LineWidth',1.6);
xline(ax,0,'--k','LineWidth',.8); yline(ax,0,':k','LineWidth',.4);
xlabel(ax,'Time (ms)'); ylabel(ax,'Amplitude (uV)'); xlim(ax,[t_ms(1) t_ms(end)]);
title(ax,sprintf('Grand average ERP +/- SEM (N=%d)',N), ...
      'FontSize',8,'FontWeight','bold','HorizontalAlignment','center','Color','k');
set_ax(ax); plabel(ax,'A');

% B: PSD
ax = subplot(2,2,2);
semilogy(ax,fw,psd_avg,'Color',[.75 .75 .75],'LineWidth',.6); hold(ax,'on');
semilogy(ax,sp.freqs,10.^sp.modeled,  'Color',col,'LineWidth',1.5);
semilogy(ax,sp.freqs,10.^sp.aperiodic,'Color',col,'LineWidth',1.0,'LineStyle','--');
xlim(ax,[2 40]); xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power');
lgd = legend(ax,sprintf('exp = %+.2f',sp.exponent), ...
       'Location','northeast','FontSize',7,'Box','off','TextColor','k');

title(ax,'PSD', ...
      'FontSize',8,'FontWeight','bold','Color','k');

text(ax,0.93,0.82,sprintf('R^2 = %.2f',sp.r_squared), ...
     'Units','normalized', ...
     'HorizontalAlignment','right', ...
     'VerticalAlignment','top', ...
     'FontSize',7,'Color','k');
set_ax(ax); plabel(ax,'B');

% C: trial-averaged wavelet TFR, baseline-normalised (dB)
ax = subplot(2,2,3); set(ax,'Color','w');
t_bin_ms = res.times*1000;
fr_msk   = res.freqs >= 2 & res.freqs <= 40;
h = imagesc(ax, t_bin_ms, log2(res.freqs(fr_msk)), tfr_db(fr_msk,:));
set(h,'AlphaData', ~isnan(tfr_db(fr_msk,:)));
set(ax,'YDir','normal','Color','w');
colormap(ax,n1f.redblue(256));
clim_v = max(abs(prctile(tfr_db(:),[0.5 99.5])),[],'omitnan');
if isempty(clim_v) || ~isfinite(clim_v) || clim_v==0, clim_v = 6; end
clim(ax,[-clim_v clim_v]);
yticks_hz = [2 4 8 16 32 40];
set(ax,'YTick',log2(yticks_hz),'YTickLabel',arrayfun(@num2str,yticks_hz,'UniformOutput',false));
cb = colorbar(ax);
cb.FontSize = 7;
cb.Color = 'k';              % tick labels + axis color
cb.Label.String = 'Power (dB re baseline)';
cb.Label.Color = 'k';        % label textxline(ax,0,'--k','LineWidth',.8);
xlabel(ax,'Time (ms)'); ylabel(ax,'Frequency (Hz)'); xlim(ax,[t_ms(1) t_ms(end)]);
title(ax,'Trial-averaged wavelet TFR', ...
      'FontSize',8,'FontWeight','bold','HorizontalAlignment','center','Color','k');
set_ax(ax); plabel(ax,'C');

% D: exponent timecourse (y-axis on right)
ax = subplot(2,2,4); set(ax,'Color','w');
v = ~isnan(res.exponent);
plot(ax,t_bin_ms(v),res.exponent(v),'Color',[.56 .27 .68],'LineWidth',1.5);
hold(ax,'on');
xline(ax,0,'--k','LineWidth',.8);
pre_v = v & res.times < 0;
if any(pre_v)
    exp_bl = mean(res.exponent(pre_v),'omitnan');
    yline(ax,exp_bl,':','Color',[.56 .27 .68],'LineWidth',.7);
    text(ax,t_ms(end)-5,exp_bl+0.006,sprintf('BL=%.2f',exp_bl), ...
         'HorizontalAlignment','right','FontSize',7,'Color',[.56 .27 .68]);
end
xlabel(ax,'Time (ms)'); xlim(ax,[t_ms(1) t_ms(end)]);
set_ax(ax);
set(ax,'YAxisLocation','right','Box','off','TickDir','out','FontSize',8,'LineWidth',.7, ...
       'XColor','k');
ylabel(ax,'Aperiodic exponent'); ax.YColor = 'k';
title(ax,'Aperiodic exponent(t)', ...
      'FontSize',8,'FontWeight','bold','HorizontalAlignment','center','Color','k');
plabel(ax,'D');

save_fig(fig, cfg, 'fig6_oddball');

%% ── Caption statistics ───────────────────────────────────────────────────
[~, n1_idx] = min(avg);
fprintf('\n--- Fig6 Caption statistics ---\n');
fprintf('N1 amplitude = %.1f uV at %.0f ms\n', avg(n1_idx), t_ms(n1_idx));
fprintf('PSD: exp = %.2f  R2 = %.2f  npeaks = %d\n', sp.exponent, sp.r_squared, sp.n_peaks);
if sp.n_peaks > 0
    [~,pk_idx] = max([sp.peaks.pw]);
    fprintf('Dominant peak: %.1f Hz\n', sp.peaks(pk_idx).cf);
end

v    = ~isnan(res.exponent);
bl_v = v & res.times < 0;
fprintf('ERP-SPRiNT unmasked range: %.0f to %.0f ms\n', ...
    min(res.times(v))*1000, max(res.times(v))*1000);
if any(bl_v)
    exp_bl   = nanmean(res.exponent(bl_v));
    n1_win   = res.times >= 0.050 & res.times <= 0.100;
    exp_n1   = nanmean(res.exponent(n1_win));
    fprintf('ERP-SPRiNT exponent: baseline = %.2f  N1 window (50-100ms) = %.2f  delta = %.2f\n', ...
        exp_bl, exp_n1, exp_n1 - exp_bl);
    fprintf('ERP-SPRiNT mean R2: %.2f\n', nanmean(res.r_squared(v)));
end
fprintf('-------------------------------\n\n');
end


%% ── Local helpers ──────────────────────────────────────────────────────────
function set_ax(ax)
set(ax,'Box','off','TickDir','out','FontSize',8,'LineWidth',.7,'Color','w', ...
    'XColor','k','YColor','k');
ax.Title.Color = 'k';
end
function plabel(ax,letter)
text(ax,-0.11,1.05,letter,'Units','normalized','FontSize',12,'FontWeight','bold', ...
     'VerticalAlignment','bottom','Color','k');
end
function save_fig(fig,cfg,name)
fname=fullfile(cfg.out_dir,[name '.png']);
exportgraphics(fig,fname,'Resolution',150); fprintf('  Saved: %s\n',fname);% close(fig);
end