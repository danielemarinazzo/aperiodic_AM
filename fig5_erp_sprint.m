function fig5_erp_sprint()
% alpha lateralisation simulation (specparam on trial-averaged
% spectrogram). 2x2 layout, contra (ERD) vs ipsi (ERS):
%   A  trial-averaged signal
%   B  aperiodic exponent timecourse
%   C  specparam-modelled spectrum (contra), dB re baseline
%   D  specparam-modelled spectrum (ipsi),   dB re baseline
%   C/D overlay: fitted peak centre frequency (black) vs true carrier (dotted)
% ERD = event-related desynchronisation (amplitude decrease);
% ERS = event-related synchronisation (amplitude increase).
%
% Signals: pure 10 Hz oscillation, no noise, no frequency jitter, carrier
% phase drawn per trial and shared across conditions. Long time base
% (t_long) to avoid wavelet edge effects; display window [-0.8, 2.0] s.
% A constant-amplitude control is still computed (res_ct) for the printed
% statistics but is not plotted.
%
% The modelled-spectrum panels (C/D) show how the specparam fit reorganises
% under strong modulation: the fitted alpha peak drifts in centre frequency
% and its power does not track the true amplitude change (see §3.4/§4.3).

cfg = n1f.config();
fs  = cfg.fs;

txt_col = [0.10 0.10 0.10];
cc  = [.16 .50 .73];   % contra / ERD : blue
ic  = [.77 .20 .21];   % ipsi         : red
ctl = [.50 .50 .50];   % control      : grey

%% ── Parameters ───────────────────────────────────────────────────────────
n_trials  = 60;
f0        = 10;
erd_start = 0.35;
erd_end   = 1.25;
t_pad     = 2.0;
t_ep      = (-0.8 : 1/fs : 2.0)';
t_long    = (-(0.8+t_pad) : 1/fs : (2.0+t_pad))';
n_long    = numel(t_long);

%% ── Amplitude envelopes ──────────────────────────────────────────────────
idx = t_long >= erd_start & t_long <= erd_end;
hw  = hann(sum(idx));

amp_contra = ones(n_long,1); amp_contra(idx) = 1 - 0.60*hw;  % ERD
amp_ipsi   = ones(n_long,1); amp_ipsi(idx)   = 1 + 0.60*hw;  % ERS
amp_ctrl   = ones(n_long,1);                                   % constant

%% ── Generate trials ──────────────────────────────────────────────────────
rng(0);
sig_contra = zeros(n_trials, n_long);
sig_ipsi   = zeros(n_trials, n_long);
sig_ctrl   = zeros(n_trials, n_long);

for tr = 1:n_trials
    ph = 2*pi*rand;
    sig_contra(tr,:) = (amp_contra .* cos(2*pi*f0*t_long + ph))';
    sig_ipsi(tr,:)   = (amp_ipsi   .* cos(2*pi*f0*t_long + ph))';
    sig_ctrl(tr,:)   = (amp_ctrl   .* cos(2*pi*f0*t_long + ph))';
end

%% ── ERP-SPRiNT ───────────────────────────────────────────────────────────
sprint_args = {'n_cycles',6,'time_step',0.020,'freq_range',[4 44],'freqs',(4:44)'};

fprintf('  Running ERP-SPRiNT (contra)...\n');
res_co = n1f.run_erp_sprint(sig_contra, t_long', fs, sprint_args{:});
fprintf('  Running ERP-SPRiNT (ipsi)...\n');
res_ip = n1f.run_erp_sprint(sig_ipsi,   t_long', fs, sprint_args{:});
fprintf('  Running ERP-SPRiNT (control)...\n');
res_ct = n1f.run_erp_sprint(sig_ctrl,   t_long', fs, sprint_args{:});

%% ── Crop to display window ───────────────────────────────────────────────
for rname = {'res_co','res_ip','res_ct'}
    r  = eval(rname{1});
    dm = r.times >= t_ep(1) & r.times <= t_ep(end);
    r.times     = r.times(dm);
    r.exponent  = r.exponent(dm);
    r.r_squared = r.r_squared(dm);
    r.full_fit  = r.full_fit(:,dm);
    r.aperiodic = r.aperiodic(:,dm);
    r.tfr_avg   = r.tfr_avg(:,dm);
    r.peak_cf   = r.peak_cf(dm);
    r.peak_pw   = r.peak_pw(dm);
    eval([rname{1} ' = r;']);
end

ep_idx   = t_long >= t_ep(1) & t_long <= t_ep(end);
t_ms     = t_long(ep_idx)' * 1000;
t_bin_ms = res_co.times * 1000;

avg_co = mean(sig_contra(:,ep_idx), 1);
avg_ip = mean(sig_ipsi(:,ep_idx),   1);
avg_ct = mean(sig_ctrl(:,ep_idx),   1);

%% ── specparam-modelled spectrum, dB re baseline ──────────────────────────
% full_fit is log10 modelled power (freq x time). dB re baseline =
% 10*(full_fit - mean_baseline_per_freq). Rows without usable baseline or
% inside the edge-mask stay NaN (render white).
bl_t = res_co.times < 0;
model_db = @(res) local_model_db(res, bl_t);
db_co = model_db(res_co);
db_ip = model_db(res_ip);

%% ── Figure (2x2): A traces, B exponent, C contra TFR, D ipsi TFR ─────────
fig = figure('Units','centimeters','Position',[0 0 18 13],'Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

shade = @(ax) fill(ax, ...
    [erd_start erd_end erd_end erd_start]*1000, ...
    [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)], ...
    [.95 .61 .07],'FaceAlpha',.08,'EdgeColor','none','HandleVisibility','off');

% A: lateralisation ERP (trial-averaged signals)
ax_A = nexttile; hold(ax_A,'on');
hA_co = plot(ax_A,t_ms,avg_co,'Color',[cc 0.75],'LineWidth',1.8, ...
    'LineStyle','-','DisplayName','Contra: ERD (decrease)');
hA_ip = plot(ax_A,t_ms,avg_ip,'Color',[ic 0.75],'LineWidth',1.8, ...
    'LineStyle','-','DisplayName','Ipsi: ERS (increase)');
xline(ax_A,0,'--','Color',txt_col,'LineWidth',.8,'HandleVisibility','off');
axis(ax_A,'tight'); shade(ax_A);
ylabel(ax_A,'Amplitude'); xlabel(ax_A,'Time (ms)');
title(ax_A,'Trial-averaged signal','FontWeight','bold','Color',txt_col);
lgd_A = legend(ax_A,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',7);
lgd_A.TextColor = 'k';     % below the axes: the traces fill the panel, and the labels are long
set_ax(ax_A,txt_col); plabel(ax_A,'A',txt_col);
xlim(ax_A,[t_ms(1) t_ms(end)]);

% B: aperiodic exponent (contra vs ipsi)
ax_B = nexttile; hold(ax_B,'on');
plot(ax_B,t_bin_ms,res_co.exponent,'Color',cc,'LineWidth',1.3,'DisplayName','Contra');
plot(ax_B,t_bin_ms,res_ip.exponent,'Color',ic,'LineWidth',1.3,'DisplayName','Ipsi');
xline(ax_B,0,'--','Color',txt_col,'LineWidth',.8,'HandleVisibility','off');
axis(ax_B,'tight'); shade(ax_B);
ylabel(ax_B,'Exponent'); xlabel(ax_B,'Time (ms)');
title(ax_B,'Aperiodic exponent','FontWeight','bold','Color',txt_col);
lgd_B = legend(ax_B,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',7);
lgd_B.TextColor = 'k';     % same placement as A, so A and B keep equal heights
set_ax(ax_B,txt_col); plabel(ax_B,'B',txt_col);
xlim(ax_B,[t_bin_ms(1) t_bin_ms(end)]);

% Shared colour scale for the two TFR panels
fr_msk = res_co.freqs >= 4 & res_co.freqs <= 44;
all_db = [db_co(fr_msk,:) db_ip(fr_msk,:)];
clim_v = max(abs(prctile(all_db(:),[0.5 99.5])),[],'omitnan');
if isempty(clim_v) || ~isfinite(clim_v) || clim_v==0, clim_v = 6; end
yt_hz = [4 8 12 20 30 44];

% C: contra (ERD) — specparam-modelled spectrum, dB re baseline
ax_C = nexttile;
hC = imagesc(ax_C, t_bin_ms, log2(res_co.freqs(fr_msk)), db_co(fr_msk,:));
set(hC,'AlphaData', ~isnan(db_co(fr_msk,:)));
set(ax_C,'YDir','normal','Color','w');
colormap(ax_C,n1f.redblue(256)); clim(ax_C,[-clim_v clim_v]);
set(ax_C,'YTick',log2(yt_hz),'YTickLabel',arrayfun(@num2str,yt_hz,'UniformOutput',false));
xline(ax_C,0,'--k','LineWidth',.8,'HandleVisibility','off');
ylabel(ax_C,'Frequency (Hz)'); xlabel(ax_C,'Time (ms)');
hold(ax_C,'on');
plot(ax_C, t_bin_ms, log2(res_co.peak_cf), 'k-', 'LineWidth', 1.2, 'DisplayName', 'fitted peak centre');
plot(ax_C, t_bin_ms([1 end]), log2([f0 f0]), 'k:', 'LineWidth', .8, 'DisplayName', sprintf('carrier (%g Hz)', f0));
legend(ax_C,'Location','southwest','FontSize',6.5,'Box','off','TextColor','k');
title(ax_C,'Contra (ERD): specparam model','FontWeight','bold','Color',txt_col);
set_ax(ax_C,txt_col); plabel(ax_C,'C',txt_col);
xlim(ax_C,[t_bin_ms(1) t_bin_ms(end)]);

% D: ipsi (ERS) — specparam-modelled spectrum, dB re baseline
ax_D = nexttile;
hD = imagesc(ax_D, t_bin_ms, log2(res_ip.freqs(fr_msk)), db_ip(fr_msk,:));
set(hD,'AlphaData', ~isnan(db_ip(fr_msk,:)));
set(ax_D,'YDir','normal','Color','w');
colormap(ax_D,n1f.redblue(256)); clim(ax_D,[-clim_v clim_v]);
set(ax_D,'YTick',log2(yt_hz),'YTickLabel',arrayfun(@num2str,yt_hz,'UniformOutput',false));
xline(ax_D,0,'--k','LineWidth',.8,'HandleVisibility','off');
ylabel(ax_D,'Frequency (Hz)'); xlabel(ax_D,'Time (ms)');
hold(ax_D,'on');
plot(ax_D, t_bin_ms, log2(res_ip.peak_cf), 'k-', 'LineWidth', 1.2, 'DisplayName', 'fitted peak centre');
plot(ax_D, t_bin_ms([1 end]), log2([f0 f0]), 'k:', 'LineWidth', .8, 'DisplayName', sprintf('carrier (%g Hz)', f0));
legend(ax_D,'Location','southwest','FontSize',6.5,'Box','off','TextColor','k');
title(ax_D,'Ipsi (ERS): specparam model','FontWeight','bold','Color',txt_col);
cb = colorbar(ax_D); cb.FontSize=7; cb.Color='k';
cb.Label.String='Model power (dB re baseline)'; cb.Label.Color='k';
set_ax(ax_D,txt_col); plabel(ax_D,'D',txt_col);
xlim(ax_D,[t_bin_ms(1) t_bin_ms(end)]);

%% ── Caption statistics ───────────────────────────────────────────────────
bl       = res_co.times >= -0.60 & res_co.times < -0.10;
erd_win  = res_co.times >= erd_start & res_co.times <= erd_end;
peak_win = res_co.times >= (erd_start+(erd_end-erd_start)/3) & ...
           res_co.times <= (erd_end  -(erd_end-erd_start)/3);

fprintf('\n--- fig3_erp_sprint caption statistics ---\n');
for row = {'Contra', res_co; 'Ipsi', res_ip; 'Control', res_ct}'
    label = row{1}; rex = row{2};
    fprintf('%s exponent: baseline=%.2f  ERD=%.2f  delta=%.2f\n', ...
        label, nanmean(rex.exponent(bl)), nanmean(rex.exponent(erd_win)), ...
        nanmean(rex.exponent(erd_win))-nanmean(rex.exponent(bl)));
    fprintf('%s R2:       baseline=%.2f  ERD=%.2f\n', ...
        label, nanmean(rex.r_squared(bl)), nanmean(rex.r_squared(erd_win)));
    fprintf('%s peak pw:  baseline=%.3f  peak=%.3f  delta=%.3f\n', ...
        label, nanmean(rex.peak_pw(bl)), nanmean(rex.peak_pw(peak_win)), ...
        nanmean(rex.peak_pw(peak_win))-nanmean(rex.peak_pw(bl)));
    fprintf('%s peak cf:  baseline=%.2f  peak=%.2f  delta=%.2f\n\n', ...
        label, nanmean(rex.peak_cf(bl)), nanmean(rex.peak_cf(peak_win)), ...
        nanmean(rex.peak_cf(peak_win))-nanmean(rex.peak_cf(bl)));
end
fprintf('-------------------------------\n\n');

save_fig(fig, cfg, 'fig3_erp_sprint');
end

function db = local_model_db(res, bl_t)
% specparam-modelled spectrum in dB relative to per-frequency baseline mean.
db = nan(size(res.full_fit));
for fi = 1:numel(res.freqs)
    row  = res.full_fit(fi,:);              % log10 modelled power
    bl_v = row(bl_t & ~isnan(row));
    if numel(bl_v) < 3, continue; end
    db(fi,:) = 10 * (row - mean(bl_v));     % dB re baseline (log10 -> dB)
end
end

function set_ax(ax,txt_col)
set(ax,'Box','off','TickDir','out','FontSize',8,'LineWidth',.8, ...
    'XColor',txt_col,'YColor',txt_col,'Color','w');
end
function plabel(ax,letter,txt_col)
text(ax,-0.13,1.05,letter,'Units','normalized','FontSize',11, ...   % clear of the top y-tick labels
     'FontWeight','bold','Color',txt_col);
end
function save_fig(fig,cfg,name)
fname = fullfile(cfg.out_dir,[name '.png']);
exportgraphics(fig,fname,'Resolution',300);
fprintf('  Saved: %s\n',fname); %close(fig);
end