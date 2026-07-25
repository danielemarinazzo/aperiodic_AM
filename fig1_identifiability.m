function fig1_identifiability()
% Figure 1 — The spectral identifiability problem.
% Two signals produced by different generative processes yield similar-
% looking power spectra. Panel A: a pure oscillation with amplitude
% modulation, no aperiodic component. Panel B: an oscillation superimposed
% on a stationary 1/f broadband process. Both are summarised by specparam
% as a peak plus an aperiodic exponent. The point is that the spectrum and
% the specparam decomposition do not identify which generative process is
% present; time-domain features are needed to recover the distinction.

cfg = n1f.config();
fs  = cfg.fs;
N   = cfg.duration * fs;
t   = (0:N-1)' / fs;
rng(cfg.seed);

% A: non-stationary AM (ERD-like piecewise suppression, Gaussian-smoothed)
envA = ones(N,1);
envA(3*fs : 7*fs) = 0.12;
envA = smoothdata(envA, 'gaussian', round(0.4*fs));
sigA_raw = envA .* cos(2*pi*10*t);
scaleA   = 1 / std(sigA_raw);
sigA = sigA_raw * scaleA;
envA = envA * scaleA;          % keep dashed envelope overlay tracking sigA

% B: oscillation + stationary 1/f broadband
nf   = n1f.sim_powerlaw(N, fs, -1, cfg.seed+1);
sigB = cos(2*pi*10*t) + nf/std(nf)*1.2;
sigB = sigB / std(sigB);

sigs   = {sigA,           sigB};
envs   = {envA,           []};
cols   = {[.90 .27 .21],  [.56 .27 .68]};
labels = {'A - Pure oscillation, non-stationary amplitude', ...
          'B - Oscillation + stationary 1/f broadband'};

fig = figure('Units','centimeters','Position',[0 0 18 11],'Color','w');
t5  = t(1:5*fs);

for i = 1:2
    sig = sigs{i}; c = cols{i};

    % Left: time domain
    ax = subplot(2,2,2*i-1);
    plot(ax, t5, sig(1:5*fs), 'Color',c, 'LineWidth',0.7);
    if ~isempty(envs{i})
        hold(ax,'on');
        plot(ax, t5,  envs{i}(1:5*fs), 'k--', 'LineWidth',1.5);
        plot(ax, t5, -envs{i}(1:5*fs), 'k--', 'LineWidth',1.5);
    end
    ylabel(ax,'Amplitude','FontSize',8);
    title(ax, labels{i}(4:end), 'FontSize',8,'FontWeight','bold', ...
          'HorizontalAlignment','left','Color','k');
    if i == 2, xlabel(ax,'Time (s)','FontSize',8); end
    set_ax(ax); plabel(ax, labels{i}(1));

    % Right: PSD + specparam
    ax2 = subplot(2,2,2*i);
    [freqs, psd] = n1f.compute_psd(sig, fs, fs*2);
    sp = n1f.run_specparam(freqs, psd, cfg.freq_range);
    semilogy(ax2, freqs, psd, 'Color',[.75 .75 .75], 'LineWidth',0.6);
    hold(ax2,'on');
    semilogy(ax2, sp.freqs, 10.^sp.modeled,   'Color',c, 'LineWidth',1.8);
    semilogy(ax2, sp.freqs, 10.^sp.aperiodic, 'Color',c, 'LineWidth',1.0, 'LineStyle','--');
    shade_bands(ax2);
    xlim(ax2,[1 45]);
    % Independent y-limits per panel based on data within fitting range
    p_in = psd(freqs >= 1 & freqs <= 45);
    ylim(ax2, [10^(floor(log10(min(p_in)))-0.3), 10^(ceil(log10(max(p_in)))+0.3)]);
    legend(ax2, sprintf('exp = %+.2f', sp.exponent), 'Location','northeast', ...
           'FontSize',7,'Box','off','TextColor','k');
    % title(ax2,'Power spectrum','FontSize',8,'FontWeight','bold', ...
    %       'HorizontalAlignment','left','Color','k');
    if i == 2, xlabel(ax2,'Frequency (Hz)','FontSize',8); end
    set_ax(ax2);
end

save_fig(fig, cfg, 'fig1_identifiability');
end


%% ── Local helpers ──────────────────────────────────────────────────────────
function shade_bands(ax)
hold(ax,'on');
defs = {[4 8],[.20 .60 .80]; [8 12],[.80 .20 .20]; [16 30],[.20 .70 .30]};
yl = get(ax,'YLim');
for i = 1:3
    x=defs{i,1}; c=defs{i,2};
    fill(ax,[x(1) x(2) x(2) x(1)],[yl(1) yl(1) yl(2) yl(2)], ...
         c,'FaceAlpha',.07,'EdgeColor','none','HandleVisibility','off');
end
end

function set_ax(ax)
set(ax,'Box','off','TickDir','out','FontSize',8,'LineWidth',.7, ...
    'Color','w','XColor','k','YColor','k');
ax.Title.Color = 'k';
end

function plabel(ax, letter)
text(ax,-0.11,1.05,letter,'Units','normalized','FontSize',12, ...
     'FontWeight','bold','VerticalAlignment','bottom','Color','k');
end

function save_fig(fig, cfg, name)
fname = fullfile(cfg.out_dir, [name '.png']);
exportgraphics(fig, fname, 'Resolution',150);
fprintf('  Saved: %s\n', fname);
close(fig);
end