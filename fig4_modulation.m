function fig4_modulation(f0, panel_C)
% Figure — what controls the exponent shift produced by oscillatory dynamics.
%
% Model (Rice): narrowband Gaussian oscillation at carrier f0 with spectral
% bandwidth B (n1f.sim_narrowband_osc) -> random amplitude and phase on a
% timescale ~1/B; B = 0 is a stationary sinusoid. Broadband floor: 1/f
% (n1f.sim_powerlaw, exponent -1, unit SD) scaled to power ratio rho relative
% to the oscillatory signal, plus white noise sigma 0.02 (as in the scenarios).
% PSD: n1f.compute_psd (Welch, 2 s Hann, 50%). Fit: n1f.run_specparam, fixed,
% [1 45] Hz, default settings.
%
%   A  time-varying oscillations at four bandwidths (exponent next to each label)
%   B  their spectra (white floor only)
%   C  exponent shift vs stationary oscillation, same floor, over B x rho
%   D  alpha interrupted by a beta burst: components
%   E  spectra, one realization (rho = 0 solid, rho = 1 dashed)
%   F  shift (burst - power-matched stationary) along the same rho axis as C
%
% fig_modulation()          carrier 10 Hz, panels A-F (main figure)
% fig_modulation(20)        carrier 20 Hz, panels A-C only (supplement)
% fig_modulation(f0, flag)  flag = true/false forces panels D-F on/off.
%   Default: D-F only for f0 = 10. They do not depend on f0 (always alpha 10 Hz /
%   beta 20 Hz), so other carriers would only duplicate them.
% (The brief-transient / synchrony panel moved to fig_timescale.)
%
% Note: n1f.sim_powerlaw calls rng(seed) internally. Each realization therefore
% generates the 1/f first and reseeds before drawing phases and noise.

cfg = n1f.config();
fs  = cfg.fs;
N   = cfg.duration * fs;
t   = (0:N-1)' / fs;
NR  = 30;
if nargin < 1 || isempty(f0), f0 = 10; end
if nargin < 2 || isempty(panel_C), panel_C = (f0 == 10); end

Bs   = [0.1 0.3 1 3 6 10];
rhos = [0 0.1 0.3 1 3];

%% ── B: shift map ─────────────────────────────────────────────────────────
M = nan(numel(rhos), numel(Bs)); S = M; P = M;
for i = 1:numel(rhos)
    for j = 1:numel(Bs)
        d = nan(NR,1);
        for k = 1:NR
            seed = 100000*i + 1000*j + k;
            pl = n1f.sim_powerlaw(N, fs, -1, seed); pl = pl / std(pl);
            rng(seed + 50000);
            fl = sqrt(rhos(i)) * pl + 0.02*randn(N,1);
            x0 = n1f.sim_narrowband_osc(N, fs, f0, 0);
            x1 = n1f.sim_narrowband_osc(N, fs, f0, Bs(j));
            d(k) = expo(x1 + fl, fs) - expo(x0 + fl, fs);
        end
        M(i,j) = mean(d); S(i,j) = std(d);
        if any(d ~= 0), P(i,j) = signrank(d); else, P(i,j) = 1; end
    end
    fprintf('  map row rho = %g done\n', rhos(i));
end
alpha_corr = 0.05 / numel(M);

if panel_C
    %% ── C: alpha interrupted by a beta burst ─────────────────────────────────
    eb = ramp_window(t, 4.5, 5.5, 0.05, fs);     % beta envelope, 50 ms ramps
    ea = 1 - 0.9*eb;                             % alpha suppressed to 10%
    rhosC  = [0 0.1 0.3 1 3];
    shiftC = nan(size(rhosC)); sdC = shiftC; pC = shiftC;
    ctrl   = nan(NR,2);                          % rho = 0: alpha only vs alpha+beta stationary
    for r = 1:numel(rhosC)
        d = nan(NR,1);
        for k = 1:NR
            seed = 900000 + k;
            pl = n1f.sim_powerlaw(N, fs, -1, seed); pl = pl / std(pl);
            rng(seed + 50000);
            [s_st, s_bu, s_al] = burst_signals(t, ea, eb, 2*pi*rand, 2*pi*rand);
            fl = sqrt(rhosC(r)*var(s_bu)) * pl + 0.02*randn(N,1);
            e_st = expo(s_st + fl, fs);
            d(k) = expo(s_bu + fl, fs) - e_st;
            if rhosC(r) == 0, ctrl(k,:) = [expo(s_al + fl, fs), e_st]; end
        end
        shiftC(r) = mean(d); sdC(r) = std(d); pC(r) = signrank(d);
    end
end

%% ── Figure ───────────────────────────────────────────────────────────────
n_rows = 2 + panel_C;
row_cm = 4.56; vg_cm = 2.03; tm_cm = 1.08; bm_cm = 1.62;   % 4 rows -> 27 cm, as before
H   = tm_cm + bm_cm + n_rows*row_cm + (n_rows-1)*vg_cm;
fig = figure('Units','centimeters','Position',[0 0 19 H],'Color','w');
cols = {[.17 .50 .72], [.25 .67 .36], [.99 .55 .24], [.89 .10 .11]};
cst  = [.17 .50 .72]; cbu = [.89 .10 .11];

lm = 0.08; rm = 0.03; hg = 0.07;
tm = tm_cm/H; bm = bm_cm/H; vg = vg_cm/H; rh = row_cm/H;
rb = @(r) 1 - tm - r*rh - (r-1)*vg;
W  = 1 - lm - rm;

% A1: traces
ax = axes(fig,'Position',[lm, rb(1), 0.60*W, rh]); hold(ax,'on');
seg = t >= 2 & t < 4;
lblA = {'B = 0.1 Hz (slow)','B = 1 Hz','B = 3 Hz','B = 10 Hz (brief events)'};
BA   = [0.1 1 3 10];
expA = nan(1,4); specA = cell(1,4);
for n = 1:4
    rng(4200 + n);
    x = n1f.sim_narrowband_osc(N, fs, f0, BA(n)) + 0.02*randn(N,1);
    yoff = (4-n)*5;
    [fq, pw] = n1f.compute_psd(x, fs);
    sp = n1f.run_specparam(fq, pw, [1 45]);
    expA(n) = sp.exponent; specA{n} = [fq pw];
    plot(ax, t(seg), x(seg) + yoff, 'Color', cols{n}, 'LineWidth', .6);
    text(ax, 4.03, yoff, {lblA{n}, sprintf('exponent %+.2f', expA(n))}, 'FontSize', 7, 'Color', cols{n});   % exponent here: no legend on the spectra
end
xlim(ax,[2 4.75]); set(ax,'YTick',[]); xlabel(ax,'Time (s)');
title(ax,sprintf('%g Hz oscillation with random amplitude and phase, bandwidth B', f0), ...
      'FontSize',8,'FontWeight','bold');
set_ax(ax); plabel(ax,'A');

% A2: spectra
ax = axes(fig,'Position',[lm + 0.60*W + hg, rb(1), W - 0.60*W - hg, rh]); hold(ax,'on');
for n = 1:4
    fq = specA{n}(:,1); pw = specA{n}(:,2); m = fq >= 1 & fq <= 45;
    plot(ax, fq(m), pw(m), 'Color', cols{n}, 'LineWidth', 1);
end
set(ax,'XScale','log','YScale','log'); xlim(ax,[1 45]);
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power');
title(ax,'Spectra (white floor only)','FontSize',8,'FontWeight','bold');
set_ax(ax); plabel(ax,'B');

% B: map
ax = axes(fig,'Position',[lm, rb(2), 0.72*W, rh]);
imagesc(ax, M);
vmax = max(abs(M(:)));
colormap(ax, n1f.redblue(256)); clim(ax, [-vmax vmax]);
for i = 1:numel(rhos)
    for j = 1:numel(Bs)
        star = ''; if P(i,j) < alpha_corr, star = '*'; end
        tc = 'k'; if abs(M(i,j)) > 0.6*vmax, tc = 'w'; end
        text(ax, j, i, {sprintf('%+.2f%s', M(i,j), star), sprintf('%c%.2f', char(177), S(i,j))}, ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'FontSize',7,'Color',tc,'Interpreter','none');
    end
end
set(ax,'XTick',1:numel(Bs),'XTickLabel',arrayfun(@num2str,Bs,'UniformOutput',false), ...
       'YTick',1:numel(rhos),'YTickLabel',arrayfun(@num2str,rhos,'UniformOutput',false));
xlabel(ax,'Bandwidth B of amplitude/phase fluctuations (Hz)  (faster \rightarrow)');
ylabel(ax,{'Broadband (1/f) /','oscillatory power'});   % two lines: shorter than the axis, clears the panel letter
title(ax,sprintf('Exponent shift vs stationary oscillation, same floor (mean \\pm SD, n = %d; * p < 0.05 Bonferroni)', NR), ...
      'FontSize',8,'FontWeight','bold');
cb = colorbar(ax); cb.FontSize = 7; cb.Color = 'k'; cb.Label.String = '\Delta exponent'; cb.Label.Color = 'k';
set_ax(ax); plabel(ax,'C');

if panel_C
    % C1: components
    cw = (W - 2*hg) / 3;
    ax = axes(fig,'Position',[lm, rb(3), cw, rh]); hold(ax,'on');
    rng(900001 + 50000); pa = 2*pi*rand; pb = 2*pi*rand;
    seg = t >= 4.1 & t < 5.9;
    fill(ax,[4.5 5.5 5.5 4.5],[-2.6 -2.6 3.3 3.3],[.99 .55 .24], ...
         'FaceAlpha',.12,'EdgeColor','none','HandleVisibility','off');
    plot(ax, t(seg), ea(seg).*cos(2*pi*10*t(seg)+pa) + 1.3, 'Color', cst, 'LineWidth', .7);
    plot(ax, t(seg), eb(seg).*cos(2*pi*20*t(seg)+pb) - 1.3, 'Color', cbu, 'LineWidth', .7);
    plot(ax, t(seg), ea(seg) + 1.3, ':k', 'LineWidth', .8);
    plot(ax, t(seg), eb(seg) - 1.3, ':k', 'LineWidth', .8);
    text(ax, 4.12, 2.85, 'alpha (10 Hz)', 'FontSize', 7, 'Color', cst);   % above the trace (top at 2.3)
    text(ax, 4.12, 0.0,  'beta (20 Hz)',  'FontSize', 7, 'Color', cbu);
    xlim(ax,[4.1 5.9]); ylim(ax,[-2.6 3.3]); set(ax,'YTick',[]); xlabel(ax,'Time (s)');
    title(ax,'Alpha interrupted by a beta burst','FontSize',8,'FontWeight','bold');
    set_ax(ax); plabel(ax,'D');

    % C2: spectra, one realization, rho = 0 and 1
    ax = axes(fig,'Position',[lm + cw + hg, rb(3), cw, rh]); hold(ax,'on');
    for rho = [0 1]
        pl = n1f.sim_powerlaw(N, fs, -1, 900001); pl = pl / std(pl);
        rng(900001 + 50000);
        [s_st, s_bu] = burst_signals(t, ea, eb, 2*pi*rand, 2*pi*rand);
        fl = sqrt(rho*var(s_bu)) * pl + 0.02*randn(N,1);
        ls = '-'; if rho > 0, ls = '--'; end
        [fq, pw] = n1f.compute_psd(s_st + fl, fs); m = fq >= 1 & fq <= 45;
        plot(ax, fq(m), pw(m), 'Color', cst, 'LineStyle', ls, 'LineWidth', .9, ...
             'DisplayName', sprintf('stat., \\rho=%g', rho));
        [fq, pw] = n1f.compute_psd(s_bu + fl, fs);
        plot(ax, fq(m), pw(m), 'Color', cbu, 'LineStyle', ls, 'LineWidth', .9, ...
             'DisplayName', sprintf('burst, \\rho=%g', rho));
    end
    set(ax,'XScale','log','YScale','log'); xlim(ax,[1 45]);
    yl = ylim(ax); ylim(ax, [yl(1) yl(2)*10^1.8]);   % headroom for the legend above the 20 Hz peak
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power');
    lg = legend(ax,'Location','northeast','FontSize',6,'Box','off','TextColor','k','NumColumns',2);
    lg.ItemTokenSize = [10 6];                       % short line samples: keeps the legend inside the axes
    title(ax,'Spectra (one realization)','FontSize',8,'FontWeight','bold');
    set_ax(ax); plabel(ax,'E');

    % C3: shift along the rho axis
    ax = axes(fig,'Position',[lm + 2*(cw + hg), rb(3), cw, rh]); hold(ax,'on');
    errorbar(ax, 1:numel(rhosC), shiftC, sdC, 'o-', 'Color', cbu, 'MarkerFaceColor', cbu, ...
             'MarkerSize', 4, 'LineWidth', 1.2);
    yline(ax, 0, ':', 'Color', [.6 .6 .6]);
    set(ax,'XTick',1:numel(rhosC),'XTickLabel',arrayfun(@num2str,rhosC,'UniformOutput',false));
    xlim(ax,[0.5 numel(rhosC)+0.5]);
    xlabel(ax,'Broadband (1/f) / oscillatory power'); ylabel(ax,'\Delta exponent (burst - stationary)');
    title(ax,sprintf('Shift vs broadband level (n = %d)', NR),'FontSize',8,'FontWeight','bold');
    set_ax(ax); plabel(ax,'F');
end

%% ── Caption statistics ───────────────────────────────────────────────────
fprintf('\n--- fig_modulation caption statistics (carrier %g Hz) ---\n', f0);
fprintf('C map, mean shift (rows rho = %s; cols B = %s Hz)\n', mat2str(rhos), mat2str(Bs));
disp(round(M,2)); fprintf('SD\n'); disp(round(S,2));
fprintf('Bonferroni-significant cells (signed-rank, alpha %.2g):\n', alpha_corr); disp(P < alpha_corr);
for n = 1:4, fprintf('A  B = %-4g Hz  exponent %+.2f\n', BA(n), expA(n)); end
if panel_C
    for r = 1:numel(rhosC)
        fprintf('C  rho = %-4g  shift %+.2f +/- %.2f  p = %.2g\n', rhosC(r), shiftC(r), sdC(r), pC(r));
    end
    fprintf('C control, rho = 0: alpha only %+.2f, alpha+beta stationary %+.2f, ranksum p = %.2g\n', ...
            mean(ctrl(:,1)), mean(ctrl(:,2)), ranksum(ctrl(:,1), ctrl(:,2)));
end
fprintf('Note: with n = %d paired realizations, shifts of ~0.03 reach significance;\n', NR);
fprintf('      the magnitude, not the star, carries the information.\n---\n\n');

name = 'fig_modulation';
if f0 ~= 10, name = sprintf('%s_f%g', name, f0); end
if panel_C ~= (f0 == 10)                      % non-default panel choice
    if panel_C, name = [name '_C']; else, name = [name '_AB']; end
end
save_fig(fig, cfg, name);
end

%% ── Local helpers ────────────────────────────────────────────────────────
function e = expo(x, fs)
[fq, pw] = n1f.compute_psd(x, fs);
sp = n1f.run_specparam(fq, pw, [1 45]);
e  = sp.exponent;
end

function e = ramp_window(t, t_on, t_off, ramp_s, fs)
% Boolean window [t_on, t_off) with raised-cosine transitions of ramp_s seconds.
on = t >= t_on & t < t_off;
e  = double(on);
L  = floor(ramp_s*fs);
k  = 0.5*(1 - cos(linspace(0, pi, L)'));
for i = find(diff(on) ~= 0)'
    if on(i+1), seg = k; else, seg = flipud(k); end
    a = i + 1 - floor(L/2);
    e(a:a+L-1) = seg;
end
end

function [s_st, s_bu, s_al] = burst_signals(t, ea, eb, pa, pb)
% Burst case, power-matched stationary alpha+beta, and stationary alpha alone.
s_bu = ea.*cos(2*pi*10*t + pa) + eb.*cos(2*pi*20*t + pb);
s_st = sqrt(mean(ea.^2))*cos(2*pi*10*t + pa) + sqrt(mean(eb.^2))*cos(2*pi*20*t + pb);
s_al = cos(2*pi*10*t + pa);
end

function set_ax(ax)
set(ax,'Box','off','TickDir','out','FontSize',8,'LineWidth',.7, ...
    'Color','w','XColor','k','YColor','k');
ax.Title.Color = 'k';
ax.Title.Units = 'normalized';
ax.Title.HorizontalAlignment = 'left';
ax.Title.Position(1) = 0;
end

function plabel(ax, letter)
% Letter at a fixed distance (in figure units) left of the axes, whatever its width.
dx = 0.055 / ax.Position(3);
text(ax,-dx,1.06,letter,'Units','normalized','FontSize',12, ...
     'FontWeight','bold','VerticalAlignment','bottom','Color','k');
end

function save_fig(fig, cfg, name)
fname = fullfile(cfg.out_dir,[name '.png']);
exportgraphics(fig, fname, 'Resolution', 150);
fprintf('  Saved: %s\n', fname); close(fig);
end
