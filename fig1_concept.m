function fig1_concept()
% Two readings of the same exponent change.
%
% Row A (conventional reading): the broadband generator steepens (1/f exponent
%   -1 -> e_top); the oscillation is unchanged.
% Row B (reading proposed here): the broadband generator is fixed (exponent
%   -1); the oscillation changes from constant amplitude (stationary
%   sinusoid) to fluctuating amplitude and phase (narrowband process with a
%   Lorentzian spectrum, half-width B_dyn; n1f.sim_narrowband_osc 'lorentz').
% Both rows: 10 Hz oscillation (unit variance) + 1/f with variance rho x the
%   oscillatory variance + white noise sigma 0.02. Condition 1 is the same
%   signal in both rows.
% Columns: generators (2 s; envelope dashed); summed signal (2 s); mean
%   Welch spectrum over NR realizations with the aperiodic fit drawn from the
%   mean offset and mean exponent of the NR specparam fits (legend: exponent
%   mean +/- SD). Row B: shaded = power of condition 2 above condition 1
%   outside +/-1.5 Hz of the carrier (power moved from the line into the
%   tails of the peak).
%
% e_top is set so the two rows give the same exponent change (Python
% validation, n = 30, same seeds: +0.074 vs +0.077, rank-sum p = 1; MATLAB's
% RNG differs, see printout). Lorentzian tails decay as 1/(f-f0)^2, so part
% of the moved power falls under the aperiodic fit even when the peak core
% stays narrow (Python, rho = 1: shift +0.07 for half-widths 0.25-2 Hz;
% unchanged for 10, 60, 300 s recordings).
% Previous version (rho = 0.01, Gaussian bandwidth 3 Hz, e_top = -1.3):
% archive/superseded_2026-10/fig_concept_rho001_gauss.m
%
% Note: n1f.sim_powerlaw calls rng(seed) internally; the 1/f components are
% generated first, then the stream is reseeded.

cfg = n1f.config();
fs  = cfg.fs;
N   = cfg.duration * fs;
t   = (0:N-1)' / fs;
NR  = 30;
f0  = 10;
rho = 1;
e_top = -1.07;
B_dyn = 1;              % Lorentzian half-width (Hz)
shape = 'lorentz';
c1 = [.30 .30 .30];
c2 = [.89 .10 .11];

%% ── Statistics over realizations ─────────────────────────────────────────
% params(k, condition, :) = [exponent offset peak_cf peak_pw peak_bw];
% conditions: 1 = condition 1 (both rows), 2 = A condition 2, 3 = B condition 2
params = nan(NR, 3, 5);
fq = [];
for k = 1:NR
    s = make_conditions(k, N, fs, f0, rho, e_top, B_dyn, shape);
    X = {s.x1, s.xA2, s.xB2};
    for c = 1:3
        [fq, pw] = n1f.compute_psd(X{c}, fs);
        if k == 1 && c == 1, Psum = zeros(numel(fq), 3); end
        Psum(:,c) = Psum(:,c) + pw;
        params(k,c,:) = fit_params(fq, pw, f0);
    end
end
Pm = Psum / NR;
ex = params(:,:,1);
of = params(:,:,2);

%% ── Example realization (time courses) ───────────────────────────────────
s = make_conditions(1, N, fs, f0, rho, e_top, B_dyn, shape);
seg = t >= 2 & t < 4;
m   = fq >= 1 & fq <= 45;
fa  = fq(m);

rows = struct( ...
    'x2',   {s.xA2, s.xB2}, ...
    'col2', {2, 3}, ...
    'osc2', {[], s.oB}, ...
    'bb2',  {s.b2, []}, ...
    'ttl',  {'A   Conventional reading: the broadband generator changed', ...
             'B   Proposed reading: broadband generator fixed, amplitude dynamics changed'});

%% ── Figure ───────────────────────────────────────────────────────────────
fig = figure('Units','centimeters','Position',[0 0 19 12.5],'Color','w');
lm = 0.05; rm = 0.02; tm = 0.09; bm = 0.09; vg = 0.15; hg = 0.06; gg = 0.085;
hb = (1 - tm - bm - vg) / 2;                     % height of one row block
yb = [bm + hb + vg, bm];                         % block bottoms
W  = 1 - lm - rm - 2*hg;
cw = W * [1 1 1.1] / 3.1;                        % column widths
cx = [lm, lm + cw(1) + hg, lm + cw(1) + cw(2) + 2*hg];
gh = (hb - gg) / 2;                              % generator axes height

for r = 1:2
    R = rows(r); j = R.col2;
    annotation(fig,'textbox',[0.01, yb(r)+hb+0.035, 0.98, 0.04], 'String', R.ttl, ...
        'FontSize',9,'FontWeight','bold','EdgeColor','none','Color','k', ...
        'VerticalAlignment','bottom','HorizontalAlignment','left');

    % col 1, top: oscillation generator
    ax = axes(fig,'Position',[cx(1), yb(r)+gh+gg, cw(1), gh]); hold(ax,'on');
    plot(ax, t(seg), s.o0(seg), 'Color', c1, 'LineWidth', .5);
    if isempty(R.osc2)
        ttl = 'oscillation: unchanged'; tc = [0 0 0];
    else
        e2 = abs(hilbert(R.osc2));                       % full signal, then cut
        plot(ax, t(seg), R.osc2(seg), 'Color', c2, 'LineWidth', .5);
        plot(ax, t(seg), e2(seg), '--', 'Color', c2, 'LineWidth', .8);
        ttl = 'oscillation: changed (fluctuating amplitude)'; tc = c2;
    end
    gen_ax(ax, ttl, tc, false);

    % col 1, bottom: broadband generator
    ax = axes(fig,'Position',[cx(1), yb(r), cw(1), gh]); hold(ax,'on');
    plot(ax, t(seg), s.b1(seg), 'Color', c1, 'LineWidth', .5);
    if isempty(R.bb2)
        ttl = 'broadband (1/f): unchanged'; tc = [0 0 0];
    else
        plot(ax, t(seg), R.bb2(seg), 'Color', c2, 'LineWidth', .5);
        ttl = sprintf('broadband (1/f): exponent -1 \\rightarrow %g', e_top); tc = c2;
    end
    gen_ax(ax, ttl, tc, true);

    % col 2: summed signal
    ax = axes(fig,'Position',[cx(2), yb(r), cw(2), hb]); hold(ax,'on');
    plot(ax, t(seg), s.x1(seg) + 5, 'Color', c1, 'LineWidth', .5);
    plot(ax, t(seg), R.x2(seg) - 5, 'Color', c2, 'LineWidth', .5);
    text(ax, 2.0, 9.6,  'condition 1', 'FontSize', 7, 'Color', c1);
    text(ax, 2.0, -0.6, 'condition 2', 'FontSize', 7, 'Color', c2);
    xlim(ax,[2 4]); ylim(ax,[-10.5 11]); xlabel(ax,'Time (s)');
    title(ax,'signal = oscillation + broadband','FontSize',7.5,'FontWeight','normal');
    set_ax(ax); ax.YColor = 'none'; ax.YTick = [];
    ax.TitleHorizontalAlignment = 'left';

    % col 3: mean spectra + aperiodic fit from mean parameters
    ax = axes(fig,'Position',[cx(3), yb(r), cw(3), hb]); hold(ax,'on');
    P1 = Pm(m,1); Pj = Pm(m,j);
    if r == 2
        tail = Pj > P1 & abs(fa - f0) > 1.5;
        up = P1; up(tail) = Pj(tail);
        patch(ax, [fa; flipud(fa)], [up; flipud(P1)], c2, ...
              'FaceAlpha', .15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    plot(ax, fa, P1, 'Color', c1, 'LineWidth', .9, ...
         'DisplayName', sprintf('cond. 1: exp %+.2f \\pm %.2f', mean(ex(:,1)), std(ex(:,1))));
    plot(ax, fa, Pj, 'Color', c2, 'LineWidth', .9, ...
         'DisplayName', sprintf('cond. 2: exp %+.2f \\pm %.2f', mean(ex(:,j)), std(ex(:,j))));
    plot(ax, fa, 10.^(mean(of(:,1)) - mean(ex(:,1))*log10(fa)), '--', 'Color', c1, 'LineWidth', .9, 'HandleVisibility','off');
    plot(ax, fa, 10.^(mean(of(:,j)) - mean(ex(:,j))*log10(fa)), '--', 'Color', c2, 'LineWidth', .9, 'HandleVisibility','off');
    set(ax,'XScale','log','YScale','log'); xlim(ax,[1 45]); ylim(ax,[1e-3 30]);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power');
    legend(ax,'Location','northeast','FontSize',6.5,'Box','off','TextColor','k');
    title(ax, sprintf('Mean spectrum and aperiodic fit (n = %d)', NR), 'FontSize', 7.5, 'FontWeight', 'normal');
    if r == 1
        txt = 'steeper broadband component: steeper fit';
    else
        txt = {'power moved from the line into the tails', 'falls under the aperiodic fit: steeper fit'};
    end
    text(ax, 1.1, 1.5e-3, txt, 'FontSize', 6.5, 'Color', c2, 'VerticalAlignment', 'bottom');
    set_ax(ax);
    ax.TitleHorizontalAlignment = 'left';
end

%% ── Caption statistics ───────────────────────────────────────────────────
fprintf('\n--- fig_concept caption statistics (n = %d) ---\n', NR);
fprintf('rho = %g, e_top = %g, oscillation in B: %s, half-width %g Hz\n', rho, e_top, shape, B_dyn);
fprintf('Row A: only the broadband generator changed. Row B: only the oscillation changed.\n');
pn = {'exponent','offset','peak cf (Hz)','peak power','peak bw (Hz)'};
fprintf('%-13s %16s | %-30s | %-30s\n', 'parameter', 'condition 1', 'A cond. 2 (delta, signed-rank p)', 'B cond. 2 (delta, signed-rank p)');
for q = 1:5
    a1 = params(:,1,q); a2 = params(:,2,q); b2 = params(:,3,q);
    fprintf('%-13s %+8.3f +/- %.3f | %+8.3f (%+.3f, p=%.2g) | %+8.3f (%+.3f, p=%.2g)\n', pn{q}, ...
        mean(a1,'omitnan'), std(a1,'omitnan'), ...
        mean(a2,'omitnan'), mean(a2-a1,'omitnan'), signrank_safe(a2-a1), ...
        mean(b2,'omitnan'), mean(b2-a1,'omitnan'), signrank_safe(b2-a1));
end
dA = ex(:,2) - ex(:,1); dB = ex(:,3) - ex(:,1);
fprintf('exponent change A vs B: %+.3f vs %+.3f, ranksum p = %.2g\n', mean(dA), mean(dB), ranksum(dA, dB));
fprintf('(peak = specparam peak nearest %g Hz; bw = specparam bandwidth parameter)\n---\n\n', f0);

save_fig(fig, cfg, 'fig_concept');
end

%% ── Local helpers ────────────────────────────────────────────────────────
function s = make_conditions(k, N, fs, f0, rho, e_top, B_dyn, shape)
% Same white-noise seed for both 1/f components, so only their spectral
% shaping differs; then reseed for the oscillations and white noise.
seed = 700 + k;
pl1 = n1f.sim_powerlaw(N, fs, -1,    seed + 200); pl1 = pl1 / std(pl1);
pl2 = n1f.sim_powerlaw(N, fs, e_top, seed + 200); pl2 = pl2 / std(pl2);
rng(seed);
wn   = 0.02*randn(N,1);
s.o0 = n1f.sim_narrowband_osc(N, fs, f0, 0);
s.oB = n1f.sim_narrowband_osc(N, fs, f0, B_dyn, shape);
s.b1 = sqrt(rho)*pl1;
s.b2 = sqrt(rho)*pl2;
s.x1  = s.o0 + s.b1 + wn;
s.xA2 = s.o0 + s.b2 + wn;
s.xB2 = s.oB + s.b1 + wn;
end

function p = fit_params(fq, pw, f0)
% [exponent offset peak_cf peak_pw peak_bw]; peak = the one nearest f0.
sp = n1f.run_specparam(fq, pw, [1 45]);
p = [sp.exponent sp.offset NaN NaN NaN];
if sp.n_peaks > 0
    [~, kk] = min(abs([sp.peaks.cf] - f0));
    p(3:5) = [sp.peaks(kk).cf sp.peaks(kk).pw sp.peaks(kk).bw];
end
end

function gen_ax(ax, ttl, tc, show_x)
xlim(ax,[2 4]); ylim(ax,[-4 4]);
title(ax, ttl, 'FontSize', 7, 'FontWeight', 'normal');
set_ax(ax);
ax.TitleHorizontalAlignment = 'left';
ax.Title.Color = tc;
ax.YColor = 'none'; ax.YTick = [];
if show_x, xlabel(ax,'Time (s)'); else, ax.XTickLabel = []; end
end

function p = signrank_safe(d)
d = d(~isnan(d));
if isempty(d) || all(d == 0), p = NaN; else, p = signrank(d); end
end

function set_ax(ax)
set(ax,'Box','off','TickDir','out','FontSize',8,'LineWidth',.7, ...
    'Color','w','XColor','k','YColor','k');
ax.Title.Color = 'k';
end

function save_fig(fig, cfg, name)
fname = fullfile(cfg.out_dir,[name '.png']);
exportgraphics(fig, fname, 'Resolution', 150);
fprintf('  Saved: %s\n', fname); close(fig);
end
