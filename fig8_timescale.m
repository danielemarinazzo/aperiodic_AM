function fig8_timescale()
% Two ends of the timescale axis (Reviewer 1).
%
% Top (A-D): "A slow AM process with a period much longer than the recording
% window will appear stationary within that window."
%   10 Hz oscillation, sinusoidal AM with period 40 s (depth 0.8), white noise
%   0.02. A: 2 s window centred on an envelope maximum; B: the 120 s record
%   (shaded: the 2 s window). C: spectrum of the 2 s window vs an unmodulated
%   oscillation; D: spectrum of the 120 s record around 10 Hz -- carrier plus
%   sidebands at +/- 1/40 Hz.
% Bottom (E-G): "Individual action potentials are the canonical 1/f-spectrum
% signal when averaged across many cells and time, yet each spike is itself a
% deterministic event."
%   E: one brief transient (rise 1 ms, decay 10 ms). F: the transient summed
%   over 500 Poisson cells (5 Hz each), independent vs rate-synchronised at
%   10 Hz (kappa = 0.8). G: spectra for kappa = 0, 0.2, 0.8 and the single
%   transient's |K(f)|^2 (scaled). Independent summation keeps the single-
%   transient shape -- a decaying spectrum with a corner at 1/(2*pi*tau_d),
%   not a power law (Campbell); synchrony adds a peak. fs = 2000 Hz for E-G.

cfg = n1f.config();
fs  = cfg.fs;
ce  = [.84 .19 .12]; cs = [.99 .55 .24];

%% ── Top: slow AM ─────────────────────────────────────────────────────────
T  = 120; Pm = 40; depth = 0.8; f0 = 10;
t  = (0:round(T*fs)-1)' / fs;
rng(4101);
env = 1 + depth*sin(2*pi*t/Pm);
x   = env .* cos(2*pi*f0*t) + 0.02*randn(numel(t),1);
x0  = cos(2*pi*f0*t)        + 0.02*randn(numel(t),1);      % unmodulated
tc  = Pm/4 + Pm;                                             % envelope maximum (50 s)
w   = t >= tc - 1 & t < tc + 1;
[fw,  pw ] = n1f.compute_psd(x(w),  fs, sum(w));
[fw0, pw0] = n1f.compute_psd(x0(w), fs, sum(w));
[fL,  pL ] = n1f.compute_psd(x,     fs, numel(x));

%% ── Bottom: brief transients ─────────────────────────────────────────────
fs2 = 2000; N2 = cfg.duration*fs2; t2 = (0:N2-1)'/fs2;
tau_r = 0.001; tau_d = 0.010;
tk = (0:1/fs2:0.1)';
K  = exp(-tk/tau_d) - exp(-tk/tau_r); K = K / max(K);
n_cells = 500; lam0 = 5;
kappas = [0 0.2 0.8];
cK = {[.30 .30 .30], [.99 .55 .24], [.89 .10 .11]};
pop = cell(1, numel(kappas));
for q = 1:numel(kappas)
    rng(77);
    pop{q} = pop_signal(n_cells, lam0, kappas(q), t2, fs2, K);
end
expK = nan(20, numel(kappas)); npkK = expK;              % does synchrony move the exponent?
for q = 1:numel(kappas)
    for k = 1:20
        rng(3000 + k);
        [fq, pq] = n1f.compute_psd(pop_signal(n_cells, lam0, kappas(q), t2, fs2, K), fs2);
        sp = n1f.run_specparam(fq, pq, [1 45]);
        expK(k,q) = sp.exponent; npkK(k,q) = sp.n_peaks;
    end
end

%% ── Figure ───────────────────────────────────────────────────────────────
fig = figure('Units','centimeters','Position',[0 0 19 18],'Color','w');
lm = 0.08; rm = 0.03; tm = 0.05; bm = 0.07; hg = 0.08; vg = 0.09;
W  = 1 - lm - rm;
rh = (1 - tm - bm - 2*vg) / 3;
rb = @(r) 1 - tm - r*rh - (r-1)*vg;
wL = 0.34*W; wR = W - wL - hg;

% A: 2 s window
ax = axes(fig,'Position',[lm, rb(1), wL, rh]); hold(ax,'on');
plot(ax, t(w) - tc, x(w), 'k', 'LineWidth', .8);
xlim(ax,[-1 1]); set(ax,'YTick',[]); xlabel(ax,'Time (s)');
title(ax,'2 s window','FontSize',8,'FontWeight','bold'); set_ax(ax); plabel(ax,'A');

% B: 120 s record
ax = axes(fig,'Position',[lm + wL + hg, rb(1), wR, rh]); hold(ax,'on');
yl = [-1 1]*1.40*max(env);                       % headroom for the window label
% 1200 carrier cycles cannot be resolved at this width: draw the carrier in light grey,
% the envelope (+/-) on top, and the 2 s window last so that it stays visible.
plot(ax, t, x, 'Color', [.78 .78 .78], 'LineWidth', .25);
plot(ax, t,  env, 'Color', ce, 'LineWidth', 1.2);
plot(ax, t, -env, 'Color', ce, 'LineWidth', 1.2);
fill(ax, tc + [-1 1 1 -1], yl([1 1 2 2]), cs, 'FaceAlpha', .55, 'EdgeColor', [.85 .40 .05], 'LineWidth', 1);
text(ax, tc + 2, yl(2), '\leftarrow 2 s window (A)', 'FontSize', 7, 'Color', [.70 .30 .00], ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
xlim(ax,[0 T]); ylim(ax, yl); set(ax,'YTick',[]); xlabel(ax,'Time (s)');
title(ax,sprintf('%g s record, AM period %g s (grey: carrier; red: envelope)', T, Pm),'FontSize',8,'FontWeight','bold');
set_ax(ax); plabel(ax,'B');

% C: spectrum of the 2 s window vs unmodulated
ax = axes(fig,'Position',[lm, rb(2), wL, rh]); hold(ax,'on');
m = fw >= 1 & fw <= 30;
plot(ax, fw(m),  pw(m)  / max(pw(m)),  'Color', ce, 'LineWidth', 1.2, 'DisplayName', 'slowly modulated');
plot(ax, fw0(m), pw0(m) / max(pw0(m)), 'k--', 'LineWidth', .9, 'DisplayName', 'unmodulated');
set(ax,'YScale','log'); xlim(ax,[1 30]);
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power (normalised)');
legend(ax,'Location','northeast','FontSize',6.5,'Box','off','TextColor','k');
title(ax,'Spectrum of the 2 s window','FontSize',8,'FontWeight','bold'); set_ax(ax); plabel(ax,'C');

% D: spectrum of the 120 s record around the carrier
ax = axes(fig,'Position',[lm + wL + hg, rb(2), wR, rh]); hold(ax,'on');
m = fL >= f0 - 0.15 & fL <= f0 + 0.15;
plot(ax, fL(m), pL(m) / max(pL(m)), 'Color', ce, 'LineWidth', 1);
for s = [-1 1], xline(ax, f0 + s/Pm, ':', 'Color', [.4 .4 .4], 'LineWidth', .8); end
xlim(ax, f0 + [-0.15 0.15]); xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power (normalised)');
title(ax,sprintf('Spectrum of the %g s record: carrier, sidebands at \\pm1/%g Hz', T, Pm), ...
      'FontSize',8,'FontWeight','bold');
set_ax(ax); plabel(ax,'D');

% E: one transient
w3 = (W - 2*hg) / 3;
ax = axes(fig,'Position',[lm, rb(3), w3*0.7, rh]); hold(ax,'on');
plot(ax, tk*1000, K, 'k', 'LineWidth', 1.2);
xlim(ax,[0 60]); set(ax,'YTick',[]); xlabel(ax,'Time (ms)');
title(ax,'One transient','FontSize',8,'FontWeight','bold'); set_ax(ax); plabel(ax,'E');

% F: summed over cells, independent vs synchronised
xF = lm + w3*0.7 + hg; wF = w3*1.3;
ax = axes(fig,'Position',[xF, rb(3), wF, rh]); hold(ax,'on');
seg2 = t2 >= 1 & t2 < 2;
yI = pop{1} - mean(pop{1}); yS = pop{3} - mean(pop{3});
off = max(yS(seg2)) - min(yI(seg2));
plot(ax, t2(seg2), yI(seg2) + off, 'Color', cK{1}, 'LineWidth', .5, 'DisplayName', 'independent (\kappa = 0)');
plot(ax, t2(seg2), yS(seg2),       'Color', cK{3}, 'LineWidth', .5, 'DisplayName', 'synchronised (\kappa = 0.8)');
ylv = [min(yS(seg2)), max(yI(seg2)) + off]; r = diff(ylv);
ylim(ax, [ylv(1) - 0.03*r, ylv(2) + 0.40*r]);
legend(ax,'Location','north','Orientation','horizontal','FontSize',6.5,'Box','off','TextColor','k');
xlim(ax,[1 2]); set(ax,'YTick',[]); xlabel(ax,'Time (s)');
title(ax,sprintf('Summed over %d cells', n_cells),'FontSize',8,'FontWeight','bold');
set_ax(ax); plabel(ax,'F');

% G: spectra
ax = axes(fig,'Position',[xF + wF + hg, rb(3), W - (xF - lm) - wF - hg, rh]); hold(ax,'on');
yall = [];
for q = 1:numel(kappas)
    [fq, pq] = n1f.compute_psd(pop{q}, fs2); md = fq >= 1 & fq <= 200;
    plot(ax, fq(md), pq(md), 'Color', cK{q}, 'LineWidth', .9, 'DisplayName', sprintf('\\kappa = %g', kappas(q)));
    yall = [yall; pq(md)]; %#ok<AGROW>
    if q == 1, fI = fq(md); pI = pq(md); end
end
KF = 1 ./ ((1 + (2*pi*fI*tau_r).^2) .* (1 + (2*pi*fI*tau_d).^2));
sc = 0.6 * median(pI(fI <= 3)) / median(KF(fI <= 3));
plot(ax, fI, sc*KF, ':k', 'LineWidth', 1, 'DisplayName', 'one transient (scaled)');
set(ax,'XScale','log','YScale','log'); xlim(ax,[1 200]);
ylG = [10^floor(log10(min(yall))) 10^ceil(log10(max(yall)))]; ylim(ax, ylG);
hf = fill(ax,[1 45 45 1],[ylG(1) ylG(1) ylG(2) ylG(2)],[.5 .5 .5],'FaceAlpha',.08, ...
          'EdgeColor','none','HandleVisibility','off');
uistack(hf, 'bottom');
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power');
legend(ax,'Location','southwest','FontSize',6.5,'Box','off','TextColor','k');
title(ax,'Spectra (grey: 1-45 Hz)','FontSize',8,'FontWeight','bold');
set_ax(ax); plabel(ax,'G');

%% ── Caption statistics ───────────────────────────────────────────────────
fprintf('\n--- fig_timescale caption statistics ---\n');
ew = env(w);
fprintf('A-B: AM period %g s, depth %.1f; 2 s window centred at %g s (envelope maximum)\n', Pm, depth, tc);
fprintf('     envelope change within the 2 s window: %.2f%% of its mean\n', 100*(max(ew)-min(ew))/mean(ew));
pn = pw/max(pw); pn0 = pw0/max(pw0);
m = pn > 1e-3 & pn0 > 1e-3;                 % peak region only (the floor is independent noise)
fprintf('C:   2 s spectra, peak region (%.1f-%.1f Hz): max |log10 ratio| modulated/unmodulated = %.3f\n', ...
        min(fw(m)), max(fw(m)), max(abs(log10(pn(m) ./ pn0(m)))));
fprintf('D:   120 s record, resolution %.4f Hz; sidebands expected at %.3f and %.3f Hz\n', ...
        1/T, f0 - 1/Pm, f0 + 1/Pm);
fprintf('E-G: transient corner 1/(2*pi*tau_d) = %.0f Hz\n', 1/(2*pi*tau_d));
for q = 1:numel(kappas)
    fprintf('     kappa = %-4g  exponent %+.2f +/- %.2f  mean n_peaks %.1f (n = 20)\n', ...
            kappas(q), mean(expK(:,q)), std(expK(:,q)), mean(npkK(:,q)));
end
fprintf('---\n\n');
save_fig(fig, cfg, 'fig_timescale');
end

%% ── Local helpers ────────────────────────────────────────────────────────
function y = pop_signal(n_cells, lam0, kappa, t, fs, K)
% Summed population of Poisson cells sharing the rate lam0*(1 + kappa*cos(2*pi*10*t)),
% each spike replaced by the transient K. Summed counts are Poisson with the
% population rate, so cells need not be simulated individually.
lam = n_cells * lam0 * (1 + kappa*cos(2*pi*10*t)) / fs;
y = conv(poissrnd(lam), K);
y = y(1:numel(t));
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
text(ax,-0.12,1.06,letter,'Units','normalized','FontSize',12, ...
     'FontWeight','bold','VerticalAlignment','bottom','Color','k');
end

function save_fig(fig, cfg, name)
fname = fullfile(cfg.out_dir,[name '.png']);
exportgraphics(fig, fname, 'Resolution', 150);
fprintf('  Saved: %s\n', fname); close(fig);
end
