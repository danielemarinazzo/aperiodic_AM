function fig3_toy()
% Figure — the Table 1 scenarios in feature space, and a toy discriminability
% test.
%
% A   The five illustrative scenarios (Table 1; n1f.make_scenario), NR
%     realizations each: specparam exponent vs envelope CV in the carrier band
%     (8-12 Hz for A-C, 16-30 Hz for D-E). Descriptive.
% B-D TOY 2x2 design at one 10 Hz carrier: amplitude dynamics (static
%     sinusoid vs narrowband process with a Lorentzian spectrum, half-width
%     1 Hz; n1f.sim_narrowband_osc 'lorentz') x broadband generator (1/f
%     exponent -1 vs -1.07, variance rho = 1 x oscillatory, + white noise
%     0.02). These are the conditions of fig_concept completed into a
%     factorial. Printout adds per-factor accuracy (dynamics: static vs
%     dynamic; broadband: -1 vs -1.07) read from the confusion matrices.
%     Python validation (n = 30): LOO accuracy exponent only 0.43, specparam
%     0.62, specparam + temporal 0.62, exponent + CV 0.68 (dynamics factor
%     1.00, broadband factor 0.68). Previous version (rho = 0.01, Gaussian
%     B = 3 Hz, 1/f -1.3): archive/superseded_2026-10/fig_toy_rho001_gauss.m
%     B  exponent vs envelope CV (8-12 Hz)
%     C  leave-one-out nearest-centroid accuracy (z-scored features) for four
%        feature sets
%     D  confusion matrices: exponent only vs exponent + envelope CV
%
% Class labels exist only because the signals are simulated, and only two
% generative factors vary. Nothing here licenses assigning recordings to
% generative classes, which in data are neither known nor exhaustive.

cfg = n1f.config();
fs  = cfg.fs;
N   = cfg.duration * fs;
NR  = 30;
ws  = warning('off','all');

%% ── A: Table 1 scenarios ─────────────────────────────────────────────────
skeys  = {'pure_alpha','am_alpha','alpha_plus_1f','burst_beta','burst_plus_1f'};
snames = {'A pure \alpha','B AM \alpha','C \alpha + 1/f','D bursty \beta','E bursty \beta + 1/f'};
bands  = {[8 12],[8 12],[8 12],[16 30],[16 30]};
scols  = {[.40 .40 .40],[.17 .50 .72],[.14 .55 .27],[.48 .00 .47],[.89 .10 .11]};
SA = nan(NR, numel(skeys), 2);                % [exponent, carrier-band CV]
for s = 1:numel(skeys)
    for k = 1:NR
        x = n1f.make_scenario(skeys{s}, 1000 + k);
        [fq, pw] = n1f.compute_psd(x, fs);
        sp = n1f.run_specparam(fq, pw, [1 45]);
        ev = n1f.envelope_cv(x, fs, bands{s});
        SA(k,s,:) = [sp.exponent ev.cv];
    end
end
fprintf('  Table 1 scenarios done\n');

%% ── B-D: toy 2x2 ─────────────────────────────────────────────────────────
rho = 1; e_top = -1.07; B_dyn = 1; shape = 'lorentz';     % = fig_concept
classes = struct('name', {'static, 1/f -1', sprintf('static, 1/f %g', e_top), 'dynamic, 1/f -1', sprintf('dynamic, 1/f %g', e_top)}, ...
                 'short', {'st, -1', sprintf('st, %g', e_top), 'dy, -1', sprintf('dy, %g', e_top)}, ...
                 'B', {0, 0, B_dyn, B_dyn}, 'e', {-1, e_top, -1, e_top});
ccols = {[.55 .55 .55], [.20 .20 .20], [.99 .55 .24], [.89 .10 .11]};
fnames = {'exponent','offset','peak cf','peak pw','peak bw','CV 8-12','SPRiNT SD'};
X = nan(NR*4, numel(fnames)); y = nan(NR*4, 1);
row = 0;
for c = 1:4
    for k = 1:NR
        pl = n1f.sim_powerlaw(N, fs, classes(c).e, 1900 + k); pl = pl / std(pl);
        rng(800 + k);
        x = n1f.sim_narrowband_osc(N, fs, 10, classes(c).B, shape) + sqrt(rho)*pl + 0.02*randn(N,1);
        row = row + 1;
        X(row,:) = toy_features(x, fs);
        y(row) = c;
    end
    fprintf('  toy class %d done\n', c);
end
warning(ws);
for j = 1:size(X,2)                           % rare missing peak: column median
    X(isnan(X(:,j)), j) = median(X(:,j), 'omitnan');
end
fsets = struct('name', {'exponent only', 'specparam (5)', 'specparam + temporal', 'exponent + CV'}, ...
               'cols', {1, 1:5, 1:7, [1 6]});
acc = nan(1, numel(fsets)); Cm = cell(1, numel(fsets));
for q = 1:numel(fsets)
    [acc(q), Cm{q}] = loo_centroid(X, y, fsets(q).cols);
end

%% ── Figure ───────────────────────────────────────────────────────────────
fig = figure('Units','centimeters','Position',[0 0 19 14],'Color','w');

ax = axes(fig,'Position',[0.08 0.58 0.38 0.36]); hold(ax,'on');
for s = 1:numel(skeys)
    scatter(ax, SA(:,s,2), SA(:,s,1), 12, scols{s}, 'filled', 'MarkerFaceAlpha', .6, ...
            'DisplayName', snames{s});
end
xlabel(ax,'Envelope CV (carrier band)'); ylabel(ax,'specparam exponent');
title(ax,sprintf('Table 1 scenarios (n = %d each)', NR),'FontSize',8,'FontWeight','bold');
legend(ax,'Location','southwest','FontSize',6.5,'Box','off','TextColor','k');   % lower left is empty
set_ax(ax); plabel(ax,'A');

ax = axes(fig,'Position',[0.58 0.58 0.38 0.36]); hold(ax,'on');
for c = 1:4
    m = y == c;
    scatter(ax, X(m,6), X(m,1), 12, ccols{c}, 'filled', 'MarkerFaceAlpha', .6, ...
            'DisplayName', classes(c).name);
end
xlabel(ax,'Envelope CV (8-12 Hz)'); ylabel(ax,'specparam exponent');
title(ax,'Toy 2\times2: dynamics \times broadband','FontSize',8,'FontWeight','bold');
yr = [min(X(:,1)) max(X(:,1))]; ylim(ax, [yr(1) - 0.45*diff(yr), yr(2) + 0.05*diff(yr)]);   % empty band at the bottom for the legend
legend(ax,'Location','southeast','NumColumns',2,'FontSize',6.5,'Box','off','TextColor','k');
set_ax(ax); plabel(ax,'B');

ax = axes(fig,'Position',[0.20 0.10 0.21 0.34]); hold(ax,'on');   % left margin holds the feature-set labels
barh(ax, 1:numel(fsets), acc, 0.6, 'FaceColor', [.5 .5 .5], 'EdgeColor', 'none');
xline(ax, 0.25, ':', 'chance', 'Color', [.3 .3 .3], 'FontSize', 6.5, ...
      'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');
for q = 1:numel(fsets)
    text(ax, acc(q) + 0.02, q, sprintf('%.2f', acc(q)), 'HorizontalAlignment','left', ...
         'VerticalAlignment','middle','FontSize',7,'Color','k');
end
set(ax,'YTick',1:numel(fsets),'YTickLabel',{fsets.name},'YDir','reverse');
xlim(ax,[0 1.15]); ylim(ax,[0.4 numel(fsets)+0.6]); xlabel(ax,'LOO accuracy');
title(ax,'Toy: nearest centroid, leave-one-out','FontSize',8,'FontWeight','bold');
set_ax(ax);
text(ax,-0.75,1.06,'C','Units','normalized','FontSize',12,'FontWeight','bold', ...
     'VerticalAlignment','bottom','Color','k');   % letter aligned with the figure's left edge

cm_show = [1 4];
for q = 1:2
    ax = axes(fig,'Position',[0.47 + (q-1)*0.27, 0.10, 0.20, 0.30]);
    C = Cm{cm_show(q)};
    imagesc(ax, C); colormap(ax, flipud(gray(256))); clim(ax, [0 NR]);
    for i = 1:4
        for j = 1:4
            tc = 'k'; if C(i,j) > 0.6*NR, tc = 'w'; end
            text(ax, j, i, sprintf('%d', C(i,j)), 'HorizontalAlignment','center', ...
                 'FontSize', 7, 'Color', tc);
        end
    end
    set(ax,'XTick',1:4,'YTick',1:4,'XTickLabel',{classes.short},'YTickLabel',{classes.short}, ...
           'XTickLabelRotation',35,'FontSize',7);
    xlabel(ax,'assigned'); if q == 1, ylabel(ax,'simulated'); end
    title(ax, fsets(cm_show(q)).name, 'FontSize', 8, 'FontWeight', 'bold');
    set_ax(ax);
    if q == 1, plabel(ax,'D'); end
end

%% ── Caption statistics ───────────────────────────────────────────────────
fprintf('\n--- fig_toy caption statistics ---\n');
fprintf('A: Table 1 scenarios, exponent and carrier-band CV (mean +/- SD, n = %d)\n', NR);
for s = 1:numel(skeys)
    fprintf('  %-22s exp %+.2f +/- %.2f   CV %.2f +/- %.2f\n', skeys{s}, ...
        mean(SA(:,s,1)), std(SA(:,s,1)), mean(SA(:,s,2)), std(SA(:,s,2)));
end
fprintf('B-D: toy classes, feature means +/- SD (n = %d)\n', NR);
for c = 1:4
    m = y == c;
    fprintf('  %-18s', classes(c).name);
    for j = 1:numel(fnames)
        fprintf('  %s %+.2f+/-%.2f', fnames{j}, mean(X(m,j)), std(X(m,j)));
    end
    fprintf('\n');
end
dyn = [1 1 2 2]; bb = [1 2 1 2];              % factor level of each class
for q = 1:numel(fsets)
    C = Cm{q}; nT = sum(C(:)); aD = 0; aB = 0;
    for i = 1:4
        for j = 1:4
            aD = aD + C(i,j) * (dyn(i) == dyn(j));
            aB = aB + C(i,j) * (bb(i)  == bb(j));
        end
    end
    fprintf('LOO accuracy, %-22s %.2f   (dynamics factor %.2f, broadband factor %.2f)\n', ...
        fsets(q).name, acc(q), aD/nT, aB/nT);
    disp(C);
end
fprintf('TOY: labels exist only because the signals are simulated.\n---\n\n');
save_fig(fig, cfg, 'fig_toy');
end

%% ── Local helpers ────────────────────────────────────────────────────────
function f = toy_features(x, fs)
% [exponent offset peak_cf peak_pw peak_bw envelope_CV SPRiNT_SD]
[fq, pw] = n1f.compute_psd(x, fs);
sp = n1f.run_specparam(fq, pw, [1 45]);
f = [sp.exponent sp.offset NaN 0 0 NaN NaN];
if sp.n_peaks > 0
    [~, kk] = min(abs([sp.peaks.cf] - 10));
    f(3:5) = [sp.peaks(kk).cf sp.peaks(kk).pw sp.peaks(kk).bw];
end
ev = n1f.envelope_cv(x, fs, [8 12]);
f(6) = ev.cv;
spr = n1f.run_sprint(x, fs);
f(7) = std(spr.exponent, 'omitnan');
end

function [acc, C] = loo_centroid(X, y, cols)
% Leave-one-out nearest-centroid on features z-scored with training-fold stats.
n = numel(y); K = max(y); pred = nan(n,1);
for i = 1:n
    tr = true(n,1); tr(i) = false;
    Xt = X(tr, cols); mu = mean(Xt,1); sd = std(Xt,0,1); sd(sd == 0) = 1;
    Z = (Xt - mu) ./ sd; z = (X(i,cols) - mu) ./ sd;
    yt = y(tr); d = inf(K,1);
    for c = 1:K, d(c) = sum((mean(Z(yt == c,:),1) - z).^2); end
    [~, pred(i)] = min(d);
end
acc = mean(pred == y);
C = zeros(K);
for i = 1:n, C(y(i), pred(i)) = C(y(i), pred(i)) + 1; end
end

function set_ax(ax)
set(ax,'Box','off','TickDir','out','LineWidth',.7, ...
    'Color','w','XColor','k','YColor','k');
if ax.FontSize > 8, ax.FontSize = 8; end
ax.Title.Color = 'k';
end

function plabel(ax, letter)
text(ax,-0.14,1.06,letter,'Units','normalized','FontSize',12, ...
     'FontWeight','bold','VerticalAlignment','bottom','Color','k');
end

function save_fig(fig, cfg, name)
fname = fullfile(cfg.out_dir,[name '.png']);
exportgraphics(fig, fname, 'Resolution', 150);
fprintf('  Saved: %s\n', fname); close(fig);
end
