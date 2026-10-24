function figS2_knee()
% Supplementary figure — does the knee model absorb the effect of amplitude
% dynamics? (Reviewer 2.)
%
% Same generic model as fig_modulation: static (stationary sinusoid) vs
% dynamic (narrowband process, bandwidth B; n1f.sim_narrowband_osc) 10 Hz
% oscillation on the same floor (1/f exponent -1 with variance rho x the
% oscillatory variance, + white noise 0.02). Each pair is fitted with the
% fixed and the knee aperiodic model (n1f.run_specparam, [1 45] Hz).
%   A  exponent shift (dynamic - static) under each model
%   B  knee frequency f_knee = kappa^(1/chi) for static and dynamic
% rho = 0 is omitted: without a broadband floor the knee fit to the static
% oscillation is undefined (non-positive knee or exponent).
% In the knee model the exponent describes the slope above the knee.

cfg = n1f.config();
fs  = cfg.fs;
N   = cfg.duration * fs;
NR  = 30;
f0  = 10;
Bs   = [1 3 6 10]';                % column, to match per-B column statistics
rhos = [0.1 1];
cr   = {[.84 .37 .00], [.17 .50 .72]};

% R(i, j, k, :) = [exp_fix_st exp_fix_dy exp_knee_st exp_knee_dy fk_st fk_dy r2_fix_dy r2_knee_dy]
R = nan(numel(rhos), numel(Bs), NR, 8);
for i = 1:numel(rhos)
    for j = 1:numel(Bs)
        for k = 1:NR
            seed = 20000 + 1000*i + 100*j + k;
            pl = n1f.sim_powerlaw(N, fs, -1, seed); pl = pl / std(pl);
            rng(seed + 50000);
            fl = sqrt(rhos(i)) * pl + 0.02*randn(N,1);
            xs = n1f.sim_narrowband_osc(N, fs, f0, 0)     + fl;
            xd = n1f.sim_narrowband_osc(N, fs, f0, Bs(j)) + fl;
            [efs, ~,   ~   ] = fit_one(xs, fs, 'fixed');
            [efd, ~,   r2f ] = fit_one(xd, fs, 'fixed');
            [eks, fks, ~   ] = fit_one(xs, fs, 'knee');
            [ekd, fkd, r2k ] = fit_one(xd, fs, 'knee');
            R(i,j,k,:) = [efs efd eks ekd fks fkd r2f r2k];
        end
    end
    fprintf('  rho = %g done\n', rhos(i));
end

%% ── Figure ───────────────────────────────────────────────────────────────
fig = figure('Units','centimeters','Position',[0 0 18 7],'Color','w');
ax = axes(fig,'Position',[0.08 0.18 0.38 0.68]); hold(ax,'on');
for i = 1:numel(rhos)
    dfix  = squeeze(R(i,:,:,2) - R(i,:,:,1));
    dknee = squeeze(R(i,:,:,4) - R(i,:,:,3));
    errorbar(ax, Bs, mean(dfix,2,'omitnan'),  std(dfix,0,2,'omitnan'),  '-o', 'Color', cr{i}, ...
        'MarkerSize', 3, 'MarkerFaceColor', cr{i}, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('fixed, \\rho = %g', rhos(i)));
    errorbar(ax, Bs, mean(dknee,2,'omitnan'), std(dknee,0,2,'omitnan'), '--s', 'Color', cr{i}, ...
        'MarkerSize', 3, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('knee, \\rho = %g', rhos(i)));
end
yline(ax, 0, ':', 'Color', [.6 .6 .6], 'HandleVisibility', 'off');
set(ax,'XScale','log','XTick',Bs'); xlim(ax,[0.8 12]);
xlabel(ax,'Bandwidth B of amplitude/phase fluctuations (Hz)');
ylabel(ax,'\Delta exponent (dynamic - static)');
title(ax,sprintf('Exponent shift by aperiodic model (mean \\pm SD, n = %d)', NR),'FontSize',8,'FontWeight','bold');
legend(ax,'Location','northwest','FontSize',6.5,'Box','off','TextColor','k');
set_ax(ax); plabel(ax,'A');

ax = axes(fig,'Position',[0.58 0.18 0.39 0.68]); hold(ax,'on');
for i = 1:numel(rhos)
    fks = squeeze(R(i,:,:,5)); fkd = squeeze(R(i,:,:,6));
    plot(ax, Bs, median(fks,2,'omitnan'), ':o', 'Color', cr{i}, 'MarkerSize', 3, 'LineWidth', 1, ...
        'DisplayName', sprintf('static, \\rho = %g', rhos(i)));
    q = prctile(fkd, [25 50 75], 2);
    errorbar(ax, Bs, q(:,2), q(:,2)-q(:,1), q(:,3)-q(:,2), '-o', 'Color', cr{i}, ...
        'MarkerSize', 3, 'MarkerFaceColor', cr{i}, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('dynamic, \\rho = %g', rhos(i)));
end
set(ax,'XScale','log','XTick',Bs'); xlim(ax,[0.8 12]);
xlabel(ax,'Bandwidth B of amplitude/phase fluctuations (Hz)');
ylabel(ax,'Knee frequency (Hz)');
title(ax,'Knee frequency (median, IQR for dynamic)','FontSize',8,'FontWeight','bold');
legend(ax,'Location','northwest','FontSize',6.5,'Box','off','TextColor','k');
set_ax(ax); plabel(ax,'B');

%% ── Caption statistics ───────────────────────────────────────────────────
fprintf('\n--- figS_knee caption statistics (n = %d per cell) ---\n', NR);
fprintf('%5s %4s | %14s %6s | %14s %18s %6s %5s\n', 'rho', 'B', 'fixed dExp', 'R2', ...
        'knee dExp', 'f_knee st -> dy', 'R2', 'nNaN');
for i = 1:numel(rhos)
    for j = 1:numel(Bs)
        r = squeeze(R(i,j,:,:));
        dfx = r(:,2) - r(:,1); dkn = r(:,4) - r(:,3);
        fprintf('%5g %4g | %+7.2f +/- %4.2f %6.2f | %+7.2f +/- %4.2f %7.1f -> %6.1f %6.2f %5d\n', ...
            rhos(i), Bs(j), mean(dfx,'omitnan'), std(dfx,'omitnan'), mean(r(:,7),'omitnan'), ...
            mean(dkn,'omitnan'), std(dkn,'omitnan'), median(r(:,5),'omitnan'), median(r(:,6),'omitnan'), ...
            mean(r(:,8),'omitnan'), sum(isnan(r(:,6))));
    end
end
fprintf(['Timescale implied by the dynamic knee: tau = 1/(2*pi*f_knee).\n' ...
         'nNaN: dynamic fits with non-positive knee or exponent (f_knee undefined).\n---\n\n']);
save_fig(fig, cfg, 'figS_knee');
end

%% ── Local helpers ────────────────────────────────────────────────────────
function [e, fk, r2] = fit_one(x, fs, mode)
[fq, pw] = n1f.compute_psd(x, fs);
sp = n1f.run_specparam(fq, pw, [1 45], 'aperiodic_mode', mode);
e = sp.exponent; r2 = sp.r_squared; fk = NaN;
if strcmp(mode, 'knee') && sp.knee > 0 && sp.exponent > 0
    fk = sp.knee^(1/sp.exponent);     % sp.knee is the knee parameter kappa
end
end

function set_ax(ax)
set(ax,'Box','off','TickDir','out','FontSize',8,'LineWidth',.7, ...
    'Color','w','XColor','k','YColor','k');
ax.Title.Color = 'k';
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
