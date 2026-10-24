function figS1_duty()
% Supplementary figure — burst duty cycle does not set the exponent.
%
% Bursty 20 Hz oscillation with the cycle-level Markov switching of
% n1f.make_beta_burst (copied locally with p_enter/p_leave exposed, so the
% shared generator is not modified). Mean burst length fixed (p_leave = 0.1,
% ~10 cycles); p_enter set for target duty cycles; duty 1 = continuous.
%   D: bursts + white noise 0.02            E: bursts + 1/f (exponent -1, unit SD)
% Solid: burst amplitude scaled by 1/sqrt(realized duty), so total oscillatory
%   power is constant and only its temporal concentration changes.
% Dashed: unit burst amplitude (power grows with duty).
% NR realizations per point; Spearman correlation of realized duty vs exponent.

cfg = n1f.config();
fs  = cfg.fs;
N   = cfg.duration * fs;
t   = (0:N-1)' / fs;
NR  = 30;
duties = [0.1 0.2 0.3 0.5 0.7 0.9 1.0];
conds  = {'D','E'};
norms  = [true false];
cols   = struct('D', [.48 .00 .47], 'E', [.14 .55 .27]);

% E(c, n, d, k): exponent;  RD(c, n, d, k): realized duty
E  = nan(2, 2, numel(duties), NR); RD = E;
for c = 1:2
    for n = 1:2
        for d = 1:numel(duties)
            for k = 1:NR
                seed = 10000*c + 1000*n + 100*d + k;
                pl = n1f.sim_powerlaw(N, fs, -1, seed); pl = pl / std(pl);
                rng(seed + 50000);
                [b, rd] = burst_sig(duties(d), norms(n), N, fs, t);
                if strcmp(conds{c}, 'D'), x = b + 0.02*randn(N,1); else, x = b + pl; end
                [fq, pw] = n1f.compute_psd(x, fs);
                sp = n1f.run_specparam(fq, pw, [1 45]);
                E(c,n,d,k) = sp.exponent; RD(c,n,d,k) = rd;
            end
        end
    end
    fprintf('  condition %s done\n', conds{c});
end

%% ── Figure ───────────────────────────────────────────────────────────────
fig = figure('Units','centimeters','Position',[0 0 18 7],'Color','w');
ax = axes(fig,'Position',[0.06 0.18 0.36 0.68]); hold(ax,'on');
rng(7); [b1, ~] = burst_sig(0.2, true, N, fs, t);
rng(7); [b2, ~] = burst_sig(0.9, true, N, fs, t);
seg = t < 3;
plot(ax, t(seg), b1(seg) + 4, 'Color', [.85 .37 .01], 'LineWidth', .5);
plot(ax, t(seg), b2(seg) - 4, 'Color', [.17 .50 .72], 'LineWidth', .5);
text(ax, 3.04,  4, 'duty 0.2', 'FontSize', 7, 'Color', [.85 .37 .01], 'Clipping', 'off');   % right of the traces
text(ax, 3.04, -4, 'duty 0.9', 'FontSize', 7, 'Color', [.17 .50 .72], 'Clipping', 'off');
xlim(ax,[0 3]); set(ax,'YTick',[]); xlabel(ax,'Time (s)');
title(ax,'Bursty 20 Hz oscillation, equal total power','FontSize',8,'FontWeight','bold');
set_ax(ax); plabel(ax,'A');

ax = axes(fig,'Position',[0.53 0.18 0.44 0.68]); hold(ax,'on');
for c = 1:2
    for n = 1:2
        mu = squeeze(mean(E(c,n,:,:), 4)); sd = squeeze(std(E(c,n,:,:), 0, 4));
        col = cols.(conds{c});
        ls = '-'; if ~norms(n), ls = '--'; end
        if norms(n)
            fill(ax, [duties fliplr(duties)], [mu'-sd' fliplr(mu'+sd')], col, ...
                 'FaceAlpha', .12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        lbl = sprintf('%s, %s', conds{c}, ternary(norms(n), 'equal power', 'unit amplitude'));
        plot(ax, duties, mu, ls, 'Color', col, 'LineWidth', 1.3, 'Marker', 'o', ...
             'MarkerSize', 3, 'MarkerFaceColor', col, 'DisplayName', lbl);
    end
end
yline(ax, 0, ':', 'Color', [.6 .6 .6], 'HandleVisibility', 'off');
xlabel(ax,'Duty cycle (fraction of time in burst)'); ylabel(ax,'specparam exponent');
ylim(ax, [-0.6 1.45]);   % headroom above the E lines for the legend
title(ax,sprintf('Exponent vs duty cycle (n = %d per point)', NR), ...
      'FontSize',8,'FontWeight','bold');
legend(ax,'Location','northwest','NumColumns',2,'FontSize',6.5,'Box','off','TextColor','k');
set_ax(ax); plabel(ax,'B');

%% ── Caption statistics ───────────────────────────────────────────────────
fprintf('\n--- figS_duty caption statistics (n = %d per point) ---\n', NR);
for c = 1:2
    for n = 1:2
        e = E(c,n,:,:); r = RD(c,n,:,:);
        [rho, p] = corr(r(:), e(:), 'Type', 'Spearman');
        fprintf('%s, %-14s: Spearman rho = %+.2f, p = %.2g (n = %d)\n', conds{c}, ...
                ternary(norms(n), 'equal power', 'unit amplitude'), rho, p, numel(e));
        for d = 1:numel(duties)
            v = squeeze(E(c,n,d,:));
            fprintf('    duty %.1f: %+.2f +/- %.2f\n', duties(d), mean(v), std(v));
        end
    end
end
fprintf('---\n\n');
save_fig(fig, cfg, 'figS_duty');
end

%% ── Local helpers ────────────────────────────────────────────────────────
function [sig, d] = burst_sig(duty, equal_power, N, fs, t)
% Cycle-level Markov bursts as in n1f.make_beta_burst (20 Hz), with
% p_leave = 0.1 and p_enter set for the target duty. Uses the global stream.
freq = 20; p_leave = 0.10;
if duty >= 1, p_enter = 1; p_leave = 0; else, p_enter = duty*p_leave/(1 - duty); end
cyc_len = fs / freq; n_cycles = floor(N / cyc_len);
in_burst = false; is_osc = false(n_cycles, 1);
for c = 1:n_cycles
    if in_burst
        if rand < p_leave, in_burst = false; end
    else
        if rand < p_enter, in_burst = true; end
    end
    is_osc(c) = in_burst;
end
mask = false(N, 1);
for c = 1:n_cycles
    i1 = round((c-1)*cyc_len) + 1; i2 = min(round(c*cyc_len), N);
    mask(i1:i2) = is_osc(c);
end
d = mean(mask);
amp = 1;
if equal_power && d > 0, amp = 1/sqrt(d); end
sig = amp * mask .* cos(2*pi*freq*t + 2*pi*rand);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function set_ax(ax)
set(ax,'Box','off','TickDir','out','FontSize',8,'LineWidth',.7, ...
    'Color','w','XColor','k','YColor','k');
ax.Title.Color = 'k';
end

function plabel(ax, letter)
text(ax,-0.10,1.06,letter,'Units','normalized','FontSize',12, ...
     'FontWeight','bold','VerticalAlignment','bottom','Color','k');
end

function save_fig(fig, cfg, name)
fname = fullfile(cfg.out_dir,[name '.png']);
exportgraphics(fig, fname, 'Resolution', 150);
fprintf('  Saved: %s\n', fname); close(fig);
end
