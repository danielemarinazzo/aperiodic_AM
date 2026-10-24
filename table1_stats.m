function table1_stats(NR)
% Table 1 with variability. The five illustrative scenarios of
% run_diagnostics.m (make_scenario), regenerated over NR independent
% realizations (default 30); features as in run_diagnostics.m
% (compute_features, short-signal branch). Only the seed varies between
% realizations. Prints mean +/- SD per scenario and feature, and writes every
% realization to <cfg.out_dir>/table1_realizations.csv.
%
% Scenario definitions: n1f.make_scenario (copied from run_diagnostics.m).
%
% The CV columns are only meaningful in a band that contains the carrier
% (Section 2.5); values in empty bands reflect filtered noise.

if nargin < 1 || isempty(NR), NR = 30; end
cfg = n1f.config();
fs  = cfg.fs;
ws  = warning('off','all');

names = {'A. Pure alpha', 'B. AM alpha', 'C. Alpha on 1/f', 'D. Bursty beta', 'E. Bursty beta on 1/f'};
keys  = {'pure_alpha', 'am_alpha', 'alpha_plus_1f', 'burst_beta', 'burst_plus_1f'};
fn    = {'exp', 'r2', 'cv_alpha', 'cv_beta', 'sprint_sd', 'fbosc_peak', 'fbosc_peak_f'};

V = nan(NR, numel(keys), numel(fn));
for s = 1:numel(keys)
    for k = 1:NR
        sig = n1f.make_scenario(keys{s}, 1000 + k);
        f = features(sig, fs);
        for j = 1:numel(fn), V(k,s,j) = f.(fn{j}); end
    end
    fprintf('  %s done\n', names{s});
end
warning(ws);

%% ── Print ────────────────────────────────────────────────────────────────
fprintf('\n--- Table 1: mean +/- SD over %d realizations ---\n', NR);
fprintf('%-24s %16s %14s %14s %14s %14s %22s\n', 'Scenario', 'exp', 'R^2', ...
        'CV alpha', 'CV beta', 'SPRiNT SD', 'fBOSC peak @ Hz');
for s = 1:numel(keys)
    m  = squeeze(mean(V(:,s,:), 1, 'omitnan'));
    sd = squeeze(std(V(:,s,:), 0, 1, 'omitnan'));
    fprintf('%-24s %+7.2f+/-%-5.2f %6.2f+/-%-4.2f %6.2f+/-%-4.2f %6.2f+/-%-4.2f %6.2f+/-%-4.2f %6.2f+/-%-4.2f @ %4.1f\n', ...
        names{s}, m(1), sd(1), m(2), sd(2), m(3), sd(3), m(4), sd(4), m(5), sd(5), ...
        m(6), sd(6), median(V(:,s,7), 'omitnan'));
end
fprintf('(fBOSC peak frequency: median over realizations)\n');
fprintf('Exponent range per scenario [min, max]:\n');
for s = 1:numel(keys)
    fprintf('  %-24s [%+.2f, %+.2f]\n', names{s}, min(V(:,s,1)), max(V(:,s,1)));
end
fprintf('---\n\n');

%% ── CSV of all realizations ─────────────────────────────────────────────
scen = repelem(names(:), NR, 1);
rid  = repmat((1:NR)', numel(keys), 1);
M    = reshape(V, NR*numel(keys), numel(fn));   % rows: realization fastest, then scenario
Tout = [table(scen, rid, 'VariableNames', {'scenario', 'realization'}), ...
        array2table(M, 'VariableNames', fn)];
fname = fullfile(cfg.out_dir, 'table1_realizations.csv');
writetable(Tout, fname);
fprintf('  Saved: %s\n', fname);
end

%% ── Local helpers ────────────────────────────────────────────────────────
function f = features(sig, fs)
[fq, pw] = n1f.compute_psd(sig, fs);
sp = n1f.run_specparam(fq, pw, [1 45]);
f.exp = sp.exponent;
f.r2  = sp.r_squared;
ea = n1f.envelope_cv(sig, fs, [8 12]);
eb = n1f.envelope_cv(sig, fs, [16 30]);
f.cv_alpha = ea.cv;
f.cv_beta  = eb.cv;
spr = n1f.run_sprint(sig, fs);
f.sprint_sd = std(spr.exponent, 'omitnan');
fb = n1f.run_fbosc(sig, fs, 'freqs', (4:44)', 'freq_range', [4 44]);
[f.fbosc_peak, i] = max(fb.abundance);
f.fbosc_peak_f = fb.freqs(i);
end
