function check_lorentz()
% Printout only (no figure). Exponent change produced by amplitude dynamics
% with a Lorentzian envelope spectrum (n1f.sim_narrowband_osc 'lorentz',
% half-width B), relative to a stationary sinusoid of the same power on the
% same 1/f background (exponent -1, variance rho x oscillatory) + white noise
% 0.02. Supports the qualifier on the Gaussian map in Section 3.4: with a
% long-tailed envelope spectrum, slow fluctuations also change the exponent.
% Python validation (n = 20): rho = 1, +0.07 for B = 0.25-2 Hz; rho = 0.1,
% +0.10 to +0.13.

cfg = n1f.config();
fs = cfg.fs; N = cfg.duration * fs; f0 = 10; NR = 30;
Bs = [0.25 0.5 1 2]; rhos = [0.1 1];
fprintf('\n--- check_lorentz (n = %d per cell; Lorentzian half-width B) ---\n', NR);
fprintf('%5s %6s  %24s  %18s  %18s\n', 'rho', 'B (Hz)', 'exp change (signrank p)', 'peak bw (Hz)', 'peak power');
for i = 1:numel(rhos)
    for j = 1:numel(Bs)
        d = nan(NR,1); bw = nan(NR,1); pw = nan(NR,1);
        for k = 1:NR
            pl = n1f.sim_powerlaw(N, fs, -1, 4000 + k); pl = pl / std(pl);   % reseeds internally
            rng(5000 + k);
            fl = sqrt(rhos(i)) * pl + 0.02*randn(N,1);
            x0 = n1f.sim_narrowband_osc(N, fs, f0, 0) + fl;
            x1 = n1f.sim_narrowband_osc(N, fs, f0, Bs(j), 'lorentz') + fl;
            [fq, p0] = n1f.compute_psd(x0, fs); s0 = n1f.run_specparam(fq, p0, [1 45]);
            [~,  p1] = n1f.compute_psd(x1, fs); s1 = n1f.run_specparam(fq, p1, [1 45]);
            d(k) = s1.exponent - s0.exponent;
            if s1.n_peaks > 0
                [~, kk] = min(abs([s1.peaks.cf] - f0));
                bw(k) = s1.peaks(kk).bw; pw(k) = s1.peaks(kk).pw;
            end
        end
        fprintf('%5g %6g  %+8.3f +/- %.3f (p=%.2g)  %6.2f +/- %.2f   %6.2f +/- %.2f\n', rhos(i), Bs(j), ...
            mean(d), std(d), signrank(d), mean(bw,'omitnan'), std(bw,'omitnan'), mean(pw,'omitnan'), std(pw,'omitnan'));
    end
end
fprintf('---\n\n');
end
