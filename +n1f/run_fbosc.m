function result = run_fbosc(sig, fs, varargin)
% run_fbosc  fBOSC: specparam-calibrated oscillatory burst detection.
%
%   result = n1f.run_fbosc(sig, fs)
%   result = n1f.run_fbosc(sig, fs, 'freqs',4:44, 'freq_range',[4 44], 'n_cycles',6)
%
%   Uses a Morlet TFR, fits specparam to the mean TFR power to obtain the
%   aperiodic baseline, then thresholds at the 95th percentile of the
%   chi-squared(2) distribution times the baseline power.

p = inputParser;
addParameter(p,'freqs',      (4:44)');
addParameter(p,'freq_range', [4 44]);
addParameter(p,'n_cycles',   6);
parse(p, varargin{:});
r = p.Results;

freqs_an = r.freqs(:);
n_freqs  = length(freqs_an);
N        = length(sig);
t        = (0:N-1)' / fs;

% ── 1. Morlet TFR ────────────────────────────────────────────────────────────
tfr = zeros(n_freqs, N);
for fi = 1:n_freqs
    f0    = freqs_an(fi);
    sigma = r.n_cycles / (2*pi*f0);
    t_win = (-3*sigma : 1/fs : 3*sigma)';
    wav   = exp(2i*pi*f0*t_win) .* exp(-t_win.^2 / (2*sigma^2));
    wav   = wav / norm(wav);
    c     = conv(sig(:), conj(wav), 'same');
    tfr(fi,:) = abs(c).^2;
end

% ── 2. specparam on mean TFR power -> aperiodic baseline ─────────────────────
mean_pow = mean(tfr, 2);
try
    sp = n1f.run_specparam(freqs_an, mean_pow, r.freq_range);
    ap_pow = 10.^interp1(sp.freqs, sp.aperiodic, freqs_an, 'linear', 'extrap');
catch
    % Fallback: log-log linear fit
    lf = log10(freqs_an); lp = log10(mean_pow+eps);
    c  = polyfit(lf, lp, 1);
    ap_pow = 10.^polyval(c, lf);
end

% chi-squared(2) 95th percentile threshold
chi2_thresh = 5.991;
threshold   = chi2_thresh * ap_pow / 2;

% ── 3. Detect bursts ──────────────────────────────────────────────────────────
detected  = false(n_freqs, N);
for fi = 1:n_freqs
    above     = tfr(fi,:) > threshold(fi);
    min_samps = round(3 / freqs_an(fi) * fs);
    labeled   = bwlabel(above);
    for region = 1:max(labeled)
        idx = find(labeled == region);
        if length(idx) >= min_samps
            detected(fi,idx) = true;
        end
    end
end

result.freqs     = freqs_an;
result.times     = t;
result.tfr       = tfr;
result.threshold = threshold;
result.detected  = detected;
result.abundance = mean(detected, 2);
result.sp        = sp;
end
