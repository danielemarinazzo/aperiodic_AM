function result = run_erp_sprint(trials, times, fs, varargin)
% run_erp_sprint  Trial-averaged wavelet specparam for short ERP epochs.
%
%   trials   (n_trials x n_samples) double
%   times    time axis in seconds (1 x n_samples or n_samples x 1)
%   fs       sampling frequency (Hz)
%
%   Optional parameters:
%     n_cycles    Morlet cycles            (default 6)
%     time_step   Output bin step (s)      (default 0.010)
%     freq_range  specparam fitting range  (default [2 40])
%     freqs       analysis frequencies     (default (2:40)')
%
%   result.times     output time vector (s)
%   result.freqs     analysis frequencies (Hz)
%   result.full_fit  log10 modeled PSD (n_freqs x n_times), NaN at edges
%   result.aperiodic log10 aperiodic   (n_freqs x n_times), NaN at edges
%   result.exponent  aperiodic exponent timecourse (NaN at edges)
%   result.r_squared R-squared timecourse

p = inputParser;
addParameter(p,'n_cycles',   6);
addParameter(p,'time_step',  0.010);
addParameter(p,'freq_range', [2 40]);
addParameter(p,'freqs',      (2:40)');
parse(p, varargin{:});
r = p.Results;

times    = times(:)';
freqs_an = r.freqs(:);
n_freqs  = length(freqs_an);
[n_trials, n_samp] = size(trials);

% ── 1. Trial-averaged Morlet TFR ─────────────────────────────────────────────
tfr_avg = zeros(n_freqs, n_samp);
for fi = 1:n_freqs
    f0     = freqs_an(fi);
    sigma  = r.n_cycles / (2*pi*f0);
    t_win  = (-3*sigma : 1/fs : 3*sigma)';
    wav    = exp(2i*pi*f0*t_win) .* exp(-t_win.^2 / (2*sigma^2));
    wav    = wav / norm(wav);
    pw_sum = zeros(1, n_samp);
    for tr = 1:n_trials
        c      = conv(trials(tr,:), conj(wav)', 'same');
        pw_sum = pw_sum + abs(c).^2;
    end
    tfr_avg(fi,:) = pw_sum / n_trials;
end

% ── 2. Output time grid ───────────────────────────────────────────────────────
step_n = max(1, round(r.time_step * fs));
t_idx  = 1 : step_n : n_samp;
t_out  = times(t_idx);
n_out  = length(t_out);

% ── 3. specparam at each time bin ────────────────────────────────────────────
full_fit  = nan(n_freqs, n_out);
aperiodic = nan(n_freqs, n_out);
exponent  = nan(1, n_out);
r_sq_out  = nan(1, n_out);
peak_cf   = nan(1, n_out);   % dominant-peak centre frequency (Hz), no gate
peak_pw   = nan(1, n_out);   % dominant-peak power (log10), aperiodic-corrected

for ti = 1:n_out
    psd_t = tfr_avg(:, t_idx(ti));
    try
        sp = n1f.run_specparam(freqs_an, psd_t, r.freq_range);
        % Interpolate model back to freqs_an grid
        full_fit(:,ti)  = interp1(sp.freqs, sp.modeled,   freqs_an, 'linear', 'extrap');
        aperiodic(:,ti) = interp1(sp.freqs, sp.aperiodic, freqs_an, 'linear', 'extrap');
        exponent(ti)    = sp.exponent;
        r_sq_out(ti)    = sp.r_squared;
        % Dominant fitted peak (largest pw), no frequency gate: this is the
        % oscillatory power/centre a user would read from specparam directly.
        if sp.n_peaks > 0
            [pw_max, kmax] = max([sp.peaks.pw]);
            peak_pw(ti) = pw_max;
            peak_cf(ti) = sp.peaks(kmax).cf;
        end
    catch
    end
end

% ── 4. Frequency-dependent edge masking ───────────────────────────────────────
t0 = times(1); t1 = times(end);
fwhm_s = @(f) r.n_cycles * sqrt(2*log(2)) / (pi * f);
for fi = 1:n_freqs
    hw  = fwhm_s(freqs_an(fi)) / 2;
    bad = t_out < (t0+hw) | t_out > (t1-hw);
    full_fit(fi,bad)  = NaN;
    aperiodic(fi,bad) = NaN;
end
hw_exp = fwhm_s(10) / 2;
bad_t  = t_out < (t0+hw_exp) | t_out > (t1-hw_exp);
exponent(bad_t)  = NaN;
r_sq_out(bad_t)  = NaN;
peak_cf(bad_t)   = NaN;
peak_pw(bad_t)   = NaN;

result.times      = t_out;
result.freqs      = freqs_an;
result.full_fit   = full_fit;
result.aperiodic  = aperiodic;
result.exponent   = exponent;
result.r_squared  = r_sq_out;
result.peak_cf    = peak_cf;   % dominant-peak centre frequency (Hz)
result.peak_pw    = peak_pw;   % dominant-peak power (log10)
result.tfr_avg    = tfr_avg(:, t_idx);   % subsampled to output time grid
result.params     = struct('freq_range',r.freq_range,'n_cycles',r.n_cycles);
end