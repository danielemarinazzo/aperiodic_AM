function result = envelope_cv(sig, fs, band)
% envelope_cv  Hilbert amplitude envelope CV in a frequency band [lo hi] Hz.
%
% Returns observed CV and excess CV relative to expected background.
%
% Background CV for bandpass-filtered 1/f processes is approximately
% constant regardless of aperiodic exponent:
%   theta (4-8 Hz):   CV_bg ~ 0.487
%   alpha (8-12 Hz):  CV_bg ~ 0.487
%   beta (16-30 Hz):  CV_bg ~ 0.310
%
% Excess CV = CV_obs - CV_bg:
%   ~0   : consistent with pure 1/f background
%   < 0  : sustained oscillation (suppresses envelope variability)
%   > 0  : amplitude-modulated oscillation or burst
%
% This makes CV_excess modality- and state-invariant because the reference
% is the expected value for the broadband background, not an absolute threshold.

lo = band(1); hi = band(2);
[b,a]   = butter(4, [lo hi]/(fs/2), 'bandpass');
sig_bp  = filtfilt(b, a, sig(:));
env     = abs(hilbert(sig_bp));
[b2,a2] = butter(4, 3/(fs/2), 'low');
env_sm  = filtfilt(b2, a2, env);

mu = mean(env_sm);
sd = std(env_sm);
cv = sd / (mu + eps);

% Expected background CV (from simulation of 1/f signals, N=500, T=2-10s)
% Stable across exponents 0.5-3.0 and segment lengths 2-10s
if lo >= 16
    cv_bg = 0.310;   % beta band
else
    cv_bg = 0.487;   % theta and alpha bands
end

result.cv         = cv;
result.cv_excess  = cv - cv_bg;   % key: negative=osc, ~0=1/f, positive=AM
result.cv_bg      = cv_bg;
result.mean_amp   = mu;
result.std_amp    = sd;
result.envelope   = env_sm;
end