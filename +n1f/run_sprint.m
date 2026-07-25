function result = run_sprint(sig, fs, varargin)
% run_sprint  SPRiNT: sliding-window specparam on a continuous signal.
%
%   result = n1f.run_sprint(sig, fs)
%   result = n1f.run_sprint(sig, fs, 'win_s',1.0, 'n_avg',5, 'freq_range',[1 45])
%
%   Algorithm (matches Python n1f.sprint):
%     1. Slide a Hann window of win_s seconds across the signal with 50%
%        overlap, computing a PSD per window via pwelch.
%     2. Average n_avg *consecutive* PSDs to smooth each output bin
%        (sliding average, not block average): output bin j averages
%        windows j..j+n_avg-1, so adjacent bins share windows.
%     3. Fit specparam to each smoothed PSD.
%
%   This gives n_wins - n_avg + 1 output bins, NOT floor(n_wins / n_avg).

p = inputParser;
addParameter(p,'win_s',      1.0);   % 1 s window matches Python default; 0 Hz bin stripped in run_specparam
addParameter(p,'n_avg',      5);
addParameter(p,'freq_range', [1 45]);
parse(p, varargin{:});
r = p.Results;

win_n  = round(r.win_s * fs);
step_n = round(win_n / 2);                  % 50% overlap
N      = length(sig);
starts = 1 : step_n : (N - win_n + 1);
n_wins = length(starts);

if n_wins < r.n_avg
    error('Signal too short for SPRiNT: need at least %.1f s', ...
          (r.n_avg * step_n + win_n) / fs);
end

% ── Step 1: compute a PSD per window ─────────────────────────────────────────
[freqs, ~] = n1f.compute_psd(sig(1:win_n), fs, win_n);
psds = zeros(length(freqs), n_wins);
for i = 1:n_wins
    chunk = sig(starts(i) : starts(i)+win_n-1);
    [~, psds(:,i)] = n1f.compute_psd(chunk, fs, win_n);
end

% ── Step 2: sliding average over n_avg consecutive PSDs ──────────────────────
n_bins   = n_wins - r.n_avg + 1;
exponent = nan(n_bins, 1);
offset   = nan(n_bins, 1);
r_sq     = nan(n_bins, 1);
times    = nan(n_bins, 1);

half       = floor(r.n_avg / 2);
n_failed   = 0;
first_err  = '';

for b = 1:n_bins
    cols    = b : b + r.n_avg - 1;
    psd_avg = mean(psds(:, cols), 2);
    t_mid   = (starts(b + half) + win_n/2 - 1) / fs;
    times(b) = t_mid;   % store time regardless of fit success
    try
        sp = n1f.run_specparam(freqs, psd_avg, r.freq_range);
        exponent(b) = sp.exponent;
        offset(b)   = sp.offset;
        r_sq(b)     = sp.r_squared;
    catch ME
        n_failed = n_failed + 1;
        if isempty(first_err), first_err = ME.message; end
    end
end

if n_failed > 0
    warning('n1f:sprint', ...
        'SPRiNT: %d/%d specparam fits failed. First error:\n  %s', ...
        n_failed, n_bins, first_err);
end

result.times     = times;
result.exponent  = exponent;
result.offset    = offset;
result.r_squared = r_sq;
end