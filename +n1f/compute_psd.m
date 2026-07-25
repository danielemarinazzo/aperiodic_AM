function [freqs, psd] = compute_psd(sig, fs, nperseg)
% compute_psd  Welch PSD via pwelch.
if nargin < 3, nperseg = round(fs*2); end
nperseg  = min(nperseg, length(sig));
noverlap = round(nperseg/2);
[psd, freqs] = pwelch(sig(:), hann(nperseg), noverlap, nperseg, fs);
freqs = freqs(:); psd = psd(:);
end
