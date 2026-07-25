function sig = sim_powerlaw(N, fs, exponent, seed)
% sim_powerlaw  Signal with a power-law spectrum (e.g. exponent=-1 for pink noise).
if nargin < 4, seed = 0; end
rng(seed);
wn      = randn(N, 1);
f       = (0:N-1)' / N * fs;
f(1)    = 1;
H       = f .^ (exponent / 2);
H(1)    = 0;
sig_f   = fft(wn) .* H;
sig     = real(ifft(sig_f));
sig     = sig - mean(sig);
end
