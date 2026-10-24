function x = sim_narrowband_osc(N, fs, f0, B, shape)
% sim_narrowband_osc  Oscillation with random amplitude and phase.
%
%   x = n1f.sim_narrowband_osc(N, fs, f0, B)            % shape 'gauss'
%   x = n1f.sim_narrowband_osc(N, fs, f0, B, shape)
%
%   Narrowband Gaussian process: complex Gaussian noise (the complex
%   envelope) low-passed and shifted to the carrier f0. Its envelope is
%   Rayleigh-distributed (Rice, 1944). B = 0 returns a stationary sinusoid
%   with random phase. Output has unit variance. Draws from the global RNG
%   stream.
%
%   shape = 'gauss'   (default) amplitude window exp(-0.5*(f/B)^2): B is the
%                     SD (Hz) of the window; power spectrum Gaussian with
%                     SD B/sqrt(2); envelope timescale ~1/B.
%   shape = 'lorentz' amplitude window 1/sqrt(1+(f/B)^2): power spectrum
%                     Lorentzian with half-width at half-maximum B (Hz);
%                     the complex envelope is an Ornstein-Uhlenbeck process
%                     with correlation time 1/(2*pi*B). Tails decay as
%                     1/(f-f0)^2 regardless of B.

if nargin < 5 || isempty(shape), shape = 'gauss'; end
t = (0:N-1)' / fs;
if B == 0
    x = sqrt(2) * cos(2*pi*f0*t + 2*pi*rand);
    return
end
z = randn(N,1) + 1i*randn(N,1);
f = (0:N-1)' / N * fs;
f(f >= fs/2) = f(f >= fs/2) - fs;              % two-sided frequency axis
switch lower(shape)
    case 'gauss',   W = exp(-0.5*(f/B).^2);
    case 'lorentz', W = 1 ./ sqrt(1 + (f/B).^2);
    otherwise, error('sim_narrowband_osc: unknown shape "%s"', shape);
end
env = ifft(fft(z) .* W);                        % complex envelope
x = real(env .* exp(1i*2*pi*f0*t));
x = x / std(x);
end
