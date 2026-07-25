function sig = make_signal(freqs, varargin)
% make_signal  Oscillation(s) with optional sinusoidal AM and/or 1/f background.
%
%   sig = n1f.make_signal(10)
%   sig = n1f.make_signal([10 20], 'am_freqs',10, 'true_1f',true, ...
%                         'noise_amplitude',0.3, 'seed',42)
p = inputParser;
addRequired(p,'freqs');
addParameter(p,'am_freqs',        []);
addParameter(p,'true_1f',         false);
addParameter(p,'noise_amplitude', 0);
addParameter(p,'seed',            42);
parse(p, freqs, varargin{:});
r = p.Results;

cfg = n1f.config();
N   = cfg.duration * cfg.fs;
t   = (0:N-1)' / cfg.fs;
rng(r.seed);

sig = zeros(N,1);
for f = r.freqs(:)'
    ph  = 2*pi*rand;
    if ismember(f, r.am_freqs)
        env = 0.5 + 0.5*cos(2*pi*cfg.am_freq*t);
        sig = sig + env .* cos(2*pi*f*t + ph);
    else
        sig = sig + cos(2*pi*f*t + ph);
    end
end

if r.true_1f
    pink = n1f.sim_powerlaw(N, cfg.fs, cfg.exp_1f, r.seed+1000);
    sig  = sig + pink/std(pink);
end

if r.noise_amplitude > 0
    sig = sig + r.noise_amplitude * randn(N,1);
end
end
