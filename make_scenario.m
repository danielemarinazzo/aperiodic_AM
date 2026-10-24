function sig = make_scenario(name, seed)
% make_scenario  The five illustrative scenarios of Table 1, as defined in
% run_diagnostics.m (make_scenario), with a selectable seed.
%
%   sig = n1f.make_scenario(name, seed)
%
%   'pure_alpha'     cos(2*pi*10*t) + white noise 0.02
%   'am_alpha'       (0.5 + 0.5*cos(2*pi*1*t)) .* cos(2*pi*10*t) + white noise 0.02
%                    (envelope mean 0.5, range 0-1, i.e. 100% modulation; CV = 1/sqrt(2))
%   'alpha_plus_1f'  cos(2*pi*10*t) + 1/f (exponent -1, unit SD) + white noise 0.02
%   'burst_beta'     n1f.make_beta_burst + white noise 0.02
%   'burst_plus_1f'  n1f.make_beta_burst_with_1f (no white noise, as in run_diagnostics)
%
% n1f.sim_powerlaw and n1f.make_beta_burst call rng(seed) internally; the
% stream is reseeded where needed so each scenario is reproducible per seed.

cfg = n1f.config();
fs  = cfg.fs;
N   = cfg.duration * fs;
t   = (0:N-1)' / fs;
switch name
    case 'pure_alpha'
        rng(seed);
        sig = cos(2*pi*10*t) + 0.02*randn(N,1);
    case 'am_alpha'
        rng(seed);
        env = 0.5 + 0.5*cos(2*pi*1.0*t);
        sig = env .* cos(2*pi*10*t) + 0.02*randn(N,1);
    case 'alpha_plus_1f'
        pink = n1f.sim_powerlaw(N, fs, -1, seed); pink = pink / std(pink);
        rng(seed + 50000);
        sig = cos(2*pi*10*t) + pink + 0.02*randn(N,1);
    case 'burst_beta'
        sig = n1f.make_beta_burst('seed', seed);
        sig = sig + 0.02*randn(size(sig));
    case 'burst_plus_1f'
        sig = n1f.make_beta_burst_with_1f('seed', seed);
    otherwise
        error('n1f:make_scenario', 'Unknown scenario "%s".', name);
end
end
