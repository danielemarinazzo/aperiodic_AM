function sig = make_beta_burst_with_1f(varargin)
% make_beta_burst_with_1f  Bursty beta + stationary 1/f background.
p = inputParser; addParameter(p,'seed',42); parse(p,varargin{:});
seed = p.Results.seed;
cfg  = n1f.config();
N    = cfg.duration * cfg.fs;
burst = n1f.make_beta_burst('seed', seed);
pink  = n1f.sim_powerlaw(N, cfg.fs, cfg.exp_1f, seed+2000);
sig   = burst + pink/std(pink);
end