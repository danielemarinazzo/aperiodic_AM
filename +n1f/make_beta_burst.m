function sig = make_beta_burst(varargin)
% make_beta_burst  Bursty beta (Markov on/off), no broadband background.
%
% Switching happens per CYCLE (not per sample), matching NeuroDSP's
% sim_bursty_oscillation behaviour. This prevents the on/off transitions
% from producing artefactual broadband low-frequency power.
p = inputParser; addParameter(p,'seed',42); parse(p,varargin{:});
cfg = n1f.config();
N   = cfg.duration * cfg.fs;
t   = (0:N-1)' / cfg.fs;
rng(p.Results.seed);

freq     = 20;                         % Hz
p_enter  = 0.10;                       % probability of entering burst per cycle
p_leave  = 0.10;                       % probability of leaving burst per cycle
cyc_len  = cfg.fs / freq;             % samples per cycle (12.5 at fs=250)
n_cycles = floor(N / cyc_len);

% Markov chain at cycle level
in_burst   = false;
is_osc     = false(n_cycles, 1);
for c = 1:n_cycles
    if in_burst
        if rand < p_leave, in_burst = false; end
    else
        if rand < p_enter, in_burst = true;  end
    end
    is_osc(c) = in_burst;
end

% Expand cycle-level mask to sample level
mask = false(N, 1);
for c = 1:n_cycles
    i1 = round((c-1)*cyc_len) + 1;
    i2 = min(round(c*cyc_len), N);
    mask(i1:i2) = is_osc(c);
end

sig = mask .* cos(2*pi*freq*t + 2*pi*rand);
end