function cfg = config()
% config  Global simulation parameters and paths.
%
% Expected directory layout:
%   <project_root>/
%       data/           (cfg.data_dir)  -- empirical data files
%       matlab_code/    (this directory)
%       python_code/
%
% Output figures are written to <project_root>/figures_matlab/ by default.

cfg.fs         = 250;
cfg.duration   = 10;
cfg.seed       = 42;
cfg.alpha_freq = 10;
cfg.beta_freq  = 20;
cfg.noise_amp  = 0.3;
cfg.am_freq    = 1.0;
cfg.exp_1f     = -1;
cfg.freq_range = [1 45];
cfg.bands      = struct('theta',[4 8],'alpha',[8 12],'beta',[16 30]);

% Paths (resolved relative to this file's location)
here           = fileparts(mfilename('fullpath'));    % .../matlab_code/+n1f
matlab_root    = fileparts(here);                     % .../matlab_code
project_root   = fileparts(matlab_root);              % .../<project_root>

cfg.data_dir   = fullfile(project_root, 'data');
cfg.out_dir    = fullfile(project_root, 'figures_matlab');

if ~exist(cfg.out_dir,'dir'), mkdir(cfg.out_dir); end
end
