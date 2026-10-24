function make_figures(varargin)
% make_figures  Generate all figures and tables, by script name.
%
%   make_figures                          % everything, in the order below
%   make_figures('fig_concept')           % one script
%   make_figures('fig_toy','figS_knee')   % several
%
%   Scripts are called by name; figure numbers live only in the manuscript.
%   fig_modulation_f20 runs fig_modulation(20) (supplementary 20 Hz carrier).
%
%   Prerequisites
%   -------------
%   1. fooof_mat on the MATLAB path (https://github.com/fooof-tools/fooof_mat)
%   2. Signal Processing Toolbox
%   3. Statistics and Machine Learning Toolbox (ranksum, signrank, poissrnd, corr)
%   4. Image Processing Toolbox (bwlabel, used in run_fbosc)

here = fileparts(mfilename('fullpath'));
if ~contains(path, here), addpath(here); end

% Main figures in manuscript order, then supplementary (fig1_identifiability is
% supplementary since the 2026-10 revision).
all_scripts = {'fig_concept', 'fig2_pipeline', ...
               'table1_stats', 'fig_toy', 'fig_modulation', ...
               'fig3_erp_sprint', 'fig4_oddball', 'fig5_rest', 'fig_timescale', ...
               'fig1_identifiability', 'fig_modulation_f20', 'figS_duty', 'figS_knee', ...
               'check_lorentz'};
if nargin > 0
    requested = varargin;
else
    requested = all_scripts;
end

for k = 1:numel(requested)
    name = requested{k};
    fprintf('\n=== %s ===\n', name);
    switch name
        case 'fig_modulation_f20', fig_modulation(20);
        otherwise
            if any(strcmp(name, all_scripts))
                feval(name);
            else
                fprintf('Unknown script: %s\n', name);
            end
    end
end
fprintf('\nAll done.\n');
end
