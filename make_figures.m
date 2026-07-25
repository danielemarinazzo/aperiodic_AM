function make_figures(varargin)
% make_figures  Generate all publication figures.
%
%   make_figures           % all figures
%   make_figures('fig1')   % one figure
%   make_figures('fig1','fig3')
%
%   Prerequisites
%   -------------
%   1. fooof_mat on the MATLAB path (https://github.com/fooof-tools/fooof_mat)
%   2. Signal Processing Toolbox
%   3. Image Processing Toolbox (bwlabel, used in run_fbosc)

here = fileparts(mfilename('fullpath'));
if ~contains(path, here), addpath(here); end

all_figs = {'fig1','fig2','fig3','fig4','fig5','fig6','fig7'};
if nargin > 0
    requested = varargin;
else
    requested = all_figs;
end

for k = 1:numel(requested)
    name = requested{k};
    fprintf('\n=== %s ===\n', name);
    switch name
        case 'fig1', fig1_identifiability();
        case 'fig2', fig2_pipeline();
        case 'fig3', fig3_sweep();
        case 'fig4', fig4_sprint();
        case 'fig5', fig5_erp_sprint();
        case 'fig6', fig6_oddball();
        case 'fig7', fig7_rest();
        otherwise,   fprintf('Unknown figure: %s\n', name);
    end
end
fprintf('\nAll done.\n');
end
