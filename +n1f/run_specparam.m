function sp = run_specparam(freqs, psd, freq_range, varargin)
% run_specparam  Thin wrapper around fooof_mat returning a uniform struct.
%
%   Requires fooof_mat on the MATLAB path:
%     https://github.com/fooof-tools/fooof_mat
%
%   sp = n1f.run_specparam(freqs, psd, [1 45])
%   sp = n1f.run_specparam(freqs, psd, [1 45], 'peak_width_limits',[0.5 12], ...
%                           'max_n_peaks',8, 'min_peak_height',0.05)
%
%   Output struct fields
%   --------------------
%   sp.exponent       aperiodic exponent (positive = steeper slope)
%   sp.offset         aperiodic offset (log10 units)
%   sp.knee           knee frequency (NaN for fixed model)
%   sp.r_squared      goodness of fit
%   sp.n_peaks        number of detected peaks
%   sp.peaks          struct array with fields cf, pw, bw
%   sp.freqs          frequency vector used for fitting
%   sp.modeled        log10 full model
%   sp.aperiodic      log10 aperiodic component

if nargin < 3 || isempty(freq_range)
    freq_range = [1 45];
end

p = inputParser;
addParameter(p, 'peak_width_limits', [0.5 12]);
addParameter(p, 'max_n_peaks',       inf);
addParameter(p, 'min_peak_height',   0.05);
addParameter(p, 'peak_threshold',    2.0);
addParameter(p, 'aperiodic_mode',    'fixed');
parse(p, varargin{:});
s = p.Results;

settings = struct( ...
    'peak_width_limits', s.peak_width_limits, ...
    'max_n_peaks',       s.max_n_peaks, ...
    'min_peak_height',   s.min_peak_height, ...
    'peak_threshold',    s.peak_threshold, ...
    'verbose',           0, ...
    'aperiodic_mode',    s.aperiodic_mode);

% fooof rejects freq == 0 (log10(0) = -inf) and any non-positive power.
% Strip those bins before handing over.
freqs = freqs(:); psd = psd(:);
keep  = freqs > 0 & isfinite(psd) & psd > 0;
freqs = freqs(keep);
psd   = psd(keep);

% fooof_mat expects freqs, psd as row vectors of doubles
fr = fooof(double(freqs(:)'), double(psd(:)'), double(freq_range), settings, true);

% ── Unpack ─────────────────────────────────────────────────────────────
ap = fr.aperiodic_params;
sp.offset = ap(1);
if strcmp(s.aperiodic_mode,'knee') && length(ap) == 3
    sp.knee     = ap(2);
    sp.exponent = ap(3);
else
    sp.knee     = NaN;
    sp.exponent = ap(end);
end
sp.r_squared = fr.r_squared;

% Peaks — defensive against empty / odd shapes
pk = fr.peak_params;
if isempty(pk)
    sp.n_peaks = 0;
elseif size(pk,2) == 3
    sp.n_peaks = size(pk,1);
elseif numel(pk) == 3                  % single peak as a 1×3 row vector edge case
    sp.n_peaks = 1;
    pk = reshape(pk, 1, 3);
else
    sp.n_peaks = 0;
end

sp.peaks = struct('cf',{},'pw',{},'bw',{});
for k = 1:sp.n_peaks
    sp.peaks(k).cf = pk(k,1);
    sp.peaks(k).pw = pk(k,2);
    sp.peaks(k).bw = pk(k,3);
end

sp.freqs     = fr.freqs(:);
sp.modeled   = fr.fooofed_spectrum(:);
sp.aperiodic = fr.ap_fit(:);
end