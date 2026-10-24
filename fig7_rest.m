function fig7_rest()
% Scalp EEG and ECoG — rest vs propofol sedation.
%
% Three panels (revision 2026-10, reduced from six at the reviewers' request):
%   A  full-recording PSD + specparam fits (exponent in legend)
%   B  exponent per 2 s segment (violins; rank-sum between states)
%   C  fBOSC alpha and beta burst abundance per 5 s segment (rank-sum)
% SPRiNT, R^2 and between-modality tests are computed and printed (reported in
% the text) but not plotted. Previous six-panel version:
% archive/superseded_2026-10/fig5_rest_6panels.m

cfg    = n1f.config();
fs_eeg = 256;
seg    = 2 * fs_eeg;

recs = load_all(cfg.data_dir, fs_eeg, seg);
if isempty(recs), return; end

style = struct( ...
    'scalp_rest', struct('col',[.16 .50 .73],'ls','-', 'lw',1.1,'fa',.70,'lbl','Scalp, rest'), ...
    'scalp_anes', struct('col',[.53 .75 .90],'ls','-', 'lw',1.1,'fa',.50,'lbl','Scalp, sedation'), ...
    'ecog_rest',  struct('col',[.90 .60 .13],'ls','-', 'lw',1.1,'fa',.70,'lbl','ECoG, rest'), ...
    'ecog_anes',  struct('col',[.98 .82 .52],'ls','-', 'lw',1.1,'fa',.50,'lbl','ECoG, sedation'));
keys = {'scalp_rest','scalp_anes','ecog_rest','ecog_anes'};

%% Layout: one row, three panels (A spectra, B exponent per segment, C burst abundance)
fig = figure('Units','centimeters','Position',[0 0 19 7.5],'Color','w');
bm = 0.20; tm = 0.11; h = 1 - bm - tm; lm = 0.075; gap = 0.08;
wA = 0.33; wB = 0.17; wC = 1 - lm - wA - wB - 2*gap - 0.015;
xA = lm; xB = xA + wA + gap; xC = xB + wB + gap;
n_k = numel(keys);

%% A: full-recording PSD + specparam fits
ax_A = axes(fig,'Position',[xA, bm, wA, h]);
hold(ax_A,'on'); set(ax_A,'YScale','log');
for k = keys
    r = recs.(k{1}); st = style.(k{1});
    semilogy(ax_A, r.f, r.p,'Color',[st.col .25],'LineWidth',.4,'HandleVisibility','off');
    semilogy(ax_A, r.sp.freqs, 10.^r.sp.modeled, ...
             'Color',st.col,'LineStyle',st.ls,'LineWidth',1.4, ...
             'DisplayName',sprintf('%s: exp %+.2f',st.lbl,r.sp.exponent));
    semilogy(ax_A, r.sp.freqs, 10.^r.sp.aperiodic, ...
             'Color',st.col,'LineStyle',':','LineWidth',.9,'HandleVisibility','off');
end
all_p_cell = cellfun(@(k) recs.(k).p, keys,'UniformOutput',false);
all_p = vertcat(all_p_cell{:}); all_p = all_p(all_p>0);
ylim(ax_A,[10^floor(log10(min(all_p))) 10^(ceil(log10(max(all_p)))+1)]);   % one decade of headroom for the legend
xlim(ax_A,[1 45]); shade_bands(ax_A);
legend(ax_A,'Location','northeast','FontSize',7,'Box','off','TextColor','k');
ylabel(ax_A,'Power (a.u.)'); xlabel(ax_A,'Frequency (Hz)');
title(ax_A,'Full-recording spectra and fits','FontSize',8,'FontWeight','bold');
set_ax(ax_A);

%% B: exponent per 2 s segment
ax_B = axes(fig,'Position',[xB, bm, wB, h]);
hold(ax_B,'on');
for ki = 1:n_k
    k = keys{ki}; r = recs.(k); st = style.(k);
    violin_kd(ax_B, r.ex(~isnan(r.ex)), ki, 0.35, st.col);
end
ax_B.XTick=1:n_k; ax_B.XTickLabel={'R','S','R','S'};    % rest / sedation
ax_B.XTickLabelRotation=0; ylabel(ax_B,'Exponent');
title(ax_B,'Exponent per 2 s segment','FontSize',8,'FontWeight','bold');
set_ax(ax_B);
yb=ax_B.YLim; pad=0.12*diff(yb); ax_B.YLim=[yb(1) yb(2)+pad];
annot_ranksum(ax_B, recs.(keys{1}).ex, recs.(keys{2}).ex, 1, 2, yb(2));
annot_ranksum(ax_B, recs.(keys{3}).ex, recs.(keys{4}).ex, 3, 4, yb(2));
group_labels(ax_B, {'Scalp','ECoG'}, [1.5 3.5]);

%% C: fBOSC alpha + beta abundance (5 s segments)
ax_C = axes(fig,'Position',[xC, bm, wC, h]);
hold(ax_C,'on');
band_defs = {'Alpha (8-12 Hz)',8,12; 'Beta (16-30 Hz)',16,30};
n_bands = size(band_defs,1);
group_w = 4.0; within = 0.85; bw_v = 0.42;
x_positions = nan(n_bands, n_k);
band_ab     = cell(n_bands, n_k);
for bi = 1:n_bands
    flo = band_defs{bi,2}; fhi = band_defs{bi,3};
    x_center = (bi-1) * (group_w + 1.5);
    for ki = 1:n_k
        k = keys{ki}; r = recs.(k); st = style.(k);
        mask = r.fb.freqs >= flo & r.fb.freqs <= fhi;
        ab   = nanmean(r.fb.abundance_segs(:, mask), 2);
        band_ab{bi,ki} = ab;
        xpos = x_center + (ki - (n_k+1)/2) * (within*group_w/n_k);
        x_positions(bi,ki) = xpos;
        if sum(~isnan(ab)) >= 3 && nanstd(ab) > 1e-6
            violin_kd(ax_C, ab, xpos, bw_v, st.col);
        else
            mv = nanmedian(ab);
            plot(ax_C,[xpos-bw_v*0.85 xpos+bw_v*0.85],[mv mv], ...
                 '-','Color',st.col*.7,'LineWidth',1.5,'HandleVisibility','off');
            scatter(ax_C, xpos+0.05*randn(numel(ab),1), ab, 3, st.col*.8, ...
                    'filled','MarkerFaceAlpha',.3,'HandleVisibility','off');
        end
    end
end
ax_C.XTick = reshape(x_positions',1,[]);
ax_C.XTickLabel = repmat({'R','S','R','S'},1,n_bands);   % scalp R S, ECoG R S, per band
ax_C.XLim  = [min(x_positions(:))-1, max(x_positions(:))+1];
ax_C.YLim  = [0, 1.12];      % headroom: brackets, then band names
ax_C.YTick = 0:0.2:1;
for bi = 1:n_bands
    yb = 0.74;
    annot_ranksum(ax_C, band_ab{bi,1}, band_ab{bi,2}, x_positions(bi,1), x_positions(bi,2), yb);
    annot_ranksum(ax_C, band_ab{bi,3}, band_ab{bi,4}, x_positions(bi,3), x_positions(bi,4), yb+0.08);
    text(ax_C, mean(x_positions(bi,:)), 1.11, band_defs{bi,1}, 'HorizontalAlignment','center', ...
         'VerticalAlignment','top','FontSize',7.5,'Color','k');
end
ylabel(ax_C,'Burst abundance (fBOSC)');
title(ax_C,'Bursts per 5 s segment','FontSize',8,'FontWeight','bold');
set_ax(ax_C);   % modality is given by colour (as in A and B); no group labels: they collide at this width

%% Panel letters
drawnow;
label_pairs = {ax_A,'A'; ax_B,'B'; ax_C,'C'};
for i = 1:size(label_pairs,1)
    pos = label_pairs{i,1}.Position;
    annotation(fig,'textbox',[max(0,pos(1)-0.05), min(0.95,pos(2)+pos(4)+0.02), 0.05, 0.05], ...
        'String',label_pairs{i,2},'FontSize',13,'FontWeight','bold', ...
        'EdgeColor','none','BackgroundColor','none', ...
        'VerticalAlignment','bottom','HorizontalAlignment','left','Color','k');
end

%% Caption statistics
% fBOSC: per segment, abundance averaged over the band's frequencies -- the
% same quantity shown in the violins and used in every test. 'any' = proportion
% of segments with abundance > 0 (the distributions are zero-inflated).
fprintf('\n--- fig5_rest caption statistics ---\n');
bands  = {'alpha',8,12; 'beta',16,30};
bandab = @(r, lo, hi) nanmean(r.fb.abundance_segs(:, r.fb.freqs >= lo & r.fb.freqs <= hi), 2);
for ki = 1:n_k
    k = keys{ki}; r = recs.(k);
    fprintf('[%s]  length=%.1f s  exp=%+.2f (%dpk)  SPRiNT SD=%.2f  seg exp=%.2f+/-%.2f (%d valid of %d segs)\n', ...
        k, numel(r.raw)/fs_eeg, r.sp.exponent, r.sp.n_peaks, nanstd(r.spr.exponent), ...
        nanmean(r.ex), nanstd(r.ex), sum(~isnan(r.ex)), numel(r.ex));
    for bi = 1:size(bands,1)
        vals = bandab(r, bands{bi,2}, bands{bi,3});
        fprintf('  fBOSC %-5s  median=%.3f  IQR=[%.3f %.3f]  mean=%.3f  any=%.2f  n_seg=%d\n', ...
            bands{bi,1}, nanmedian(vals), prctile(vals,25), prctile(vals,75), nanmean(vals), ...
            mean(vals(~isnan(vals)) > 0), sum(~isnan(vals)));
    end
end
fprintf('--- between states, within modality: rank-sum (independent segments) ---\n');
rs = @(a,b) ranksum(a(~isnan(a)), b(~isnan(b)));
for md = {'scalp','ecog'}
    a = recs.([md{1} '_rest']); b = recs.([md{1} '_anes']);
    fprintf('%-5s exponent p=%.3g | R2 p=%.3g', md{1}, rs(a.ex,b.ex), rs(a.r2,b.r2));
    for bi = 1:size(bands,1)
        fprintf(' | fBOSC %s p=%.3g', bands{bi,1}, ...
            rs(bandab(a,bands{bi,2},bands{bi,3}), bandab(b,bands{bi,2},bands{bi,3})));
    end
    fprintf('\n');
end
fprintf('--- between modalities, within state: signed-rank, paired by simultaneous segment ---\n');
for st = {'rest','anes'}
    a = recs.(['scalp_' st{1}]); b = recs.(['ecog_' st{1}]);
    [pe, de] = paired_signrank(a.ex, b.ex);
    [pr, dr] = paired_signrank(a.r2, b.r2);
    fprintf('%-4s exponent scalp-ECoG median diff %+.2f p=%.3g | R2 %+.3f p=%.3g', st{1}, de, pe, dr, pr);
    for bi = 1:size(bands,1)
        [pb, db] = paired_signrank(bandab(a,bands{bi,2},bands{bi,3}), bandab(b,bands{bi,2},bands{bi,3}));
        fprintf(' | fBOSC %s %+.3f p=%.3g', bands{bi,1}, db, pb);
    end
    fprintf('\n');
end
fprintf(['NOTE: scalp (T5) and ECoG channels record overlapping sources at different\n' ...
         '      measurement scale/geometry; differing exponents, peak counts and burst\n' ...
         '      abundance across modality reflect the measurement, not only the neural process.\n']);
fprintf('---\n\n');

save_fig(fig, cfg, 'fig5_rest');
end

%% ── Data loading ─────────────────────────────────────────────────────────
function r = process_rec(raw_vec, fs_eeg, seg)
r.raw = raw_vec(:);
n_s   = floor(numel(r.raw)/seg);
[r.f, r.p] = n1f.compute_psd(r.raw, fs_eeg, fs_eeg*4);
r.sp        = n1f.run_specparam(r.f, r.p, [1 45]);
r.ex = nan(n_s,1); r.r2 = nan(n_s,1); r.nk = zeros(n_s,1);
for i = 1:n_s
    idx = (i-1)*seg+1:i*seg;
    try
        [fc,pc] = n1f.compute_psd(r.raw(idx), fs_eeg, seg);
        sc = n1f.run_specparam(fc, pc, [1 45]);
        r.ex(i)=sc.exponent; r.r2(i)=sc.r_squared; r.nk(i)=sc.n_peaks;
    catch, end
end
n60   = min(numel(r.raw), 60*fs_eeg);
r.spr = n1f.run_sprint(r.raw(1:n60), fs_eeg);
seg_fb  = 5 * fs_eeg;
n_sfb   = floor(numel(r.raw) / seg_fb);
fq_fb   = (4:38)';
ab_segs = nan(n_sfb, numel(fq_fb));
for i = 1:n_sfb
    idx = (i-1)*seg_fb+1 : i*seg_fb;
    try
        fb_i = n1f.run_fbosc(r.raw(idx), fs_eeg, 'freqs',fq_fb,'freq_range',[4 38]);
        ab_segs(i,:) = fb_i.abundance(:)';
    catch, end
end
r.fb.freqs          = fq_fb;
r.fb.abundance      = nanmedian(ab_segs, 1)';
r.fb.abundance_segs = ab_segs;
end

function recs = load_all(data_dir, fs_eeg, seg)
recs  = [];
files = struct('scalp_rest','scalp_rest.txt','scalp_anes','scalp_anes.txt', ...
               'ecog_rest','ECoG_rest.txt','ecog_anes','ECoG_anes.txt');
for k = fieldnames(files)'
    fpath = fullfile(data_dir, files.(k{1}));
    if ~exist(fpath,'file')
        fprintf('  File not found: %s\n', fpath); return;
    end
    fprintf('  Loading %s...\n', k{1});
    raw = readmatrix(fpath); raw = raw(:);
    recs.(k{1}) = process_rec(raw, fs_eeg, seg);
end
end

%% ── Helpers ──────────────────────────────────────────────────────────────
function shade_bands(ax)
hold(ax,'on');
yl = ax.YLim;
defs = {[4 8],[.20 .60 .80]; [8 12],[.80 .20 .20]; [16 30],[.20 .70 .30]};
for i = 1:3
    x = defs{i,1}; c = defs{i,2};
    fill(ax,[x(1) x(2) x(2) x(1)],[yl(1) yl(1) yl(2) yl(2)], ...
         c,'FaceAlpha',.07,'EdgeColor','none','HandleVisibility','off');
end
end

function set_ax(ax)
set(ax,'Box','off','TickDir','out','FontSize',8,'LineWidth',.7, ...
    'Color','w','XColor','k','YColor','k');
ax.Title.Color = 'k';
ax.Title.Units = 'normalized';
ax.Title.HorizontalAlignment = 'left';
ax.Title.Position(1) = 0;
end

function [p, dmed] = paired_signrank(a, b)
% Signed-rank on segment-wise differences a - b. Scalp and ECoG are recorded
% simultaneously, so segment i covers the same time window in both.
n = min(numel(a), numel(b)); d = a(1:n) - b(1:n); d = d(~isnan(d));
dmed = median(d);
if isempty(d) || all(d == 0), p = NaN; else, p = signrank(d); end
end

function group_labels(ax, names, xs)
% Group names centred under pairs of ticks, below the tick labels.
yl = ax.YLim; y = yl(1) - 0.065*diff(yl);
for i = 1:numel(names)
    text(ax, xs(i), y, names{i}, 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize', 8, 'Color', 'k', 'Clipping', 'off');
end
end

function annot_ranksum(ax, a, b, x1, x2, ytop)
% Wilcoxon rank-sum between a and b; bracket + p above the pair.
a=a(~isnan(a)); b=b(~isnan(b));
if numel(a)<3 || numel(b)<3, return; end
try, p=ranksum(a,b); catch, return; end
xm=(x1+x2)/2;
plot(ax,[x1 x1 x2 x2],[ytop ytop+0.02*diff(ax.YLim) ytop+0.02*diff(ax.YLim) ytop], ...
     '-','Color',[.3 .3 .3],'LineWidth',.6,'HandleVisibility','off');
if     p<0.001, s='***'; elseif p<0.01, s='**'; elseif p<0.05, s='*'; else, s='n.s.'; end
text(ax,xm,ytop+0.03*diff(ax.YLim),s,'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom','FontSize',7,'Color','k');
end

function violin_kd(ax, data, x, hw, col)
if nargin < 4 || isempty(hw), hw = 0.35; end
if isempty(data) || numel(data) < 3, return; end
[f,xi] = ksdensity(data,'NumPoints',100);
f = f(:)'/max(f)*hw; xi = xi(:)';
fill(ax,[x+f fliplr(x-f)],[xi fliplr(xi)],col, ...
     'FaceAlpha',.50,'EdgeColor','none','HandleVisibility','off');
plot(ax,[x-hw*0.85 x+hw*0.85],[median(data) median(data)],'-', ...
     'Color',col*.7,'LineWidth',1.5,'HandleVisibility','off');
scatter(ax,x+hw*0.17*randn(size(data)),data,2,col*.8,'filled', ...
        'MarkerFaceAlpha',.25,'HandleVisibility','off');
end

function save_fig(fig, cfg, name)
fname = fullfile(cfg.out_dir,[name '.png']);
exportgraphics(fig,fname,'Resolution',150);
fprintf('  Saved: %s\n',fname); close(fig);
end