function fig5_rest()
% Scalp EEG and ECoG — rest vs propofol sedation.
%
% Panel reading order (left-to-right, top-to-bottom):
%   A  raw signals (time)         B  exponent violin
%   C  SPRiNT exponent (time)     D  R² violin
%   E  full PSD (frequency)       F  peak count histogram
%   G  fBOSC alpha+beta abundance (5 s segs)
%
% Left column (A,C,E,G): wide. Right column (B,D,F): narrow.
% A,C share time x-axis.

cfg    = n1f.config();
fs_eeg = 256;
seg    = 2 * fs_eeg;

recs = load_all(cfg.data_dir, fs_eeg, seg);
if isempty(recs), return; end

style = struct( ...
    'scalp_rest', struct('col',[.16 .50 .73],'ls','-', 'lw',1.1,'fa',.70,'lbl','Scalp rest'), ...
    'scalp_anes', struct('col',[.53 .75 .90],'ls','-', 'lw',1.1,'fa',.50,'lbl','Scalp anes'), ...
    'ecog_rest',  struct('col',[.90 .60 .13],'ls','-', 'lw',1.1,'fa',.70,'lbl','ECoG rest'), ...
    'ecog_anes',  struct('col',[.98 .82 .52],'ls','-', 'lw',1.1,'fa',.50,'lbl','ECoG anes'));
keys = {'scalp_rest','scalp_anes','ecog_rest','ecog_anes'};

%% Layout
fig = figure('Units','centimeters','Position',[0 0 26 24],'Color','w');

lm = 0.09; mg = 0.05; rm = 0.02; tm = 0.10; bm = 0.07; rg = 0.04;
w_wide   = 0.56;
w_narrow = 1 - lm - mg - w_wide - rm;
n_rows   = 4;
row_h    = (1 - tm - bm - (n_rows-1)*rg) / n_rows;
rb       = @(r) 1 - tm - (r-1)*(row_h+rg) - row_h;
x_wide   = lm;
x_narrow = lm + w_wide + mg;
rh3      = (3*(row_h+rg) - rg) / 3;
rh2      = (4*(row_h+rg) - rg) / 2;   % right col: 2 panels over full height

%% A: raw signals
ax_A = axes(fig,'Position',[x_wide, rb(1), w_wide, row_h]);
hold(ax_A,'on');
t60 = (0:60*fs_eeg-1)'/fs_eeg;
off_step = 6;                        % vertical offset between stacked traces (z units)
y_ticks = zeros(1,numel(keys)); y_lbls = cell(1,numel(keys));
for ki = 1:numel(keys)
    k = keys{ki}; r = recs.(k); st = style.(k);
    n_plot = min(numel(r.raw), 60*fs_eeg);
    yoff = (numel(keys)-ki)*off_step;         % top-to-bottom in key order
    plot(ax_A, t60(1:n_plot), zscore(r.raw(1:n_plot))+yoff, ...
         'Color',st.col,'LineStyle','-','LineWidth',.4,'DisplayName',st.lbl);
    y_ticks(ki)=yoff; y_lbls{ki}=st.lbl;
end
xlim(ax_A,[0 60]);
[y_ticks,ord]=sort(y_ticks); y_lbls=y_lbls(ord);
set(ax_A,'YTick',y_ticks,'YTickLabel',y_lbls);
ylabel(ax_A,'Amplitude (z, offset)');
set(ax_A,'XTickLabel',[]);
title(ax_A,'Raw signals, first 60 s','FontSize',8,'FontWeight','bold');
set_ax(ax_A);

%% C: SPRiNT exponent
ax_C = axes(fig,'Position',[x_wide, rb(2), w_wide, row_h]);
hold(ax_C,'on');
for k = keys
    r = recs.(k{1}); st = style.(k{1});
    v = ~isnan(r.spr.exponent);
    plot(ax_C, r.spr.times(v), r.spr.exponent(v), ...
         'Color',st.col,'LineStyle',st.ls,'LineWidth',st.lw, ...
         'DisplayName',sprintf('%s  SD=%.2f',st.lbl,nanstd(r.spr.exponent)));
end
yline(ax_C,0,':k','LineWidth',.4,'HandleVisibility','off');
xlim(ax_C,[0 60]); ylabel(ax_C,'Exponent'); xlabel(ax_C,'Time (s)');
legend(ax_C,'Location','northeast','FontSize',7,'Box','off','TextColor','k','NumColumns',2);
title(ax_C,'SPRiNT exponent, first 60 s','FontSize',8,'FontWeight','bold');
set_ax(ax_C);

drawnow;
ax_A.Position(1) = ax_C.Position(1);
ax_A.Position(3) = ax_C.Position(3);
linkaxes([ax_A ax_C],'x');

%% E: full PSD
ax_E = axes(fig,'Position',[x_wide, rb(3), w_wide, row_h]);
hold(ax_E,'on'); set(ax_E,'YScale','log');
for k = keys
    r = recs.(k{1}); st = style.(k{1});
    semilogy(ax_E, r.f, r.p,'Color',[st.col .20],'LineWidth',.4,'HandleVisibility','off');
    semilogy(ax_E, r.sp.freqs, 10.^r.sp.modeled, ...
             'Color',st.col,'LineStyle',st.ls,'LineWidth',1.4, ...
             'DisplayName',sprintf('%s  exp=%+.2f (%dpk)',st.lbl,r.sp.exponent,r.sp.n_peaks));
    semilogy(ax_E, r.sp.freqs, 10.^r.sp.aperiodic, ...
             'Color',st.col,'LineStyle',':','LineWidth',.8,'HandleVisibility','off');
end
all_p_cell = cellfun(@(k) recs.(k).p, keys,'UniformOutput',false);
all_p = vertcat(all_p_cell{:}); all_p = all_p(all_p>0);
ylim(ax_E,[10^floor(log10(min(all_p))) 10^ceil(log10(max(all_p)))]);
xlim(ax_E,[1 45]); shade_bands(ax_E);
legend(ax_E,'Location','northeast','FontSize',7,'Box','off','TextColor','k','NumColumns',2);
ylabel(ax_E,'Power (a.u.)'); xlabel(ax_E,'Frequency (Hz)');
title(ax_E,'Power spectral density (full 300 s) + specparam fits','FontSize',8,'FontWeight','bold');
set_ax(ax_E);

%% G: fBOSC alpha + beta abundance (5 s segments, no theta)
ax_G = axes(fig,'Position',[x_wide, rb(4), w_wide, row_h]);
hold(ax_G,'on');
band_defs = {'Alpha (8-12)',8,12; 'Beta (16-30)',16,30};
n_bands = size(band_defs,1);
n_k     = numel(keys);
group_w = 4.0;
within  = 0.85;
bw_v    = 0.42;
x_positions = nan(n_bands, n_k);
band_ab     = cell(n_bands, n_k);

for bi = 1:n_bands
    flo = band_defs{bi,2}; fhi = band_defs{bi,3};
    x_center = (bi-1) * (group_w + 1.5);
    for ki = 1:n_k
        k = keys{ki}; r = recs.(k); st = style.(k);
        fq   = r.fb.freqs;
        mask = fq >= flo & fq <= fhi;
        ab   = nanmean(r.fb.abundance_segs(:, mask), 2);
        band_ab{bi,ki} = ab;
        off  = (ki - (n_k+1)/2) * (within*group_w/n_k);
        xpos = x_center + off;
        x_positions(bi,ki) = xpos;
        if sum(~isnan(ab)) >= 3 && nanstd(ab) > 1e-6
            violin_kd(ax_G, ab, xpos, bw_v, st.col);
        else
            mv = nanmedian(ab);
            plot(ax_G,[xpos-bw_v*0.85 xpos+bw_v*0.85],[mv mv], ...
                 '-','Color',st.col*.7,'LineWidth',1.5,'HandleVisibility','off');
            scatter(ax_G, xpos+0.05*randn(numel(ab),1), ab, 3, st.col*.8, ...
                    'filled','MarkerFaceAlpha',.3,'HandleVisibility','off');
        end
    end
end

for bi = 1:n_bands
    text(ax_G, mean(x_positions(bi,:)), -0.045, band_defs{bi,1}, ...
         'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color','k');
end
ax_G.XTick = [];
ax_G.XLim  = [min(x_positions(:))-1, max(x_positions(:))+1];
ax_G.YLim  = [0, 0.6];
for bi = 1:n_bands
    yb = 0.46;  % annotate below the 0.6 ceiling
    annot_ranksum(ax_G, band_ab{bi,1}, band_ab{bi,2}, x_positions(bi,1), x_positions(bi,2), yb);
    annot_ranksum(ax_G, band_ab{bi,3}, band_ab{bi,4}, x_positions(bi,3), x_positions(bi,4), yb+0.055);
end
ylabel(ax_G,'Abundance'); xlabel(ax_G,'Band');
for ki = 1:n_k
    st = style.(keys{ki});
    plot(ax_G, NaN, NaN, 's','MarkerFaceColor',st.col,'MarkerEdgeColor','none', ...
         'MarkerSize',8,'DisplayName',st.lbl);
end
legend(ax_G,'Location','northeast','FontSize',7,'Box','off','TextColor','k','NumColumns',2);
title(ax_G,'fBOSC burst abundance — alpha and beta (5 s segs)','FontSize',8,'FontWeight','bold');
set_ax(ax_G);

%% B: exponent violin
ax_B = axes(fig,'Position',[x_narrow, rb(1), w_narrow, rh2]);
hold(ax_B,'on');
for ki = 1:n_k
    k = keys{ki}; r = recs.(k); st = style.(k);
    violin_kd(ax_B, r.ex(~isnan(r.ex)), ki, 0.35, st.col);
end
ax_B.XTick=1:n_k; ax_B.XTickLabel={'Sr','Sa','Er','Ea'};
ax_B.XTickLabelRotation=30; ylabel(ax_B,'Exponent');
title(ax_B,'Exponent per seg','FontSize',8,'FontWeight','bold');
set_ax(ax_B);
yb=ax_B.YLim; pad=0.10*diff(yb); ax_B.YLim=[yb(1) yb(2)+pad];
annot_ranksum(ax_B, recs.(keys{1}).ex, recs.(keys{2}).ex, 1, 2, yb(2));
annot_ranksum(ax_B, recs.(keys{3}).ex, recs.(keys{4}).ex, 3, 4, yb(2));

%% D: R² violin
ax_D = axes(fig,'Position',[x_narrow, rb(1)-rh2-rg, w_narrow, rh2]);
hold(ax_D,'on');
for ki = 1:n_k
    k = keys{ki}; r = recs.(k); st = style.(k);
    violin_kd(ax_D, r.r2(~isnan(r.r2)), ki, 0.35, st.col);
end
ax_D.XTick=1:n_k; ax_D.XTickLabel={'Sr','Sa','Er','Ea'};
ax_D.XTickLabelRotation=30; ylabel(ax_D,'R^2');
title(ax_D,'Fit quality R^2','FontSize',8,'FontWeight','bold');
set_ax(ax_D);
yd=ax_D.YLim; pad=0.10*diff(yd); ax_D.YLim=[yd(1) yd(2)+pad];
annot_ranksum(ax_D, recs.(keys{1}).r2, recs.(keys{2}).r2, 1, 2, yd(2));
annot_ranksum(ax_D, recs.(keys{3}).r2, recs.(keys{4}).r2, 3, 4, yd(2));

%% Panel letters
drawnow;
label_pairs = {ax_A,'A'; ax_C,'C'; ax_E,'E'; ax_G,'F'; ...
               ax_B,'B'; ax_D,'D'};
for i = 1:size(label_pairs,1)
    ax  = label_pairs{i,1};
    ltr = label_pairs{i,2};
    pos = ax.Position;
    box_left = max(0, pos(1) - 0.04);
    box_bot  = min(0.97, pos(2) + pos(4) + 0.005);
    annotation(fig,'textbox',[box_left, box_bot, 0.05, 0.03], ...
        'String',ltr,'FontSize',13,'FontWeight','bold', ...
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