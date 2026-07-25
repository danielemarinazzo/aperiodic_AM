function fig7_rest()
% Figure 7: Scalp EEG and ECoG — rest vs propofol sedation.
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
    'scalp_anes', struct('col',[.53 .75 .90],'ls','--','lw',1.1,'fa',.50,'lbl','Scalp anes'), ...
    'ecog_rest',  struct('col',[.90 .60 .13],'ls','-', 'lw',1.1,'fa',.70,'lbl','ECoG rest'), ...
    'ecog_anes',  struct('col',[.98 .82 .52],'ls','--','lw',1.1,'fa',.50,'lbl','ECoG anes'));
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

%% A: raw signals
ax_A = axes(fig,'Position',[x_wide, rb(1), w_wide, row_h]);
hold(ax_A,'on');
t60 = (0:60*fs_eeg-1)'/fs_eeg;
for k = keys
    r = recs.(k{1}); st = style.(k{1});
    n_plot = min(numel(r.raw), 60*fs_eeg);
    plot(ax_A, t60(1:n_plot), zscore(r.raw(1:n_plot)), ...
         'Color',[st.col st.fa],'LineStyle',st.ls,'LineWidth',.4,'DisplayName',st.lbl);
end
xlim(ax_A,[0 60]); ylabel(ax_A,'Amplitude (z)');
legend(ax_A,'Location','northeast','FontSize',7,'Box','off','TextColor','k','NumColumns',2);
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
set(ax_E,'XTickLabel',[]); ylabel(ax_E,'Power');
title(ax_E,'Full 300 s PSD + specparam fits','FontSize',8,'FontWeight','bold');
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

for bi = 1:n_bands
    flo = band_defs{bi,2}; fhi = band_defs{bi,3};
    x_center = (bi-1) * (group_w + 1.5);
    for ki = 1:n_k
        k = keys{ki}; r = recs.(k); st = style.(k);
        fq   = r.fb.freqs;
        mask = fq >= flo & fq <= fhi;
        ab   = nanmean(r.fb.abundance_segs(:, mask), 2);
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
         'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8);
end
ax_G.XTick = [];
ax_G.XLim  = [min(x_positions(:))-1, max(x_positions(:))+1];
ax_G.YLim  = [0, 0.6];
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
ax_B = axes(fig,'Position',[x_narrow, rb(1), w_narrow, rh3]);
hold(ax_B,'on');
for ki = 1:n_k
    k = keys{ki}; r = recs.(k); st = style.(k);
    violin_kd(ax_B, r.ex(~isnan(r.ex)), ki, 0.35, st.col);
end
ax_B.XTick=1:n_k; ax_B.XTickLabel={'Sr','Sa','Er','Ea'};
ax_B.XTickLabelRotation=30; ylabel(ax_B,'Exponent');
title(ax_B,'Exponent per seg','FontSize',8,'FontWeight','bold');
set_ax(ax_B);

%% D: R² violin
ax_D = axes(fig,'Position',[x_narrow, rb(1)-rh3-rg, w_narrow, rh3]);
hold(ax_D,'on');
for ki = 1:n_k
    k = keys{ki}; r = recs.(k); st = style.(k);
    violin_kd(ax_D, r.r2(~isnan(r.r2)), ki, 0.35, st.col);
end
ax_D.XTick=1:n_k; ax_D.XTickLabel={'Sr','Sa','Er','Ea'};
ax_D.XTickLabelRotation=30; ylabel(ax_D,'R^2');
title(ax_D,'Fit quality R^2','FontSize',8,'FontWeight','bold');
set_ax(ax_D);

%% F: peak histogram
ax_F = axes(fig,'Position',[x_narrow, rb(1)-2*(rh3+rg), w_narrow, rh3]);
hold(ax_F,'on');
all_nk_cell = cellfun(@(k) recs.(k).nk, keys,'UniformOutput',false);
all_nk = vertcat(all_nk_cell{:});
max_pk = max(all_nk,[],'omitnan');
edges  = -0.5:1:max_pk+0.5;
for ki = 1:n_k
    k = keys{ki}; r = recs.(k); st = style.(k);
    histogram(ax_F, r.nk, edges,'Normalization','probability', ...
        'FaceColor',st.col,'FaceAlpha',.5,'EdgeColor','none', ...
        'DisplayName',sprintf('%s m=%.1f',st.lbl,mean(r.nk,'omitnan')));
end
xlim(ax_F,[-0.5 max_pk+1.5]);
xlabel(ax_F,'Peaks / 2 s'); ylabel(ax_F,'Proportion');
legend(ax_F,'FontSize',6,'Box','off','TextColor','k','NumColumns',1);
title(ax_F,'Peak count','FontSize',8,'FontWeight','bold');
set_ax(ax_F);

%% Panel letters
drawnow;
label_pairs = {ax_A,'A'; ax_C,'C'; ax_E,'E'; ax_G,'G'; ...
               ax_B,'B'; ax_D,'D'; ax_F,'F'};
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
fprintf('\n--- Fig5 caption statistics ---\n');
bands = {'alpha',8,12; 'beta',16,30};
for ki = 1:n_k
    k = keys{ki}; r = recs.(k);
    fq = r.fb.freqs;
    fprintf('[%s]  exp=%+.2f (%dpk)  SPRiNT SD=%.2f  seg exp=%.2f+/-%.2f (n=%d)\n', ...
        k, r.sp.exponent, r.sp.n_peaks, nanstd(r.spr.exponent), ...
        nanmean(r.ex), nanstd(r.ex), sum(~isnan(r.ex)));
    for bi = 1:size(bands,1)
        m = fq >= bands{bi,2} & fq <= bands{bi,3};
        vals = nanmedian(r.fb.abundance_segs(:,m), 2);
        fprintf('  %-5s  median=%.3f  IQR=[%.3f %.3f]  n_seg=%d\n', ...
            bands{bi,1}, nanmedian(vals), prctile(vals,25), prctile(vals,75), ...
            sum(~isnan(vals)));
    end
end
fprintf('---\n\n');

save_fig(fig, cfg, 'fig7_rest');
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