function fig2_pipeline()
% Figure 2 — Multi-method pipeline on Condition 15 (beta burst + 1/f).
%
% Manual axes layout so time panels (A, C, E) are left-edge aligned and
% frequency panels (B, D) are left-edge aligned in the left column.
%
%   Left column  (x=freq): B specparam, D fBOSC
%   Right column (x=time): A signal,    C envelope,  E SPRiNT

cfg = n1f.config();
fs  = cfg.fs;
rng(cfg.seed);

sig = n1f.make_beta_burst_with_1f('seed', cfg.seed);
t   = (0:length(sig)-1)' / fs;

[freqs, psd] = n1f.compute_psd(sig, fs);
sp  = n1f.run_specparam(freqs, psd, cfg.freq_range);
spr = n1f.run_sprint(sig, fs);
fb  = n1f.run_fbosc(sig, fs, 'freqs',(4:44)', 'freq_range',[4 44]);

[b1,a1] = butter(4,[16 30]/(fs/2),'bandpass');
[b2,a2] = butter(4,3/(fs/2),'low');
env_beta = filtfilt(b2,a2, abs(hilbert(filtfilt(b1,a1,sig))));
cv = std(env_beta)/(mean(env_beta)+eps);

cn = [.17 .24 .31]; cg = [.15 .68 .38]; cp = [.56 .27 .68]; cr = [.90 .27 .21];

fig = figure('Units','centimeters','Position',[0 0 18 20],'Color','w');

% ── Manual layout in normalised figure units ──────────────────────────────
% [left  bottom  width  height]
lm = 0.08;   % left margin
mg = 0.10;   % gap between columns
rm = 0.03;   % right margin
tm = 0.04;   % top margin
bm = 0.07;   % bottom margin
rg = 0.05;   % row gap

col_w  = (1 - lm - mg - rm) / 2;   % width of each column
row_h  = (1 - tm - bm - 2*rg) / 3; % height of each row

TEMP_left_PLACEHOLDER = lm;
left_x  = lm;
right_x = lm + col_w + mg;

% row bottom edges (row 1 = top)
row_b = @(r) 1 - tm - r*row_h - (r-1)*rg;  % bottom of row r

% Left column: B (row 1), D (row 2) — frequency domain
pos_A = [left_x,  row_b(1), col_w, row_h];
pos_C = [left_x,  row_b(2), col_w, row_h];
pos_E = [left_x,  row_b(3), col_w, row_h];

% Right column: B (row 1), D (row 2) — frequency domain
pos_B = [right_x, row_b(1), col_w, row_h];
pos_D = [right_x, row_b(2), col_w, row_h];

% ── Panel B: specparam (freq) ─────────────────────────────────────────────
ax_B = axes(fig,'Position',pos_B); hold(ax_B,'on');
semilogy(ax_B,freqs,psd,'Color',[.75 .75 .75],'LineWidth',.6);
semilogy(ax_B,sp.freqs,10.^sp.modeled,  'Color',cr,'LineWidth',1.5);
semilogy(ax_B,sp.freqs,10.^sp.aperiodic,'Color',cr,'LineWidth',1.0,'LineStyle','--');
shade_bands(ax_B); xlim(ax_B,[1 45]);
legend(ax_B,sprintf('exp = %+.2f',sp.exponent),'Location','northeast','FontSize',7,'Box','off','TextColor','k');
set(ax_B,'XTickLabel',[]);
ylabel(ax_B,'Power');
title(ax_B,'specparam','FontSize',8,'FontWeight','bold');
set_ax(ax_B); fix_title(ax_B); plabel(ax_B,'B');

% ── Panel D: fBOSC (freq) ─────────────────────────────────────────────────
ax_D = axes(fig,'Position',pos_D); hold(ax_D,'on');
for k = 1:length(fb.freqs)
    if fb.abundance(k) > 0.10, fc = cr; else, fc = [.75 .75 .75]; end
    bar(ax_D,fb.freqs(k),fb.abundance(k),.8,'FaceColor',fc,'EdgeColor','none');
end
yline(ax_D,.10,'--','Color',[.6 .6 .6],'LineWidth',.7);
shade_bands(ax_D); xlim(ax_D,[1 45]);
xlabel(ax_D,'Frequency (Hz)'); ylabel(ax_D,'Abundance');
title(ax_D,'fBOSC burst abundance','FontSize',8,'FontWeight','bold');
set_ax(ax_D); fix_title(ax_D); plabel(ax_D,'D');

% ── Panel A: raw signal (time) ────────────────────────────────────────────
ax_A = axes(fig,'Position',pos_A); hold(ax_A,'on');
plot(ax_A,t,sig,'Color',cn,'LineWidth',.5);
xlim(ax_A,[0 10]);
set(ax_A,'XTickLabel',[]);
ylabel(ax_A,'Amplitude');
title(ax_A,'Beta burst + 1/f','FontSize',8,'FontWeight','bold');
set_ax(ax_A); fix_title(ax_A); plabel(ax_A,'A');

% ── Panel C: beta envelope (time) ─────────────────────────────────────────
ax_C = axes(fig,'Position',pos_C); hold(ax_C,'on');
plot(ax_C,t,env_beta,'Color',cg,'LineWidth',.9);
xlim(ax_C,[0 10]);
text(ax_C,.57,.93,sprintf('CV = %.2f',cv),'Units','normalized', ...
     'HorizontalAlignment','right','FontSize',8,'Color',cg);
set(ax_C,'XTickLabel',[]);
ylabel(ax_C,'Envelope');
title(ax_C,'Beta amplitude envelope','FontSize',8,'FontWeight','bold');
set_ax(ax_C); fix_title(ax_C); plabel(ax_C,'C');

% ── Panel E: SPRiNT (time) ────────────────────────────────────────────────
ax_E = axes(fig,'Position',pos_E); hold(ax_E,'on');
v = ~isnan(spr.exponent);
plot(ax_E,spr.times(v),spr.exponent(v),'o-','MarkerSize',4,'LineWidth',1.2,'Color',cp);
area(ax_E,spr.times(v),max(spr.exponent(v),0),'FaceColor',cp,'FaceAlpha',.10,'EdgeColor','none');
yline(ax_E,0,':k','LineWidth',.5);
xlim(ax_E,[0 10]);
xlabel(ax_E,'Time (s)'); ylabel(ax_E,'Exponent');
title(ax_E,sprintf('SPRiNT exponent  (SD = %.2f)',nanstd(spr.exponent)), ...
      'FontSize',8,'FontWeight','bold');
set_ax(ax_E); fix_title(ax_E); plabel(ax_E,'E');

% ── Link axes ──────────────────────────────────────────────────────────────
linkaxes([ax_A ax_C ax_E], 'x');
linkaxes([ax_B ax_D],      'x');

save_fig(fig, cfg, 'fig2_pipeline');
end


%% ── Local helpers ──────────────────────────────────────────────────────────
function shade_bands(ax)
hold(ax,'on');
defs = {[4 8],[.20 .60 .80];[8 12],[.80 .20 .20];[16 30],[.20 .70 .30]};
yl = get(ax,'YLim');
for i = 1:3, x=defs{i,1}; c=defs{i,2};
    fill(ax,[x(1) x(2) x(2) x(1)],[yl(1) yl(1) yl(2) yl(2)], ...
         c,'FaceAlpha',.07,'EdgeColor','none','HandleVisibility','off'); end
end
function set_ax(ax)
set(ax,'Box','off','TickDir','out','FontSize',8,'LineWidth',.7, ...
    'Color','w','XColor','k','YColor','k');
ax.Title.Color = 'k';
end
function fix_title(ax)
ax.Title.Units            = 'normalized';
ax.Title.HorizontalAlignment = 'left';
ax.Title.Position(1)      = 0;
end
function plabel(ax,letter)
text(ax,-0.11,1.05,letter,'Units','normalized','FontSize',12,'FontWeight','bold', ...
     'VerticalAlignment','bottom','Color','k');
end
function save_fig(fig,cfg,name)
fname=fullfile(cfg.out_dir,[name '.png']);
exportgraphics(fig,fname,'Resolution',150); fprintf('  Saved: %s\n',fname); close(fig);
end