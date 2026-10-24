function fig2_pipeline()
% Multi-method pipeline on scenario E (bursty beta + 1/f).
%
% Manual axes layout:
%   Left column  (x=time): A signal, C envelope, E SPRiNT
%   Right column: B specparam, D fBOSC (x=freq); F zoom on one burst (x=time)

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

fig = figure('Units','centimeters','Position',[0 0 18 22],'Color','w');

% ── Manual layout in normalised figure units ──────────────────────────────
% [left  bottom  width  height]
lm = 0.08;   % left margin
mg = 0.10;   % gap between columns
rm = 0.03;   % right margin
tm = 0.04;   % top margin
bm = 0.07;   % bottom margin
rg = 0.085;  % row gap: room for x labels + next title (row height unchanged, ~5.3 cm)

col_w  = (1 - lm - mg - rm) / 2;   % width of each column
row_h  = (1 - tm - bm - 2*rg) / 3; % height of each row

TEMP_left_PLACEHOLDER = lm;
left_x  = lm;
right_x = lm + col_w + mg;

% row bottom edges (row 1 = top)
row_b = @(r) 1 - tm - r*row_h - (r-1)*rg;  % bottom of row r

% Left column: A, C, E — time domain
pos_A = [left_x,  row_b(1), col_w, row_h];
pos_C = [left_x,  row_b(2), col_w, row_h];
pos_E = [left_x,  row_b(3), col_w, row_h];

% Right column: B, D — frequency domain; F — burst zoom (time)
pos_B = [right_x, row_b(1), col_w, row_h];
pos_D = [right_x, row_b(2), col_w, row_h];
pos_F = [right_x, row_b(3), col_w, row_h];

% ── Panel B: specparam (freq) ─────────────────────────────────────────────
ax_B = axes(fig,'Position',pos_B); hold(ax_B,'on');
plot(ax_B,freqs,psd,'Color',[.75 .75 .75],'LineWidth',.6);
hB = plot(ax_B,sp.freqs,10.^sp.modeled,'Color',cr,'LineWidth',1.5);
plot(ax_B,sp.freqs,10.^sp.aperiodic,'Color',cr,'LineWidth',1.0,'LineStyle','--');
set(ax_B,'YScale','log');          % hold was on before plotting, so set log scale explicitly
xlim(ax_B,[1 45]);
pB = psd(freqs >= 1 & freqs <= 45);
ylim(ax_B,[10^floor(log10(min(pB))) 10^ceil(log10(max(pB)))]);
shade_bands(ax_B);                 % after ylim: shading spans the final y-range
legend(ax_B,hB,sprintf('exp = %+.2f',sp.exponent),'Location','northeast','FontSize',7,'Box','off','TextColor','k');
xlabel(ax_B,'Frequency (Hz)');
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
xlabel(ax_A,'Time (s)');
ylabel(ax_A,'Amplitude');
title(ax_A,'Beta burst + 1/f','FontSize',8,'FontWeight','bold');
set_ax(ax_A); fix_title(ax_A); plabel(ax_A,'A');
% Zoom window for panel F: 0.8 s around the largest beta-envelope value
[~, i_pk] = max(env_beta);
t_c  = min(max(t(i_pk), 0.4), t(end) - 0.4);
w_in = t >= t_c - 0.4 & t <= t_c + 0.4;
yl_A = ax_A.YLim;
rectangle(ax_A,'Position',[t_c-0.4, yl_A(1), 0.8, diff(yl_A)],'EdgeColor',[.89 .10 .11],'LineWidth',.8);

% ── Panel C: beta envelope (time) ─────────────────────────────────────────
ax_C = axes(fig,'Position',pos_C); hold(ax_C,'on');
plot(ax_C,t,env_beta,'Color',cg,'LineWidth',.9);
xlim(ax_C,[0 10]);
text(ax_C,.57,.93,sprintf('CV = %.2f',cv),'Units','normalized', ...
     'HorizontalAlignment','right','FontSize',8,'Color',cg);
xlabel(ax_C,'Time (s)');
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

% ── Panel F: zoom on one burst (window = red box in A) ─────────────────────
ax_F = axes(fig,'Position',pos_F); hold(ax_F,'on');
beta_sig = filtfilt(b1, a1, sig);                  % 16-30 Hz component (same filter as the envelope)
plot(ax_F, t(w_in), sig(w_in),      'Color', [.70 .70 .70], 'LineWidth', .6, 'DisplayName', 'signal');
plot(ax_F, t(w_in), beta_sig(w_in), 'Color', cr,            'LineWidth', 1.0, 'DisplayName', 'beta band (16-30 Hz)');
xlim(ax_F,[t_c-0.4 t_c+0.4]);
ylF = [min(sig(w_in)) max(sig(w_in))]; ylim(ax_F, ylF + [-0.05 0.35]*diff(ylF));   % headroom for the legend
legend(ax_F,'Location','northeast','FontSize',7,'Box','off','TextColor','k');
xlabel(ax_F,'Time (s)'); ylabel(ax_F,'Amplitude');
title(ax_F,'Zoom on one burst (red box in A)','FontSize',8,'FontWeight','bold');
set_ax(ax_F); fix_title(ax_F); plabel(ax_F,'F');

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