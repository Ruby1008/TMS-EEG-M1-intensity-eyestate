function figS5_oldS4_1108n_v2(variant)
% v2 (user 2026-10-28): legend and peak labels set to the tick-label size (10.25 pt), output FigureS5_1108n_v2.png; otherwise as figS5_oldS4_1108n.m.
% _1108n copy (user 2026-10-28, Figures_2026-11-08n) of qc/figures_a3/figS4_a3.m (old Figure S4 = new Figure S5, IMF instantaneous-frequency distributions).
% Kept from figS4_a3.m: the source .fig (opened read-only, not saved back), Times New Roman, line widths (solid 1.5 pt, dashed 1.0 pt), print at 600 dpi.
% Changes (D): x label 'Frequency' -> 'Frequency (Hz, log2 scale)' (ticks unchanged), the title 'IMF Distribution' removed, axis font sizes as Figure S1
% (Figure S1: axis labels 9 pt, tick labels 8 pt on a 7.315-in wide figure; this figure is 9.375 in wide):
%   variant 'scaled' (default, main output): sizes x 9.375 / 7.315, so that both figures show the same text size when placed at the same width
%                                           (labels 11.53 pt, ticks 10.25 pt) -> Figures_2026-11-08n/FigureS5_1108n.png
%   variant 'pt'     (alternative):         the same nominal sizes as Figure S1 (labels 9 pt, ticks 8 pt) -> _check/alt_fontsize_pt/FigureS5_1108n_pt.png
% Values, peak labels (12 pt), dashed mean lines, legend (14 pt, manual size) unchanged.
data_dir = 'path/to/data';   % release: set this path
if nargin < 1, variant = 'scaled'; end
set(groot, 'defaultAxesFontName', 'Times New Roman', 'defaultTextFontName', 'Times New Roman', 'defaultLegendFontName', 'Times New Roman');
src = [data_dir '\Revision_2026\Manuscript_Figures\Fig4_new_main17.fig'];
N = [data_dir '\Revision_2026\Manuscript_Figures\Figures_2026-11-08n'];
fig = openfig(src, 'invisible');
set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
L = findall(fig, 'Type', 'line'); w0 = get(L, 'LineWidth'); if iscell(w0), w0 = cell2mat(w0); end
fprintf('line widths before: %s\n', mat2str(unique(w0)'));
for k = 1:numel(L)
    if strcmp(get(L(k), 'LineStyle'), '--') || strcmp(get(L(k), 'LineStyle'), ':'), set(L(k), 'LineWidth', 1.0); else, set(L(k), 'LineWidth', 1.5); end
end
% ---- 1108n changes
ax = findall(fig, 'Type', 'axes'); assert(numel(ax) == 1); lg = findall(fig, 'Type', 'legend'); lgs = lg.FontSize;
tx = findall(ax, 'Type', 'text'); tx = tx(~ismember(tx, [ax.XLabel; ax.YLabel; ax.Title])); txs = get(tx, 'FontSize');   % peak labels only
r = 9.375 / 7.315; if strcmp(variant, 'pt'), r = 1; end
ax.FontSize = 8 * r;                                   % tick labels
ax.XLabel.String = 'Frequency (Hz, log2 scale)'; ax.XLabel.FontSize = 9 * r; ax.YLabel.FontSize = 9 * r;
ax.Title.String = '';                                  % title removed
lg.FontSize = 8 * r; set(tx(isvalid(tx)), 'FontSize', 8 * r);   % v2: legend and peak labels = tick-label size
fprintf('variant %s: ticks %.2f pt, labels %.2f pt, legend %.0f pt\n', variant, ax.FontSize, ax.XLabel.FontSize, lg.FontSize);
if strcmp(variant, 'pt'), d = fullfile(N, '_check', 'alt_fontsize_pt'); if ~exist(d, 'dir'), mkdir(d); end; out = fullfile(d, 'FigureS5_1108n_v2_pt');
else, out = fullfile(N, 'FigureS5_1108n_v2'); end
print(fig, [out '.png'], '-dpng', '-r600'); close(fig);
disp(['saved ' out]);
end
