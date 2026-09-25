function run_all(lessons)
%RUN_ALL Run the tutorial lessons headlessly and save their figures.
%
% run_all            runs every lessonNN_*.m in order
% run_all([5 6])     runs only lessons 5 and 6
%
% Figures are written to figs/lessonNN_k.png. From a shell:
%   matlab -batch "run_all"
%
% The lessons themselves are meant to be stepped through section by
% section in the MATLAB editor; this runner is for checking that they all
% still run end to end.

here = fileparts(mfilename('fullpath'));
files = dir(fullfile(here, 'lesson*.m'));
files = sort({files.name});
if nargin >= 1
  keep = false(size(files));
  for k = 1:numel(files)
    keep(k) = ismember(sscanf(files{k}, 'lesson%d'), lessons);
  end
  files = files(keep);
end

outdir = fullfile(here, 'figs');
if ~exist(outdir, 'dir'), mkdir(outdir); end

for k = 1:numel(files)
  name = files{k}(1:end-2);
  fprintf('\n==================== %s ====================\n', name);
  t0 = tic;
  run_one(fullfile(here, files{k}));
  figs = findobj(groot, 'Type', 'figure');
  [~, ord] = sort([figs.Number]);
  figs = figs(ord);
  for j = 1:numel(figs)
    if isprop(figs(j), 'Theme'), figs(j).Theme = 'light'; end   % R2025a+ follows the OS theme
    for ax = findall(figs(j), 'Type', 'axes').'
      ax.Toolbar.Visible = 'off';
    end
    exportgraphics(figs(j), fullfile(outdir, sprintf('%s_%d.png', name(1:8), j)), ...
      'Resolution', 110);
  end
  close all;
  fprintf('---- %s done in %.1f s, %d figure(s)\n', name, toc(t0), numel(figs));
end
end

function run_one(path)
% Run in a separate workspace so each lesson's clear does not touch ours.
run(path);
end
