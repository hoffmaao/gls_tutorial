function root = tutorial_setup()
%TUTORIAL_SETUP Put the tutorial's own code on the MATLAB path.
%
% Everything the tutorial needs ships in this folder, so it runs on a
% personal computer with only MATLAB installed:
%   vendor/+ptt   the fabric_anisotropy functions the lessons read and
%                 compare against (see vendor/README.md for the version)
%   +apres        ApRES reading, range processing and quad-pol helpers
%   project/      the build-your-own-estimator project
%
% Call it at the top of any lesson or script; calling it twice is harmless.

root = fileparts(mfilename('fullpath'));
dirs = {root, fullfile(root, 'vendor'), fullfile(root, 'project'), ...
        fullfile(root, 'project', 'student')};
% Add only what is missing, and never reorder: an instructor grading a
% student's own folder puts it first on the path, and it must stay first.
onpath = strsplit(path, pathsep);
for k = 1:numel(dirs)
  if exist(dirs{k}, 'dir') && ~any(strcmp(onpath, dirs{k}))
    addpath(dirs{k}, '-end');
  end
end
end
