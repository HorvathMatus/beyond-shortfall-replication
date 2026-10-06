function base = bs_paths()
% BS_PATHS  Folder of the Beyond Shortfall MATLAB code, with a trailing file
% separator, and set-up of the MATLAB path. Called at the start of the
% numbered scripts; the scripts read their inputs from, and write their
% outputs to, this folder.
%
% Third-party code (see thirdparty/README.md) is added in this order, later
% folders taking precedence:
%   1. the folders returned by local_paths.m, if that file exists (an
%      optional, machine-specific list of toolbox folders; see
%      local_paths_example.m),
%   2. the subfolders of thirdparty/ other than thirdparty/patton,
%   3. thirdparty/patton, the replication files of Patton, Ziegel and Chen
%      (2019), which must take precedence over the files of the same name in
%      the UCSD GARCH toolbox (garchpq.m, skewtdis_*.m) and the Econometrics
%      Toolbox (nwest.m, mprint.m).
base = [fileparts(mfilename('fullpath')), filesep];
dirs = {};
if exist(fullfile(base, 'local_paths.m'), 'file')
    dirs = local_paths();
end
tp = fullfile(base, 'thirdparty');
if exist(tp, 'dir')
    s = dir(tp);
    names = {s([s.isdir]).name};
    names = names(~startsWith(names, '.') & ~strcmp(names, 'patton'));
    dirs = [dirs, fullfile(tp, names)];
    if exist(fullfile(tp, 'patton'), 'dir'), dirs{end+1} = fullfile(tp, 'patton'); end
end
for k = 1:numel(dirs)
    if exist(dirs{k}, 'dir'), addpath(genpath(dirs{k})); end
end
rehash path
end
