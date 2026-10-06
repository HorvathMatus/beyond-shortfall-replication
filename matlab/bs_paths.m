function base = bs_paths()
% BS_PATHS  Folder of the Beyond Shortfall MATLAB code, with a trailing file
% separator, and set-up of the MATLAB path. Called at the start of the
% numbered scripts; the scripts read their inputs from, and write their
% outputs to, this folder.
%
% Third-party code (see thirdparty/README.md) is added from
%   1. the subfolder thirdparty/ of this folder, if it exists, and
%   2. the folders returned by local_paths.m, if that file exists (an
%      optional, machine-specific list of toolbox folders; see
%      local_paths_example.m).
base = [fileparts(mfilename('fullpath')), filesep];
dirs = {fullfile(base, 'thirdparty')};
if exist(fullfile(base, 'local_paths.m'), 'file')
    dirs = [dirs, local_paths()];
end
for k = 1:numel(dirs)
    if exist(dirs{k}, 'dir'), addpath(genpath(dirs{k})); end
end
rehash path
end
