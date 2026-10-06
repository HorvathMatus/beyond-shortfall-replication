function d = local_paths()
% LOCAL_PATHS  Optional, machine-specific folders of third-party MATLAB code,
% added to the path by bs_paths.m. Copy this file to local_paths.m and edit
% the folders if the code is not placed in matlab/thirdparty/. The order
% matters: folders listed later take precedence.
d = {'/path/to/mfe-toolbox', ...
     '/path/to/ucsd_garch', ...
     '/path/to/etoolbox'};
end
