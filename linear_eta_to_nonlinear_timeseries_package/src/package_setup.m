function package_setup()
%PACKAGE_SETUP Add package-local and repo dependency paths.

root_dir = fileparts(fileparts(mfilename('fullpath')));
src_dir = fullfile(root_dir, 'src');
config_dir = fullfile(root_dir, 'config');
examples_dir = fullfile(root_dir, 'examples');
deps_dir = fullfile(root_dir, 'deps');
vwa_local_dir = fullfile(deps_dir, 'vwa');
mf12_local_dir = fullfile(deps_dir, 'mf12');
vwa_fallback_dir = fullfile(fileparts(root_dir), 'test functions for VWA Opensource');
mf12_fallback_dir = fullfile(fileparts(root_dir), 'irregularWavesMF12', 'Source');

addpath(src_dir);
addpath(config_dir);
addpath(examples_dir);

if isfolder(vwa_local_dir)
    addpath(vwa_local_dir);
elseif isfolder(vwa_fallback_dir)
    addpath(vwa_fallback_dir);
else
    error(['Missing VWA helper directory. Expected either package-local ' ...
        'dependency folder or fallback repo folder:\n  %s\n  %s'], ...
        vwa_local_dir, vwa_fallback_dir);
end

if isfolder(mf12_local_dir)
    addpath(mf12_local_dir);
elseif isfolder(mf12_fallback_dir)
    addpath(mf12_fallback_dir);
else
    error(['Missing MF12 source directory. Expected either package-local ' ...
        'dependency folder or fallback repo folder:\n  %s\n  %s'], ...
        mf12_local_dir, mf12_fallback_dir);
end
end
