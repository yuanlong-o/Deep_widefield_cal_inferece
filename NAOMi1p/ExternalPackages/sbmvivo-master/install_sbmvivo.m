function install_sbmvivo
%update the path
basep = fileparts(mfilename('fullpath'));
if basep(end) ~= filesep, basep = [basep filesep]; end
p = { ...
    '' ...  % base directory    
    ['mex' filesep 'sbm_cuda'] ...
    ['mex' filesep 'sbm_fit'] ...
    'c2swrapper' ...
    'MLspikewrapper' ...
    'util' ...
    'analysis' ...
    'gui' ...
    'optelecdataformat'};
addpath(basep);
for k = 1:numel(p)
    
    addpath([basep p{k}]);
    
end
savepath

%Try to add precompiled mex files to the path if mex files have not already
%been copied or compiled:
if isempty(mexext)
    
    warning('unrecognized system architecture: %s, cannot install precompiled mex files', computer);
    return;
    
end

precomp_dir = [basep 'mex' filesep 'sbm_cuda' filesep 'precompiled' filesep];
d = dir([precomp_dir '*.' mexext]);
d = {d.name};
for k = 1:numel(d)
    
    target = [basep 'mex' filesep 'sbm_cuda' filesep d{k}];
    if exist(target, 'file') == 0  % file does not exist
        
        copyfile([precomp_dir d{k}], target);
        
    end
    
end

%FIXME check whether mex files are compiled, and test them if GPU is
%present
