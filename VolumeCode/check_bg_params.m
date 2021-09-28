function bg_params = check_bg_params(bg_params)

% bg_params = check_bg_params(bg_params)
%  
% This function checks the elements of the struct bg_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - bg_params   - Struct containing parameters for background generation
%       .flag        - Flag for generation of background dendrites
%       .distsc      - Parameter to determine how directed the random walk 
%                      to generate background processes is Higher values
%                      are more directed/less random (default = 0.5)
%       .fillweight  - Maximum length for a single process branch (default
%                      = 100 um)
%       .maxlength   - Maximum length for background processes (default =
%                      200 um) 
%       .minlength   - Minimum length for background processes (default =
%                      10 um) 
%       .maxdist     - Maximum distance to the end of a process (default =
%                      100 um) 
%       .maxel       - Max number of axons per voxel (default = 8)
%
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(bg_params)                                                      % Make sure that bg_params is a struct
    clear bg_params
    bg_params = struct;
end

dParams.flag       = 1;
dParams.distsc     = 0.5;
dParams.fillweight = 100;
dParams.maxlength  = 200;
dParams.minlength  = 10;
dParams.maxdist    = 100;
dParams.maxel      = 1;

bg_params = setParams(dParams,bg_params);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
