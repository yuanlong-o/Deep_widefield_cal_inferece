function vol_params = check_vol_params(vol_params)

% vol_params = check_vol_params(vol_params)
%  
% This function checks the elements of the struct vol_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - vol_params      - Struct with parameters for the volume generation
%       .vol_sz       - 3-element vector with the size (in um) of the 
%                       volume to generate (default = 100x100x30um)
%       .min_dist     - Minimum distance between neurons (default = 15um)
%       .N_neur       - Number of neurons to generate (default = 50)
%       .vres         - resolution to simulate volume at (default = 2
%                       samples/um)
%       .N_den        - Width of apical dendrites (default = 10)
%       .N_bg         - Number of background/neuropil components to 
%                       simulate (default = 50)
%       .vol_depth    - Depth of the volume under the brain surface
%                       (default = 100um)
%       .dendrite_tau - Dendrite decay strength exponential distance
%                       (default = 5)
%       .verbose      - Level of verbosity in the output during the volume
%                       generation. Can be 0,1,2. 0 = no text updates, 1 =
%                       some text outputs. 2 = detailed text outputs
%                       (default = 1)
% 
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(vol_params)                                                     % Make sure that vol_params is a struct
    clear vol_params
    vol_params = struct;
end

dParams.vol_sz       = [100,100,50];                                      % Default Volume size is 100x100x50um
dParams.min_dist     = 16;                                                % Default minimum distance between neuron centers (Default = 16)
dParams.vres         = 2;                                                 % Default volume resolution(Default 10 points per unit)
dParams.N_bg         = 1e6;                                               % Default number of background processes
dParams.vol_depth    = 200;                                               % Default depth of the volume within the tissue is 100um
dParams.dendrite_tau = 5;                                                 % Default dendrite decay strength exponential distance is 5
dParams.verbose      = 1;                                                 % Default verbose level is 1


vol_params = setParams(dParams, vol_params);
if mod(vol_params.vol_sz(3),10)~=0
    vol_params.vol_sz(3) = 10*ceil(vol_params.vol_sz(3)/10);               % Ensure that the volume depth is a multiple of 10 (for indexing purposes in the dendrite code)
end

if (~isfield(vol_params,'N_neur'))||isempty(vol_params.N_neur)             % Default number of neurons
    if(isfield(vol_params,'neur_density')&&(~isempty(vol_params.neur_density)))
        vol_params.N_neur = ceil(vol_params.neur_density*prod(vol_params.vol_sz)/(1e9));
    else
        vol_params.neur_density = 1e5;
        vol_params.N_neur = ceil(vol_params.neur_density*prod(vol_params.vol_sz)/(1e9));
    end
else
    if(isfield(vol_params,'neur_density')&&(~isempty(vol_params.neur_density)))
%         warning('Both number and density of neurons provided. Simulation will ignore the number of neurons and stay with the neural density...')
%         vol_params.N_neur = ceil(vol_params.neur_density*prod(vol_params.vol_sz)/(1e9));
    else
        vol_params.neur_density = 1e9*vol_params.N_neur...
                                                /prod(vol_params.vol_sz);  % Calculate the neural density from the number of neurons
    end
end

if (~isfield(vol_params,'N_den'))||isempty(vol_params.N_den)               % Default number of dendrites
  if(isfield(vol_params,'AD_density')&&(~isempty(vol_params.AD_density)))
    vol_params.N_den = vol_params.AD_density*prod(vol_params.vol_sz(1:2))/(1e6);
  else
    vol_params.AD_density = 2e3;
    vol_params.N_den = vol_params.AD_density*prod(vol_params.vol_sz(1:2))/(1e6);
  end
end

if (~isfield(vol_params,'N_den'))||isempty(vol_params.N_den)               % Default number of apical dendrites
    vol_params.N_den = 10;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
