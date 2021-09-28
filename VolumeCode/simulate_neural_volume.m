function [vol_out,vol_params,neur_params,vasc_params,dend_params,bg_params,...
    axon_params] = simulate_neural_volume(vol_params, neur_params, ...
    vasc_params, dend_params, bg_params, axon_params, psf_params, varargin)

% [vol_out,vol_params,neur_params,vasc_params,dend_params,axon_params] ...
%  = simulate_neural_volume(vol_params, neur_params, vasc_params, ...
%                                 dend_params, bg_params, axon_params, debug_opt) 
%
% Function to create a volume. This is the main function to create a neural
% volume, inclusive of neural somas and dendrites, vasculature, and
% background/neuropil components. The volume created here can be used in
% conjunction with genearted time-traces to create a full volume of
% activity to be scanned in simulation. 
% 
% The inputs to this function are:
%   - vol_params  - Struct containing parameters for the volume generation
%       .vol_sz   - 3-element vector with the size (in um) of the volume to
%                   generate (default = 100x100x30um)
%       .min_dist - Minimum distance between neurons (default = 15um)
%       .N_neur   - Number of neurons to generate (default = 50)
%       .vres     - resolution to simulate volume at (default = 2
%                   samples/um)
%       .N_den    - Width of apical dendrites (default = 10)
%       .N_bg     - Number of background/neuropil components to simulate
%                   (default = 50)
%       .vol_depth- Depth of the volume under the brain surface (default =
%                   100um)
%       .verbose  - Level of verbosity in the output during the volume
%                   generation. Can be 0,1,2. 0 = no text updates, 1 = some
%                   some text outputs. 2 = detailed text outputs (default =
%                   1)
%   - neur_params - Struct containing parameters for neuron generation
%       .n_samps  - Number of sphere samples to use in creating the mesh
%                   for generating soma and nucleus shapes (default = 1000)
%       .l_scale  - length-scale for the isotropic GP of the soma shapes.
%                   This controls the shape `bumpiness' (default = 5)
%       .p_scale  - Overall variance of the isotropic GP of the soma shape.
%                   (default = 3.4) 
%       .avg_rad  - Average radius of each neuron in um (default = 5.5um)
%       .nuc_fluorsc - Potential fluorescence in the nucleus (default =
%                      0.3)
%       .min_thic - Minimum cytoplasmic thickness (default = 1)
%       .eccen    - Maximum eccentricity of neuron (default = 0.25)
%   - vasc_params - Struct containing parameters for vasculature simulation
%       .ves_shift       - 3-vector of the amount of wobble allowed for
%                          blood vessels (default = [30 15 15] um)
%       .depth_vasc      - Depth into tissue for which vasculature is
%                          simulated (default = 200um) 
%       .depth_surf      - Depth into tissue of surface vasculature
%                          (default = 15um) 
%       .distWeightScale - Scaling factor for weight of node distance
%                          scaling (default = 2) 
%       .randWeightScale - scaling factor for weight of nodes (additional
%                          variability) (default = 0.1) 
%       .cappAmpScale    - scaling factor for weights (capillaries,
%                          lateral) (default = 0.5) 
%       .cappAmpZscale   - scaling factor for weights (capillaries, axial)
%                          (default = 6) 
%       .vesSize         - vessel radius (surface, axial, capillaries) in
%                          um (default = [15 15 5] um) 
%       .vesFreq         - blood vessel freqency in microns (default = 
%                          [250 150 65] um) 
%       .vesNumScale     - blood vessel number random scaling factor
%                          (default = 0.2) 
%   - dend_params - Struct containing parameters for dendrite simulation
%       .dtParams        - dendritic tree number,radius in um of branches
%                          (uniform across circle),radius in um of branches
%                          (uniform in z) (default = [35 100 50 1])
%       .atParams        - Apical dendrite number,radius in um of branches
%                          (uniform across circle),radius in um of branches
%                          (uniform in z),offset from center in um (default
%                          = [1 5 2 2 4]) 
%       .dweight         - Weight for path planning randomness in the
%                          dendrites (default = 10) 
%       .bweight         - Weight for obstruction (default = 50)
%       .thicknessScale  - Scaling for dendrite thickness (int). Should be
%                          1 for 1um sampling,(4 for 0.5um sampling)
%                          (default = 0.75) 
%       .dims            - dims set at 10 um per space (default = [20 20 20])
%       .dimsSS          - dims subsampling factor (10 samples per dim
%                          grid) (default = [10 10 10]) 
%   - axon_params   - axon_params   - Struct containing parameters for background generation
%       .distvar         - (default = 5)
%       .distscale       - (default = 2)
%       .nanchors        - (default = 5)
%       .numptssc        - (default = 20)
%   - debug_opt   - [Optional] true/false flag to set debug mode (extra
%                   outputs and maximum verbosity - default = false)
%
% The outputs of this function are:
%   - vol_out     - Output struct containint the following quantitites:
%       .neur_vol - Fluorescence over the entire volume
%       .gp_vals  - Cellular fluorescence distribution
%       .gp_back  - Extra-cellular fluorescence distribution
%       .neur_all - Binary index of all neruron spatial locations
%       .neur_idx - Binary index of each neuron's locations
%       .Vcell    - Mesh grid points for the neurons
%       .Vnuc     - Mesh grid points for the neural nucleus
%       .Tri      - Triangulation of the mesh grid point .Vcell and .Vnuc.
% 
% 2016 - Adam Charles & Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Checking

if nargin > 7
    debug_opt = varargin{1};
else
    debug_opt = false;                                                     % Default to no debug mode
end
%%
vol_params  = check_vol_params(vol_params);                                % Check volume parameters
neur_params = check_neur_params(neur_params);                              % Check neuron parameters
vasc_params = check_vasc_params(vasc_params);                              % Check vasculature parameters
dend_params = check_dend_params(dend_params);                              % Check dendrite parameters
bg_params   = check_bg_params(bg_params);                                  % Check background parameters
axon_params = check_axon_params(axon_params);                              % Check axon parameters

if debug_opt
    vol_params.verbose = 2;                                                % If in debug mode, set the maximum verbosity level
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Blood Vessel Simulation

if(~isfield(vol_params,'vasc_sz'))||isempty(vol_params.vasc_sz)
  vol_params.vasc_sz = gaussianBeamSize(psf_params,vol_params.vol_depth+ ...
   vol_params.vol_sz(3)/2)+vol_params.vol_sz+[0 0 1]*vol_params.vol_depth; % Set up vasculature size (i.e. add in some space on top for the surface vasculature)
end
if(vasc_params.flag)
  [neur_ves,vasc_params,neur_ves_all] = simulatebloodvessels(vol_params,...
                                                             vasc_params); % Create a 3D volume with blood vessels (Large verticle vessels & small horizontal vessels)
else
  neur_ves     = false(vol_params.vres*(vol_params.vol_sz...
                                          +[0 0 1]*vol_params.vol_depth)); % If no vasculature is required, ...
  neur_ves_all = false(vol_params.vres*vol_params.vasc_sz);                %        ... then set up a fully empty volume
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample a set of Neurons

[neur_locs, Vcell, Vnuc, Tri, rotAng] = sampleDenseNeurons(neur_params, vol_params,...
                                                                neur_ves); % Sample shapes and locations for each neuron in the volume
vol_params.N_neur = size(Vcell,3);                                         % Extract the number of neurons, just in case not all of them fit into the volume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up volume via a 3-D grid of points

[neur_soma, neur_vol, gp_nuc, gp_soma] = generateNeuralVolume(neur_params,...
                            vol_params,neur_locs, Vcell, Vnuc, neur_ves);  % Place Somas
neur_ves(find(neur_soma(:)>0)+(vol_params.vol_depth*vol_params.vres)*...
                         prod(vol_params.vol_sz(1:2)*vol_params.vres)) = 0;% Set the location of somas to zero in the neural vessel location array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grow dendrites

[neur_num, cellVolumeAD, dend_params,gp_soma] = growNeuronDendrites(vol_params,...
                     dend_params, neur_soma, neur_ves, neur_locs, gp_nuc, gp_soma, rotAng); % Grow out dendrites

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add apical dendrites from bottom of volume to top

[neur_num, neur_num_AD, dend_params] = growApicalDendrites(vol_params, ...
                  dend_params, neur_num, cellVolumeAD, gp_nuc, gp_soma); % Grow out apical dendrites
clear cellVolumeAD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up non-uniform fluorescence for each neuron

[gp_vals, neur_vol] = setCellFluoresence(vol_params, neur_params, dend_params, ...
                   neur_num, neur_soma, neur_num_AD, neur_locs, neur_vol); % Set interior cell fluorescence distributions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up background non-uniform fluorescence

if(bg_params.flag)
[neur_num,neur_vol,vol_params,gp_vals,neur_locs] = ...
      generate_bgdendrites(vol_params, bg_params, dend_params, neur_vol, neur_num,...
                                                 gp_vals, gp_nuc,neur_locs); % Simulate the background/neuropil 
end

if(axon_params.flag)                                                     
  [neur_vol,gp_bgvals, axon_params] = ...
      generate_axons(vol_params, axon_params, neur_vol, neur_num, ...
                                                         gp_vals, gp_nuc); % Simulate the background/neuropil 
  axon_params.N_proc = size(gp_vals,1);
  bg_proc = sort_axons(vol_params, axon_params, gp_bgvals, neur_locs*vol_params.vres);                                                     
%   bg_proc = correlate_background(vol_params, axon_params, neur_ves, ...
%                                               gp_bgvals, gp_vals, gp_nuc); % Correlate the background components into a smaller number of processes
else
  gp_bgvals = [];
  bg_proc = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up output struct

if vol_params.verbose >= 1
    fprintf('Setting up output struct...')
end
vol_out.neur_vol     = neur_vol;                                            % Position and base fluoresence of neural volume 
clear neur_vol
vol_out.gp_nuc       = gp_nuc;                                              % Gaussian Process fluoresence distribution inside each nucleus
clear gp_nuc
vol_out.gp_soma      = gp_soma;                                             % Gaussian Process fluoresence distribution inside each soma
clear gp_soma
vol_out.gp_vals      = gp_vals;                                             % Gaussian Process fluoresence distribution inside each cell (soma+dendrites)
clear gp_vals
vol_out.neur_ves     = neur_ves;                                            % Locations of neural blood vessels
clear neur_ves
vol_out.bg_proc      = bg_proc;                                             % The full background processes (locations and values)
clear bg_proc
vol_out.neur_ves_all = neur_ves_all;                                        % Full vasculature volume simulated, for use in optical propagation
clear neur_ves_all
vol_out.locs         = neur_locs;                                           % Locations belonging to each neuron (in microns)
clear neur_locs

if debug_opt                                                               % --- These outputs are not needed, but can be output for debugging purposes
    vol_out.neur_num_AD = neur_num_AD;                                     %  | -
    clear neur_num_AD
    vol_out.neur_soma   = neur_soma;                                       %  | - 
    clear neur_soma
    vol_out.neur_num    = neur_num;                                        %  | - Numbered location of what neuron each voxel is a part of
    clear neur_num
    vol_out.gp_bgvals   = gp_bgvals;                                       %  | - The locations and values of all individual processes
    clear gp_bgvals
    vol_out.Vcell       = Vcell;                                           %  | - Cell soma vertices
    clear Vcell
    vol_out.Vnuc        = Vnuc;                                            %  | - Cell nucleus vertices
    clear Vnuc
    vol_out.Tri         = Tri;                                             %  | - Cell soma and nucleus triangulation
    clear Tri
end                                                                        % ---

if vol_params.verbose >= 1
    fprintf('done.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%