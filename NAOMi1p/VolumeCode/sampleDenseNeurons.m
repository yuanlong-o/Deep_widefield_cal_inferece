function [neur_locs,Vcell,Vnuc,Tri,rotAng] = sampleDenseNeurons(neur_params,vol_params,neur_ves)

% [neur_locs,Vcell,Vnuc,Tri,rotAng] = sampleNeurons(neur_params,vol_params,neur_ves)
%
% Function to sample the shapes and locations for all somas in a volume.
% The inputs to this function are:
%   - neur_params - Struct containing parameters for neuron generation
%       .n_samps     - Number of sphere samples to use in creating the mesh
%                      for generating soma and nucleus shapes (default =
%                      200) 
%       .l_scale     - length-scale for the isotropic GP of the soma
%                      shapes. This controls the shape `bumpiness' (default
%                      = 105) 
%       .p_scale     - Overall variance of the isotropic GP of the soma 
%                      shape. (default = 95) 
%       .avg_rad     - Average radius of each neuron in um (default =
%                      6.6 um) 
%       .nuc_fluorsc - Potential fluorescence in the nucleus (default =
%                      0)
%       .min_thic    - Minimum cytoplasmic thickness (default = 0.25)
%       .eccen       - Maximum eccentricity of neuron (default = 0.25)
%       .exts        - Parameters dictating the max/min of the soma radii
%                      (Default = [0.75,1.7])
%       .nexts       - Parameters dictating the extent to shrink and smooth
%                      the nucleus (Default = )
%       .neur_type   - Option for neuron type (Default 'pyr')
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
%   - neur_ves    - An array deliniating off-limit areas in the volume that
%                   are already occupied by vasculature
%
% The ouputs of this function are:
%   - neur_locs - A Kx3 array of the 3D locations for all K neurons
%   - Vcell     - A length K cell array where the k^th element contains the
%                 locations for the verticies defining the soma shape for
%                 each neuron
%   - Vnuc      - A length K cell array where the k^th element contains the
%                 locations for the verticies defining the nucleus shape 
%                 for each neuron
%   - Tri       - The triangulation for the surface mesh grids connecting
%                 the vertics in both Vcell and Vnuc
%   - rotAng    - Nx3 Vector with the rotation angle of the cell (Rx,Ry,Rz)
%
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing

vol_params  = check_vol_params(vol_params);                                % Check volume parameters
neur_params = check_neur_params(neur_params);                              % Check neuron parameters

eta = 1.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample neural soma shapes and locations

[x,y,z] = ndgrid(-ceil(vol_params.min_dist/2):ceil(vol_params.min_dist/2));% Set up a 3D grid to sample over
se      = strel(sqrt(x.^2 + y.^2 + z.^2) <= vol_params.min_dist/2);        % Find all the points that are within the minimum distance to the above points
neur_ves_trunc = imdilate(neur_ves,se);                                    % Dilate the above set to get a group of points to avoid later (based on the blodd vessel locations)

[V,Tri]=SpiralSampleSphere(neur_params.n_samps,false);                     % This function uses a spiral method to sample uniformly on a sphere
neur_params.S_samp = V;                                                    % Store sampling of the sphere
neur_params.Tri    = Tri;                                                  % Store triangulation of the sampling on the sphere

Vcell  = [];                                                               % Initialize the vertex points for each neuron soma's surface mesh 
Vnuc   = [];                                                               % Initialize the vertex points for each neuron nuceus' surface mesh
rotAng = [];
vol_sz = vol_params.vol_sz;                                                % Extract size of volume from parameter struct
[mesh_x, mesh_y, mesh_z] = meshgrid(...
              single(linspace(0,vol_sz(1),vol_sz(1)*vol_params.vres)),... 
              single(linspace(0,vol_sz(2),vol_sz(2)*vol_params.vres)),...
              single(linspace(0,vol_sz(3),vol_sz(3)*vol_params.vres)));    % Set up grid of points in the volume
mesh_x = permute(mesh_x,[2 1 3]);                                          % Rotate (x) mesh to be in line with how the volume dimensions are set up
mesh_y = permute(mesh_y,[2 1 3]);                                          % Rotate (y) mesh to be in line with how the volume dimensions are set up
mesh_z = permute(mesh_z,[2 1 3]);                                          % Rotate (z) mesh to be in line with how the volume dimensions are set up

vol_depth = vol_params.vol_depth*vol_params.vres;                          % Total volume depth
idx_good  = neur_ves_trunc(:,:,...
                          1+vol_depth:vol_depth+vol_sz(3)*vol_params.vres);% Isolate points in the vasculature array that are also in the volume
idx_good  = ~idx_good;                                                     % Remove points too close to the blood vessels from the list of allowed neural locations
idx_bad   = idx_good;                                                      % Remove points too close to the blood vessels from the list of allowed neural locations

if vol_params.verbose == 1    
    fprintf('Sampling random locations for the neurons...')
elseif vol_params.verbose > 1
    fprintf('Sampling random locations for the neurons...\n')
end
neur_locs = [Inf,Inf,Inf];                                                 % Initialize the list of neuron locations
kk        = 0;                                                             % Initialize the number of neurons currently in the volume to 0
while (sum(idx_good(:))>1)&&(size(Vcell,3)-isempty(Vcell)<vol_params.N_neur)% WHILE the requested number of neurons is not reached AND there is still room to place neurons....
    if vol_params.verbose >1                                               % Optional outputs
        tic
    end
    kk = kk+1;                                                             % Increment th enumber of neurons created by one

    % Draw random shape
    [V_tmp,Vnuc_tmp,~,rotAng_tmp] = generateNeuralBody(neur_params);           % Get the K^th neuron soma shape
    Vcell = cat(3,Vcell,V_tmp);                                            % Normalize so that each cell is ~10 pixels
    Vnuc  = cat(3,Vnuc,Vnuc_tmp);                                          % Add the nucleus to the list
    rotAng = cat(1,rotAng,rotAng_tmp);
    
    % Sample new location
    if sum(idx_good(:)) ~= 0
        idx_now = find(cumsum(idx_good(:)) == ...
                                randsample(sum(idx_good(:)),1),1,'first'); % Pick a new random point ...
    else
        idx_now = find(cumsum(idx_bad(:)) == ...
                                randsample(sum(idx_bad(:)),1),1,'first'); % Pick a new random point ...
    end
    new_pt    = [mesh_x(idx_now), mesh_y(idx_now), mesh_z(idx_now)];       % ... and get it's meshgrid location
    
    
    if(vol_params.N_neur==1) %center single neuron volumes
      new_pt = round(vol_params.vol_sz/2);
    end

    
    tmp_dist  = min(sqrt(sum(bsxfun(@minus,new_pt,neur_locs).^2,2)));      % Check the minimum distance to other points
    neur_locs = cat(1,neur_locs,new_pt);                                   % Add the new point to the list
    if vol_params.verbose > 1                                              % Optional outputs
        fprintf('%d, tmp_dist is %f...',size(Vcell,3), tmp_dist);              
    end
    idx_good(sqrt((mesh_x-new_pt(1)).^2 + (mesh_y-new_pt(2)).^2 + ...
        (mesh_z-new_pt(3)).^2) <= eta*vol_params.min_dist) = false;        % Exclude seed points too close to other neurons
    idx_bad(sqrt((mesh_x-new_pt(1)).^2 + (mesh_y-new_pt(2)).^2 + ...
        (mesh_z-new_pt(3)).^2) <= vol_params.min_dist) = false;            % Exclude seed points too close to other neurons
    idx_good = ~(idx_good|~idx_bad);
    if vol_params.verbose == 1                                             % Optional outputs
        fprintf('.');
    elseif vol_params.verbose >1
        Tdone = toc;
        fprintf('done (%f seconds).\n',Tdone);
    end
end 

neur_locs = single(neur_locs(2:end,:));                                    % Remove the initialization row and cast as a single
N_neur    = size(Vcell,3);                                                 % Store the number of neurons that actually fit into the volume
Vcell     = single(bsxfun(@plus,Vcell,reshape(neur_locs.',[1,3,N_neur]))); % Shift neurons to their new locations and cast as a single
Vnuc      = single(bsxfun(@plus,Vnuc,reshape(neur_locs.',[1,3,N_neur])));  % Shift nuclei to their new locations and cast as a single

if vol_params.verbose >= 1                                                 % Optional outputs
    fprintf('done.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
