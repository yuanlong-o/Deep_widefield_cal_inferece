function [neur_soma, neur_vol, gp_nuc, gp_soma] = generateNeuralVolume(neur_params,vol_params,neur_locs,Vcell,Vnuc,neur_ves)

% function [neur_soma, neur_vol, gp_nuc] = generateNeuralVolume(neur_params,vol_params,neur_locs,Vcell,Vnuc)
% 
% Function to place the neural soma in a volume. The inputs to this
% function are:
%   - neur_params - Struct containing parameters for neuron generation
%       .n_samps     - Number of sphere samples to use in creating the mesh
%                      for generating soma and nucleus shapes (default =
%                      1000) 
%       .l_scale     - length-scale for the isotropic GP of the soma
%                      shapes. This controls the shape `bumpiness' (default
%                      = 105) 
%       .p_scale     - Overall variance of the isotropic GP of the soma 
%                      shape. (default = 95) 
%       .avg_rad     - Average radius of each neuron in um (default =
%                      5.5um) 
%       .nuc_fluorsc - Potential fluorescence in the nucleus (default =
%                      0.3)
%       .min_thic    - Minimum cytoplasmic thickness (default = 0.25)
%       .eccen       - Maximum eccentricity of neuron (default = 0.25)
%       .exts        - Parameters dictating the max/min of the soma radii
%                      (Default = [0.75,1.7])
%       .nexts       - Parameters dictating the extent to shrink and smooth
%                      the nucleus (Default = )
%       .neur_type   - Option for neuron type (Default 'pyr')
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
%   - neur_locs - Nx3 array containing the locations in the volume for each
%                 cell 
%   - Vcell     - Surface points of each cell's somas
%   - Vnuc      - Surface points of each cell's nuclei
%   - neur_ves  - Binary array indicating what pixels are occupied by
%                 blood vessels
% 
% The outputs for this function are:
%   - neur_soma - vol_params.vol_sz-sized array where the integers at each
%                 location indicate which neuron exists in that voxel
%   - neur_vol  - vol_params.vol_sz-sized array containing the base,
%                 relative, fluorescence levels at each voxel
%   - gp_nuc    - number-of-neurons-by-2 Cell array containing, for each
%                 neuron, the locations and fluorescence values for that
%                 neuron. THis array will be augmented with dendrite values
%                 as a third column later on
%   - gp_soma    - number-of-neurons-by-1 Cell array containing, for each
%                 neuron, the locations for the cell body of that neuron.
%
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

vol_params  = check_vol_params(vol_params);                                % Check volume parameters
neur_params = check_neur_params(neur_params);                              % Check neuron parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample soma locations in the volume

if vol_params.verbose >= 1
    fprintf('Setting up volume...')
end

vol_sz = vol_params.vol_sz;
vres = vol_params.vres;
% [mesh_x, mesh_y, mesh_z] = meshgrid(...
%                 single(linspace(0+1/vres/2,vol_sz(1)-1/vres/2,vol_sz(1)*vres)),... 
%                 single(linspace(0+1/vres/2,vol_sz(2)-1/vres/2,vol_sz(2)*vres)),...
%                 single(linspace(0+1/vres/2,vol_sz(3)-1/vres/2,vol_sz(3)*vres)));  % Set up grid of points in the volume
% mesh_x = permute(mesh_x,[2 1 3]);
% mesh_y = permute(mesh_y,[2 1 3]);
% mesh_z = permute(mesh_z,[2 1 3]);

neur_soma = zeros(vol_sz*vres, 'uint16');                                 % Initialize the array that stores where all the different somas are
neur_vol  = zeros(vol_sz*vres, 'single');                                 % Initialize the full neural voume array
gp_nuc    = cell(vol_params.N_neur,2);                                     % Set up a cell array to hold the GPs for the nuclei
gp_soma   = cell(vol_params.N_neur,1);                                     % Set up a cell array to hold the GPs for the nuclei
taken_pts = neur_ves;                                                      % Initialize an array with the locations not taken up by cells or vasculature
vol_depth = vol_params.vol_depth*vres;                          % Total volume depth
taken_pts = taken_pts(:,:,vol_depth+1:vol_depth+vol_sz(3)*vres);% Only keep the portion of the vasculature volume that intersects with the fully simulated volume
[~,Tri]   = SpiralSampleSphere(neur_params.n_samps,false);                 % Get a sampling of points on a sphere

if vol_params.verbose == 1
    fprintf('done.\nFinding interior points...')
elseif vol_params.verbose > 1
    fprintf('done.\nFinding interior points...\n')
end

for kk = 1:vol_params.N_neur
    tic
    max_ext     = ceil(max(sqrt(sum(bsxfun(@plus,Vcell(:,:,kk),...
                                                -neur_locs(kk,:)).^2,2)))); % Find the furthest point of the soma from the "origin"
                                              
    idx_pos = [round(vres*neur_locs(kk,1)), round(vres*neur_locs(kk,2)), round(vres*neur_locs(kk,3))];
    idxX = max(1,idx_pos(1)-vres*max_ext):min(idx_pos(1)+vres*max_ext,vol_sz(1)*vres);
    idxY = max(1,idx_pos(2)-vres*max_ext):min(idx_pos(2)+vres*max_ext,vol_sz(2)*vres);
    idxZ = max(1,idx_pos(3)-vres*max_ext):min(idx_pos(3)+vres*max_ext,vol_sz(3)*vres);

    [mesh_x, mesh_y, mesh_z] = meshgrid(...
                single((idxX-idx_pos(1))/vres),... 
                single((idxY-idx_pos(2))/vres),... 
                single((idxZ-idx_pos(3))/vres));                                              
    mesh_x = permute(mesh_x,[2 1 3]);
    mesh_y = permute(mesh_y,[2 1 3]);
    mesh_z = permute(mesh_z,[2 1 3]);
                                    
    idx_to_test = sqrt((mesh_x).^2+(mesh_y).^2 + (mesh_z).^2)  <= max_ext; % Search a subset of points to improve speed
    
    idx_tri = [mesh_x(idx_to_test)+idx_pos(1)/vres+1/vres/2,...
               mesh_y(idx_to_test)+idx_pos(2)/vres+1/vres/2,...
               mesh_z(idx_to_test)+idx_pos(3)/vres+1/vres/2];               
    TMP1 = intriangulation(Vcell(:,:,kk),Tri,idx_tri);                % Find all the points inside the K^th neuron triangulation                           
    TMP  = intriangulation(Vnuc(:,:,kk),Tri,idx_tri);                % Find the points inside the K^th nucleus to remove (no fluorescence)
               
%     max_ext     = max(sqrt(sum(bsxfun(@plus,Vcell(:,:,kk),...
%                                                 -neur_locs(kk,:)).^2,2))); % Find the furthest point of the soma from the "origin"
%     idx_to_test = sqrt((mesh_x-neur_locs(kk,1)).^2 ...
%         + (mesh_y-neur_locs(kk,2)).^2 + (mesh_z-neur_locs(kk,3)).^2) ...
%                                                                <= max_ext; % Search a subset of points to improve speed
%     TMP1 = intriangulation(Vcell(:,:,kk),Tri,[mesh_x(idx_to_test),...
%                  mesh_y(idx_to_test),mesh_z(idx_to_test)]);                % Find all the points inside the K^th neuron triangulation                           
%     TMP  = intriangulation(Vnuc(:,:,kk),Tri,[mesh_x(idx_to_test),...
%                  mesh_y(idx_to_test),mesh_z(idx_to_test)]);                % Find the points inside the K^th nucleus to remove (no fluorescence)
    TMP1a = false(size(idx_to_test));
    TMP1a(idx_to_test) = TMP1;
    TMPa = false(size(idx_to_test));
    TMPa(idx_to_test) = TMP;
               
    neur_idx               = TMP1a&(~TMPa);                                  % Remove nucleus points
    neur_idx               = neur_idx&(~(taken_pts(idxX,idxY,idxZ)));      % Remove points from potentially overlapping neurons
    taken_pts(idxX,idxY,idxZ) = taken_pts(idxX,idxY,idxZ)|neur_idx;        % Take the union of the points that are taken with the current index points
    [iX,iY,iZ] = ind2sub(size(neur_idx),find(neur_idx));
    iX = iX+idxX(1)-1;
    iY = iY+idxY(1)-1;
    iZ = iZ+idxZ(1)-1;
    neur_idx2              = sub2ind(vol_sz*vres,iX,iY,iZ);
    neur_soma(neur_idx2)    = uint16(kk);                                   % Number where the kth neuron is saved to the points inside the soma
    gp_soma{kk,1}          = int32(neur_idx2);                   % Find points in the nucleus

    neur_idx               = TMPa;                                          % Look at nucleus points
    [iX,iY,iZ] = ind2sub(size(neur_idx),find(neur_idx));
    iX = iX+idxX(1)-1;
    iY = iY+idxY(1)-1;
    iZ = iZ+idxZ(1)-1;
    neur_idx2              = sub2ind(vol_sz*vres,iX,iY,iZ);
    gp_nuc{kk,1}           = int32(neur_idx2);                   % Find points in the nucleus
    gp_nuc{kk,2}           = neur_params.nuc_fluorsc;                      % Set some fluorescence inside the nucleus
    neur_vol(gp_nuc{kk,1}) = gp_nuc{kk,2};                                 % Save the fluoresence values only (for memory considerations)

%     neur_idx               = false(vol_sz*vres);                          % Set up an array indexing the points inside the neuron
%     neur_idx(idx_to_test)  = TMP1&(~TMP);                                  % Remove nucleus points
%     neur_idx               = neur_idx&(~(taken_pts));                      % Remove points from potentially overlapping neurons
%     taken_pts              = taken_pts|neur_idx;                           % Take the union of the points that are taken with the current index points
%     neur_soma(neur_idx)    = uint16(kk);                                   % Number where the kth neuron is saved to the points inside the soma
%     gp_soma{kk,1}          = int32(vec(find(neur_idx)));                   % Find points in the nucleus
    
%     neur_idx               = false(vol_sz*vres);                          % Set up an array indexing the points inside the neuron
%     neur_idx(idx_to_test)  = TMP;                                          % Look at nucleus points
%     gp_nuc{kk,1}           = int32(vec(find(neur_idx)));                   % Find points in the nucleus
%     gp_nuc{kk,2}           = neur_params.nuc_fluorsc;                      % Set some fluorescence inside the nucleus
%     neur_vol(gp_nuc{kk,1}) = gp_nuc{kk,2};                                 % Save the fluoresence values only (for memory considerations)
    
    if vol_params.verbose == 1
        fprintf('.');
    elseif vol_params.verbose >1
        Tdone = toc;
        fprintf('%d done (%f seconds).\n',kk,Tdone);
    end
end

if vol_params.verbose >= 1
    fprintf('done.\n')
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
