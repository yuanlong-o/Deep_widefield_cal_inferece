function [gp_vals, neur_vol] = setCellFluoresence(vol_params, neur_params, dend_params, neur_num, neur_soma, neur_num_AD, neur_locs, neur_vol)

% [gp_vals, neur_vol] = setCellFluoresence(vol_params, neur_params, ...
%                    neur_num, neur_soma, neur_num_AD, neur_locs, neur_vol)
%
% This function populates the interior points of the various cells with
% non-uniform local fluorescences. This basically models the non-uniform
% distributoin of fluorescing proteins in the cell. The Inputs to this
% function are:
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
%   - dend_params - Struct containing parameters for dendrite simulation
%       .dtParams    - dendritic tree number,radius in um of branches
%                      (uniform across circle),radius in um of branches
%                      (uniform in z) (default = [40 150 50 1 10])
%       .atParams    - Apical dendrite number,radius in um of branches
%                      (uniform across circle),radius in um of branches
%                      (uniform in z),offset from center in um (default
%                      = [1 5 2 2 4]) 
%       .atParams2   - Through-volume apical dendrite number,radius in
%                      um of branches (uniform across circle),radius in
%                      um of branches (uniform in z),offset from center
%                      in um (default = = [1 5 2 2 6])
%       .dweight     - Weight for path planning randomness in the
%                      dendrites (default = 10) 
%       .bweight     - Weight for obstruction (default = 50)
%       .thicknessScale  - Scaling for dendrite thickness (int). Should be
%                      1 for 1um sampling,(4 for 0.5um sampling)
%                      (default = 0.75) 
%       .dims        - dims set at 10 um/space (default = [30 30 30])
%       .dimsSS      - dims subsampling factor (10 samples per dim
%                      grid) (default = [10 10 10]) 
%       .rallexp     - Rall exponent that controls the cahnge in size
%                      of dendrites at branching locations (default =
%                      1.5) 
%   - neur_num    - An array that has the value 'k' at each location in the
%                   volume occupied by the k^th neuron
%   - neur_soma   - An array that has the value 'k' at each location in the
%                   volume occupied by the k^th neural soma
%   - neur_num_AD - An array that has the value 'k' at each location in the
%                   volume occupied by the k^th apical dendrite
%   - neur_locs   - Kx3 array where each row is the location (in um) of one
%                   of the neurons in the volume
%   - neur_vol    - Array with each element containing the absolute base
%                   fluorescence (either soma, dendrite, vasculature, or
%                   background/neuropil) 
% 
% The ouputs of this function are:
%   - gp_vals  - Kx2 cell arraay where for the k^th "row", the first
%                element has the list of all the interior points inside the
%                k^th cell, and the second element has the local
%                fluorescence values for that cell
%   - neur_vol - Updated array with each element containing the absolute
%                base fluorescence (either soma, dendrite, vasculature, or
%                background/neuropil) 
%
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

vol_params  = check_vol_params(vol_params);                                % Check volume parameters
neur_params = check_neur_params(neur_params);                              % Check neuron parameters
dend_params = check_dend_params(dend_params);                              % Check dendrite parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set the fluorescence at each cell

if vol_params.verbose == 1
    fprintf('Setting Fluorescence Distribution.')
elseif vol_params.verbose >1
    fprintf('Setting Fluorescence Distribution...\n')
end
if vol_params.verbose >1
    fprintf('(part 1)...\n')
end

vol_sz  = vol_params.vol_sz;                                               % Extract the volume size 
N_neur  = vol_params.N_neur;                                               % Get the number of neurons in the volume
vres    = vol_params.vres;                                                 % Get the volume resolution (in samples/um)
gp_vals = cell(N_neur+vol_params.N_den,3);                                 % Initialize the array that will contain the locations of the cell's interior points and the local relative fluorescences
wtSc    = dend_params.weightScale;
flSc    = neur_params.fluor_dist;
%%
if vol_params.verbose >1
        tic
end
numcomps = N_neur+vol_params.N_den;
numvox = zeros(numcomps,1);
currvox = zeros(numcomps,1);
for kk = 1:numel(neur_num)
  if(neur_num(kk)>=1 && neur_num(kk)<= numcomps)
    numvox(neur_num(kk)) = numvox(neur_num(kk))+1;
  end
end
for kk = 1:numcomps
  gp_vals{kk,1} = zeros(numvox(kk),1,'int32');
end
for kk = 1:numel(neur_num)
  if(neur_num(kk)>=1 && neur_num(kk)<= numcomps)
    currvox(neur_num(kk)) = currvox(neur_num(kk))+1;
    gp_vals{neur_num(kk),1}(currvox(neur_num(kk))) = int32(kk);
  end
end
if vol_params.verbose >1
  toc
end

%%
for kk = 1:N_neur
    if vol_params.verbose >1
        tic
    end
    TMP      = masked_3DGP_v2(round(neur_params.avg_rad*6*vres*ones(1,3)),...
                                                flSc(1)*vres, flSc(2), 0); % Set up a GP for the fluorophore distribution inside each neuron
    TMP_loc  = gp_vals{kk,1};                                  % Get all the points in the cell. Save memory by only storing the locations...
    TMP_soma = TMP_loc(neur_soma(TMP_loc)==kk);
    TMP_AD = TMP_loc(neur_num_AD(TMP_loc)==kk);
%     TMP_soma = (vec(find(neur_soma==kk)));                                 % Get all the points in the soma. Save memory by only storing the locations...
%     TMP_AD   = (vec(find(neur_num_AD==kk)));                               % Get all the points in the dendrites. Save memory by only storing the locations...
    [~,TMP_idxs,~]    = intersect(TMP_loc,TMP_soma);                       % Find where the soma and interior points intersect
    [~,TMP_idxs_AD,~] = intersect(TMP_loc,TMP_AD);                         % Find where the dendrites and interior points intersect
    [lx,ly,lz]        = ind2sub(vol_sz*vres,TMP_soma);                     % Get index triples for all soma locations
    
    TMP_dist = bsxfun(@minus,cat(2,lx,ly,lz),floor(vres*neur_locs(kk,:)))...
                                                    + floor(size(TMP,1)/2);%
    TMP_dist = bsxfun(@max,TMP_dist,[1 1 1]);
    TMP_dist = bsxfun(@min,TMP_dist,size(TMP));
    TMP_vals = TMP(sub2ind(size(TMP),TMP_dist(:,1),TMP_dist(:,2),...
                                                           TMP_dist(:,3)));% 
    TMP_vals = 0.5*(TMP_vals - mean(TMP_vals))/max(abs(TMP_vals ...
                                                    - mean(TMP_vals))) + 1;% normalize and shift
    TMP_vals(isnan(TMP_vals)) = 1;
                                                  
    [rx,ry,rz]    = ind2sub(vol_sz*vres,TMP_loc);                          % Translate interior point locations of each neuron from single index to index triples
    TMP_sep       = sqrt((rx-vres*neur_locs(kk,1)).^2+(ry-vres*...
                         neur_locs(kk,2)).^2+(rz-vres*neur_locs(kk,3)).^2);% 
%     gp_vals{kk,1} = int32(TMP_loc);                                        % Save memory by only storing the locations...
    gp_vals{kk,2} = (wtSc(2)*exp(-TMP_sep/(vres*wtSc(1)))+(1-wtSc(2))).* ...
        (1-wtSc(3));                                                          % ... and the values
%         (1-wtSc(3)*rand(size(gp_vals{kk,1})));                             % ... and the values
    gp_vals{kk,2}(TMP_idxs)    = TMP_vals;                                 % Store the values at the somas
    gp_vals{kk,2}(TMP_idxs_AD) = ones(size(TMP_AD));                       % Store the values at the dendrites
    gp_vals{kk,3}              = false(size(gp_vals{kk,1}));               % Initialize an array that stores where the soma is
    gp_vals{kk,3}(TMP_idxs)    = true;                                     % Store where the soma only is
    if(~isempty(neur_vol))
      neur_vol(gp_vals{kk,1})    = gp_vals{kk,2};                            % Store the fluoresence values in the volume array
    end
    if vol_params.verbose == 1
        fprintf('.')
    elseif vol_params.verbose >1
        Tdone = toc;
        fprintf('%d (%f seconds).\n',kk,Tdone);
    end
end

if vol_params.verbose >1
    fprintf('(part 2)...\n')
end

for kk = (N_neur+1):(N_neur+vol_params.N_den)
    if vol_params.verbose >1
        tic
    end
%     gp_vals{kk,1}           = int32(vec(find(neur_num==kk)));              % Save memory by only storing the locations...
    gp_vals{kk,2}           = ones(size(gp_vals{kk,1}),'single');          % ... and the values
    gp_vals{kk,3}           = false(size(gp_vals{kk,1}));               
    if(~isempty(neur_vol))
      neur_vol(gp_vals{kk,1}) = gp_vals{kk,2}; 
    end
    if vol_params.verbose == 1
        fprintf('.')
    elseif vol_params.verbose >1
        Tdone = toc;
        fprintf('%d (%f seconds).\n',kk,Tdone);
    end
end

if vol_params.verbose >= 1
    fprintf('done.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%