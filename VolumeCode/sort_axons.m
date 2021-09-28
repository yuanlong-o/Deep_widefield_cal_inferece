function bg_proc = sort_axons(vol_params, axon_params, gp_bgvals, cell_pos)

% bg_proc = sort_axons(vol_params, axon_params, gp_bgvals, cell_pos)
% 
% This function ranndomly sorts length(gp_bgvals) axons into N_proc bins
%
% The inputs to this function are:
%   - vol_params - Struct containing parameters for the volume generation
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
%                   text outputs. 2 = detailed text outputs (default = 1)
%   - axon_params   - Struct containing parameters for background generation
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
%       .numbranches - Number of allowable branches for a single process
%                      (default = 20) 
%       .varbranches - Standard deviation of the number of branches per
%                      process (default = 5) 
%       .maxfill     - Voxel maximum occupation: fraction of volume that is
%                      to be filled by background processes (default = 0.7)
%       .N_proc      - Number of background components (default = 10)
%       .l           - Gaussian process length scale for correllation
%                      background processes over the image (default = 25) 
%       .rho         - Gaussian process variance parameter for correllation
%                      background processes over the image (default = 0.1) 
%   - gp_bgvals     - Cell array with number of rows equal to the number of
%                     background processes, which contains the locations and
%                     values of the processes at those locations
%   - cell_pos         - Position of cells within the volume
%
% 
% The ouptut to this function is:
%     bg_proc - The volume structure provided by simulate_neural_vol,
%               modified to have the background components.
% 
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs

vol_params  = check_vol_params(vol_params);                                % Check volume parameters
axon_params   = check_axon_params(axon_params);                            % Check background parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate background

N_proc    = axon_params.N_proc;                                            % Extract the number of correlated background components to create
vol_sz    = vol_params.vol_sz*vol_params.vres;
if vol_params.verbose > 0                                                  % Optional verbose output
    fprintf('Sorting axons...')
end

bg_proc  = cell(N_proc,2);                                                 % Initialize the cell array to store the background processes
if(size(bg_proc,1)>vol_params.N_neur+vol_params.N_den)
    N_comps = vol_params.N_neur+vol_params.N_den;
    gp_bgpos = zeros(size(gp_bgvals,1),3);
    for kk = 1:size(gp_bgvals,1)
        clear TMP_pos
        if(~isempty(gp_bgvals{kk,1}))
          [TMP_pos(:,1),TMP_pos(:,2),TMP_pos(:,3)] = ind2sub(vol_sz,gp_bgvals{kk,1});
          gp_bgpos(kk,:) = mean(TMP_pos);
        end
    end
  
    cell_pos2 = cell_pos(1:N_comps,:);
    dist_mat = sqrt((bsxfun(@minus,cell_pos2(:,1),gp_bgpos(:,1)').^2)+...
          (bsxfun(@minus,cell_pos2(:,2),gp_bgpos(:,2)').^2)+...
          (bsxfun(@minus,cell_pos2(:,3),gp_bgpos(:,3)').^2));
    
    idxlist = zeros(N_comps,1);
    for ii = 1:N_comps
        [~,idx] = min(dist_mat(ii,:));
        dist_mat(:,idx) = inf;
        bg_proc{ii,1} = gp_bgvals{idx,1};
        bg_proc{ii,2} = gp_bgvals{idx,2};
        idxlist(ii) = idx;
    end

    for kk = 1:size(gp_bgvals, 1)
        if(~ismember(kk,idxlist))
            index = N_comps+ceil((N_proc-N_comps)*rand);
            bg_proc{index,1} = cat(1,bg_proc{index,1},gp_bgvals{kk,1});
            bg_proc{index,2} = cat(1,bg_proc{index,2},gp_bgvals{kk,2});
        end
    end
else
    for kk = 1:size(gp_bgvals, 1)
        index = ceil(N_proc*rand);
        bg_proc{index,1} = cat(1,bg_proc{index,1},gp_bgvals{kk,1});
        bg_proc{index,2} = cat(1,bg_proc{index,2},gp_bgvals{kk,2});
    end
end
if vol_params.verbose > 0    
    fprintf('done.\n')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%