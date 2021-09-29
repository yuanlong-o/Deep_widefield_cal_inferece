function [neur_vol,gp_bgvals, axon_params, vol_params] = ...
                       generate_axons(vol_params, axon_params, neur_vol, ...
                                                neur_num, gp_vals, gp_nuc, neur_vol_flag)

% [bg_pix,neur_vol,gp_bgvals, axon_params, vol_params] = ...
%                      generate_axons(vol_params, axon_params, neur_vol, ...
%                                               neur_num, gp_vals, gp_nuc)
% This function simulates background components. The inputs to this
% function are:
%   - vol_params  - Struct with parameters for the volume generation
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
%   - neur_vol    - An 3D array with the overall base fluorescence at each
%                   3D location (includes cells/vasculature)
%   - neur_num    - An array where the k^th neuron's locations in the
%                   volume (both the soma and dendrites) are deliniated by
%                   the value 'k'
%   - gp_vals     - Cell array that contains the locations and fluorescence
%                   values 
%   - gp_nuc      - Locations of the nucleus voxels in the volume
%
% The ouptut to this function is:
%   - neur_vol   - Updated 3D array with the overall base fluorescence at
%                  each 3D location (includes cells/vasculature and now also
%                  background/neuropil)
%   - gp_bgvals  - Cell array with number of rows equal to the number of
%                  background processes, which contains the locations and
%                  values of the processes at those locations
%   - axon_params  - Potentially updated background parameter struct
%   - vol_params - Potentially updated volume parameter struct
%
% 2016 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

vol_params  = check_vol_params(vol_params);                                % Check volume parameters
axon_params   = check_axon_params(axon_params);                                  % Check background parameters
if(nargin<7)
  neur_vol_flag = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate background

if vol_params.verbose == 1
  fprintf('Generating background fluorescence.')
elseif vol_params.verbose >1
  fprintf('Generating background fluorescence...\n')
end

bg_pix = (neur_num == 0);                                                  % Get the locations where the background can be
for kk = 1:size(gp_nuc,1)
  bg_pix(gp_nuc{kk,1}) = 0;                                                % Remove any zeros that may be in the nuclii
end
fillnum   = round((axon_params.maxfill)*(axon_params.maxvoxel)*sum(bg_pix(:)));   % Set the fill-number for filling the background with processes
volsize   = vol_params.vol_sz*vol_params.vres;                             % Extract the size of the neural volume
N_bg      = vol_params.N_bg;                                               % Get the number of neurons in the volume
gp_bgvals = cell(N_bg,2);                                                  % Initialize an array of background components

if vol_params.verbose >1                                                   % Optional verbose output
  fprintf('Initializing volume')
end
if(neur_vol_flag)
  neur_vol  = zeros(size(neur_vol),'single');                                % Initialize baseline neural volume
  for kk = 1:size(gp_vals,1)
    neur_vol(gp_vals{kk,1}) = gp_vals{kk,2};                                 % Initialize the new full neural volume to the set neural fluorescences. Do this for somas...
    if kk <= size(gp_nuc,1)
      neur_vol(gp_nuc{kk,1})  = gp_nuc{kk,2};                              % ... and nuclei (if needed)
    end
    if vol_params.verbose >=1                                                % Optional verbose output
      fprintf('.')
    end
  end
end
if vol_params.verbose >1                                                   % Optional verbose output
  fprintf('\n')
end

padsize = axon_params.padsize;
volpad = volsize+2*padsize;

M            = rand(volpad,'single');                           % 
M(padarray(bg_pix==0,padsize*[1 1 1],false,'both')) = realmax('single');   % 
clear bg_pix

if vol_params.verbose >1                                                   % Optional verbose output
  tic
end

j      = 1;                                                                % Initialize background process count
numit2 = 0;
nummax = 10000;
while((fillnum>0)&&(j<=N_bg)&&(numit2<nummax))
    bgpts  = [];                                                           % Initialize the list of points in the background as an empty vector
    numit2 = 0;                                                            % Set up a counter to test for stuck processes with nowhere to grow
    while((length(bgpts)<axon_params.minlength)&&(numit2<nummax))
        numit2 = numit2+1;                                                 % Increment counters of number of trials where there is nowhere to grow
        root   = ceil((volpad-2).*rand(1,3)+1);
        while(M(root(1),root(2),root(3))>(axon_params.fillweight*...
                                                         axon_params.maxvoxel))
            root = ceil((volpad-2).*rand(1,3)+1);
        end
        ends   = ceil(root + 2*axon_params.maxdist*vol_params.vres*...
                                                         (rand(1,3)-0.5)); % 
        if(ends(1)>volpad(1)); ends(1) = volpad(1); end
        if(ends(2)>volpad(2)); ends(2) = volpad(2); end
        if(ends(3)>volpad(3)); ends(3) = volpad(3); end
        if(ends(1)<1); ends(1) = 1; end
        if(ends(2)<1); ends(2) = 1; end
        if(ends(3)<1); ends(3) = 1; end
        bgpts = dendrite_randomwalk2(M,root,ends,axon_params.distsc,...
                         axon_params.maxlength,axon_params.fillweight,...
                                   axon_params.maxvoxel,axon_params.minlength);   % 
    end
    if ~isempty(bgpts)
      nbranches  = max(0,round(axon_params.numbranches + ...
                                           axon_params.varbranches*randn));
      for i = 1:nbranches
        bgpts2   = [];
        numit = 0;
        while(length(bgpts2)<axon_params.minlength && numit<100)
          numit = numit+1;
          root  = bgpts(ceil(rand*length(bgpts)),:);
          while(root(1) == 1 || root(1) == volpad(1) || root(2) == 1 ||...
                                root(2) == volpad(2) || root(3) == 1 ||...
                                                     root(3) == volpad(3))
            root = bgpts(ceil(rand*length(bgpts)),:);
          end
          ends   = ceil(root + int32(2*axon_params.maxdist*...
                                        vol_params.vres*(rand(1,3)-0.5)));
          if(ends(1)>volpad(1)); ends(1) = volpad(1); end
          if(ends(2)>volpad(2)); ends(2) = volpad(2); end
          if(ends(3)>volpad(3)); ends(3) = volpad(3); end
          if(ends(1)<1); ends(1) = 1; end
          if(ends(2)<1); ends(2) = 1; end
          if(ends(3)<1); ends(3) = 1; end
          bgpts2 = dendrite_randomwalk2(M,root,ends,axon_params.distsc,...
                            axon_params.maxlength,axon_params.fillweight,...
                                     axon_params.maxvoxel,axon_params.minlength);
        end
        bgpts = cat(1,bgpts,bgpts2);
      end
      bgpts = bgpts-padsize;
      TMPidxs = (bgpts(:,1)<=0)|(bgpts(:,1)>volsize(1)) ...
          |(bgpts(:,2)<=0)|(bgpts(:,2)>volsize(2)) ...
          |(bgpts(:,3)<=0)|(bgpts(:,3)>volsize(3));
      
      bgpts(TMPidxs,:) = [];      
      if(~isempty(bgpts))
        gp_bgvals{j,1} = bgpts(:,1)+(bgpts(:,2)-1)*volsize(1)+...
                                         (bgpts(:,3)-1)*volsize(1)*volsize(2);
        gp_bgvals{j,2} = single((1/axon_params.maxel)*ones(size(bgpts,1),1)* ... 
                          max(0,1+axon_params.varfill*randn));
        fillnum        = fillnum-size(bgpts,1);
        if(neur_vol_flag)
          neur_vol(gp_bgvals{j,1}) = neur_vol(gp_bgvals{j,1})+gp_bgvals{j,2};
        end
        j      = j+1;                                                          % Increment background process count
      end
      if vol_params.verbose >1                                             % Optional verbose output
        if(mod(j,1000)==0)
          Tdone = toc;
          fprintf('%d (%f seconds).\n',j,Tdone);
        end
      end
    end
end
if(j>N_bg)
  j = N_bg;
end
vol_params.N_bg = j;                                                       % Save the number of background pieces generated
gp_bgvals       = gp_bgvals(1:j,:);                                        % Only output the number of components actually generated (initialized to a much more optomistic number)

if(~neur_vol_flag)
  neur_vol = [];
end
if vol_params.verbose >= 1
  fprintf('done.\n')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%