function [neur_num,neur_vol,vol_params,gp_vals,neur_locs] = ...
  generate_bgdendrites(vol_params, bg_params, dend_params, neur_vol, neur_num, gp_vals, gp_nuc, neur_locs,neur_vol_flag)

% [neur_num,neur_vol,vol_params] = generate_bgdendrites(vol_params, ...
%              dend_params, neur_vol, neur_num, gp_vals, gp_nuc, neur_locs)
%
% This function simulates background components. The inputs to this
% function are:
%   - vol_params  - Struct containing parameters for the volume generation
%      .vol_sz    - 3-element vector with the size (in um) of the volume to
%                   generate (default = 100x100x30um)
%      .min_dist  - Minimum distance between neurons (default = 15um)
%      .N_neur    - Number of neurons to generate (default = 50)
%      .vres      - resolution to simulate volume at (default = 2
%                   samples/um)
%      .N_den     - Width of apical dendrites (default = 10)
%      .vol_depth - Depth of the volume under the brain surface (default =
%                   100um)
%      .verbose   - Level of verbosity in the output during the volume
%                   generation. Can be 0,1,2. 0 = no text updates, 1 = some
%                   some text outputs. 2 = detailed text outputs (default =
%                   1)
%   - bg_params   - Struct containing parameters for background generation
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
%   - gp_vals     - Cell array of index locations and cellular fluorescence
%                   distributions of all the neurons in the volume
%   - gp_nuc      - Cell array contianing the index locations of the nuclei
%   - neur_locs   - Neuron locations (in units of microns)
%
% The ouptut to this function is:
%   - neur_num   - vol_params.vol_sz-sized array containing the number of
%                  the neuron occupying each voxel
%   - neur_vol   - vol_params.vol_sz-sized array containing the base,
%                  relative, fluorescence levels at each voxel
%   - vol_params - Updated 3D array with the overall base fluorescence at
%                  each 3D location (includes cells/vasculature and now also
%                  background/neuropil)
%
% 2016 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

vol_params  = check_vol_params(vol_params);                                % Check volume parameters
dend_params = check_dend_params(dend_params);                            % Check dendrite parameters
bg_params   = check_bg_params(bg_params);                                  % Check background parameters
if(nargin<9)
  neur_vol_flag = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate background

if(nargin<8)
  neur_locs = [];
end
  
if vol_params.verbose == 1
  fprintf('Generating background fluorescence.')
elseif vol_params.verbose >1
  fprintf('Generating background fluorescence...\n')
end

bg_pix = (neur_num == 0);                                                  % Get the locations where the background can be
for i = 1:vol_params.N_neur
  bg_pix(gp_nuc{i,1}) = 0;                                                % Remove any zeros that may be in the nuclii
end

vres           = vol_params.vres;                                          % Extract the delta volume size
dtParams       = dend_params.dtParams;                                     % Extract the delta node parameter 
thicknessScale = dend_params.thicknessScale;                               % Extract the dendrite thickness
dtParams(2:3)  = dtParams(2:3)*vres;                                       % Adjust the delta node parameter to the resolution
thicknessScale = thicknessScale*vres*vres;                                 % Adjust the dendrite thickness to the resolution
volsize        = vol_params.vol_sz*vres;                                   % Extract the size of the neural volume

if vol_params.verbose >1                                                   % Optional verbose output
  fprintf('Initializing volume')
end
if(neur_vol_flag)
  neur_vol       = zeros(size(neur_vol),'single');                           % Initialize baseline neural volume
  for i = 1:size(gp_vals,1)
    neur_vol(gp_vals{i,1}) = gp_vals{i,2};                                 % Initialize the new full neural volume to the set neural fluorescences. Do this for somas...
    if i <= size(gp_nuc,1)
      neur_vol(gp_nuc{i,1})  = gp_nuc{i,2};                              % ... and nuclei (if needed)
    end
    if vol_params.verbose >=1                                                % Optional verbose output
      fprintf('.')
    end
  end
end
if vol_params.verbose >1                                                   % Optional verbose output
  fprintf('\n')
end

M            = rand(volsize,'single');                                     %
M(bg_pix==0) = realmax('single');                                          %
M(1,:,:)     = realmax('single');                                          %
M(:,1,:)     = realmax('single');                                          %
M(:,:,1)     = realmax('single');                                          %
M(end,:,:)   = realmax('single');                                          %
M(:,end,:)   = realmax('single');                                          %
M(:,:,end)   = realmax('single');                                          %


if vol_params.verbose >1                                                   % Optional verbose output
  tic
end

if (~isfield(dend_params,'dendVar'))||isempty(dend_params.dendVar)
  dendVar = 0.25;
else
  dendVar = dend_params.dendVar;
end

idxvol = zeros(volsize,'uint16');
numvol = zeros(volsize,'single');

maxlength  = bg_params.maxlength;                                          % Set the maximum dendrite length to 300um
distsc     = bg_params.distsc;                                             % 
fillweight = bg_params.fillweight;                                         % 
maxel      = bg_params.maxel;                                              % 
minlength  = bg_params.minlength;                                          % Set the minimum dendrite length to 20um 
dtSize     = [dtParams(2) dtParams(2) dtParams(3)];                        % 
numpts     = 0;                                                            % 
idx        = 0;                                                            % Initialize the dendrite counter to 0
shiftdist  = 3;                                                            % 
for j = 1:((prod(volsize+2*dtSize)/prod(volsize))-1)*vol_params.N_neur
    dendpts = [];
    root    = floor(rand(1,3).*(volsize+2*dtSize)-dtSize);
    while(root(1)>0&&root(2)>0&&root(3)>0&&root(1)<=volsize(1)&&root(2)<=volsize(2)&&root(3)<=volsize(3))
        root = floor(rand(1,3).*(volsize+2*dtSize)-dtSize);
    end
    neur_locs = cat(1,neur_locs,single(root/vres));
    
    for i = 1:dtParams(1)
      theta = rand*2*pi;                                                   % Get a random orientation
      r     = sqrt(rand)*dtParams(2);                                      % Get a random distance
      dends = floor([r.*cos(theta)+root(1);r.*sin(theta)+root(2); ...
                                      2*dtParams(3)*(rand-0.5)+root(3)])'; % Propose an end-point based on the orientation and distance drawn
      if(dends(1)>0&&dends(2)>0&&dends(3)>0&&dends(1)<=volsize(1)&&dends(2)<=volsize(2)&&dends(3)<=volsize(3))
          [maxShift,shiftLoc] = max([(root<1).*(ones(1,3)-root)./(dends-root) (root>volsize).*(volsize-root)./(dends-root)]);
          bgpts = [];
          numit = 0;
          while(isempty(bgpts) && numit<30)
              numit = numit+1;
              root2 = round(maxShift*(dends-root)+root);
              switch shiftLoc
                  case 1; root2 = root2+[0 randi(shiftdist) randi(shiftdist)];
                  case 2; root2 = root2+[randi(shiftdist) 0 randi(shiftdist)];
                  case 3; root2 = root2+[randi(shiftdist) randi(shiftdist) 0];
                  case 4; root2 = root2+[0 randi(shiftdist) randi(shiftdist)];
                  case 5; root2 = root2+[randi(shiftdist) 0 randi(shiftdist)];
                  case 6; root2 = root2+[randi(shiftdist) randi(shiftdist) 0];            
              end
              if(root2(3)>volsize(3)); root2(3)=volsize(3); end              
              if(root2(2)>volsize(2)); root2(2)=volsize(2); end                
              if(root2(1)>volsize(1)); root2(1)=volsize(1); end                  
              if(root2(3)<1); root2(3)=1; end                        
              if(root2(2)<1); root2(2)=1; end                             
              if(root2(1)<1); root2(1)=1; end  
                  bgpts = dendrite_randomwalk2(M,root2,dends,distsc,maxlength,...
                                                  fillweight,maxel,minlength);
              if(~isempty(bgpts))
                  bgpts = cat(1,root2,bgpts);
                  try
                      dendSz = max(0,normrnd(1,dendVar))^2;
                      if(size(bgpts,1)>2)
                        bgptsW = dendSz*single(1-(1-1/sqrt(2))*[0;sum(abs(diff(abs(diff(bgpts)))),2)/2;0]);
                      else
                        bgptsW = dendSz*ones(size(bgpts,1),1);
                      end
                      bgptsI = sub2ind(volsize,bgpts(:,1),bgpts(:,2),bgpts(:,3));
                      dendpts = cat(1,dendpts,bgptsI);
                      numvol(bgptsI) = bgptsW;
                  catch ME
                      pause(0.5);
                      rethrow(ME)
                  end
              end
          end
      end
    end
    if(~isempty(dendpts))
      idx    = idx+1;
      numpts = numpts+length(dendpts);
      idxvol(dendpts) = idx;
      numvol(dendpts) = numvol(dendpts)*thicknessScale*dtParams(4);         % Adjust thickness of generated dendrites to follow variation of dendrite diameters
    end
end
[~,pathnum]        = dilateDendritePathAll(numvol,idxvol,bg_pix==0);       % Give the dendrites some width
vol_params.N_den2  = idx;                                                  % Save the number of dendrites made
Ncomps             = vol_params.N_neur + vol_params.N_den;                 % Get full number of components (neurons + dendrites)
pathnum(pathnum>0) = pathnum(pathnum>0) + Ncomps;                          % Create indecis that start at the total number of components
neur_num           = neur_num + pathnum;                                   % Add dendrite indices to numberically indexed volume
wtSc               = dend_params.weightScale;

for i = (Ncomps+(1:idx))
    gp_vals{i,1}           = int32(vec(find(neur_num==i)));                % Save memory by only storing the locations...
    gp_vals{i,2}          = (wtSc(2)*exp(-((dtParams(2)/vres)/wtSc(1))) ... 
                             +(1-wtSc(2)))*(1-wtSc(3)*rand(size(gp_vals{i,1})));  % ... and the values    
    gp_vals{i,3}          = false(size(gp_vals{i,1}));
    if(neur_vol_flag)
      neur_vol(gp_vals{i,1}) = gp_vals{i,2};                                 % Add dendrites to fluoresence volume (with constant fluoresence 1 that will change later)
    end
end

if(~neur_vol_flag)
  neur_vol = [];
end
if vol_params.verbose >= 1
  fprintf('done.\n')
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%