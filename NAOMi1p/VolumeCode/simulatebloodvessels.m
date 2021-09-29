function [neur_ves,vasc_params,neur_ves_all] = simulatebloodvessels(vol_params,vasc_params)

% [neur_ves,vasc_params] = simulatebloodvessels(vol_params,vasc_params)
%
% Function to simulate vasculature in and above a neural volume. The inputs
% to this function are:
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
%   - vasc_params - Struct containing parameters for vasculature simulation
%       .ves_shift       - 3-vector of the amount of wobble allowed for
%                          blood vessels (default = [5 15 5] um)
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
%                          (default = 0.5) 
%       .vesSize         - vessel radius (surface, axial, capillaries) in
%                          um (default = [15 6 2] um)
%       .vesFreq         - blood vessel freqency in microns (default = 
%                          [125 200 50] um) 
%       .vesNumScale     - blood vessel number random scaling factor
%                          (default = 0.2) 
%       .sourceFreq      - Rate of generation of source nodes (default =
%                          1000 um/node) 
%       .sepweight       - Set weight that guides how far, on average, the
%                          vasculature nodes are placed (value is between 0
%                          and 1; default = 0.75)
%       .distsc          - How strongly local capillary connections are.
%                          Higher numbers indicate more local connections
%                          (default = 4)
%       .node_params     - Set of parameters to guide node placement:
%          .maxit        - Maximum iteration to place nodes (default = 25)
%          .lensc        - Average distance between vasculature branch
%                          points (default = 50 um) 
%          .varsc        - Standard deviation of distances between
%                          vascualure branch points (default = 15 um) 
%          .mindist      - Minimum inter-node distance (default = 10 um)
%          .varpos       - Standard deviation of vasculature placement
%                          (default = 5 um) 
%          .dirvar       - The maximum branching angle (default = pi/8)
%          .branchp      - Probability of branching surface vasculature
%                          (default = 0.02) 
%          .vesrad       - Radius of surface vasculature (default = 25 um)
%
% The outputs of this function are
%   - neur_ves     - Array of the blood vessel locations
%   - vasc_params  - Updated scruct of parameters (with defaults filled in
%                    if needed)
%
% 2017 - Alex Song and Adam Charles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

vol_params  = check_vol_params(vol_params);                                % Check volume parameters
vasc_params = check_vasc_params(vasc_params);                              % Check vasculature parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulating Blood Vessels

if vol_params.verbose == 1
  fprintf('Generating in-volume blood vessels...');
elseif vol_params.verbose >1
  fprintf('Generating in-volume blood vessels...\n');
  tic
end

% Setup parameters for simulation
vres       = vol_params.vres;                                              % Pull out volume resolution

vp = vasc_params;
vp.depth_surf = vasc_params.depth_surf*vres;                                  % Calculate depth beneath the surface
vp.mindists    = vasc_params.vesFreq*vres/2;                                  % Set minimum distance between nodes
vp.maxcappdist = 2*vasc_params.vesFreq(3)*vres;                               % Maximum capilliary distance
vp.vesSize = vasc_params.vesSize*vres;                                        % Calculate the blood-vessel size

node_params = vasc_params.node_params;
np = node_params;
np.lensc   = node_params.lensc*vres;
np.varsc   = node_params.varsc*vres;
np.mindist = node_params.mindist*vres;
np.varpos  = node_params.varpos*vres;
np.vesrad  = node_params.vesrad*vres;

if (~isfield(vol_params,'vasc_sz'))||isempty(vol_params.vasc_sz)||nargout<3
  nv.vol_sz  = vol_params.vol_sz+[0 0 1]*vol_params.vol_depth;             % Extract volume size
else
  nv.vol_sz  = vol_params.vasc_sz;
end
nv.size = nv.vol_sz*vres;                                                   %
nv.szum = nv.vol_sz;                                                   %
nv.nsource = max(round((2*(nv.vol_sz(1)+nv.vol_sz(2))/(vp.sourceFreq))*...
                                 abs(1+vp.vesNumScale*randn)),0); % 
nv.nvert = max(round((nv.vol_sz(1)*nv.vol_sz(2)/(vp.vesFreq(2)^2))*...
                                 abs(1+vp.vesNumScale*randn)),0); %
nv.nsurf = max(round((nv.vol_sz(1)*nv.vol_sz(2)/(vp.vesFreq(1)^2))*...
                                 abs(1+vp.vesNumScale*randn)),0); % 
nv.ncapp = max(round((prod(nv.vol_sz)/(vp.vesFreq(3)^3))*...
                                 abs(1+vp.vesNumScale*randn)),0); %
                               
%% Initialize a few points for vertical vessels, Initialize some points in surface for surface vessels
[nodes,nv] = growMajorVessels(nv,np,vp);

%% Convert node structure to a connection structure
conn = nodesToConn(nodes);
nv.nconn = length(conn);

% Shift surface vessel location to adjust for vessel diameter
for i = 1:length(conn)
  if(sum(strcmp(nodes(conn(i).start).type,{'edge','surf','sfvt'})))
    nodes(conn(i).start).pos(3) = min(nodes(conn(i).start).pos(3)+ ...
      ceil(conn(i).weight/length(nodes(conn(i).start).conn)),nv.size(3));
  end
  if(sum(strcmp(nodes(conn(i).ends).type,{'edge','surf','sfvt'})))
    nodes(conn(i).ends).pos(3) = min(nodes(conn(i).ends).pos(3)+ ...
      ceil(conn(i).weight/length(nodes(conn(i).start).conn)),nv.size(3));
  end
end

%% Create initial volume with major vessels
[neur_ves,conn] = connToVol(nodes,conn,nv);

%% Initialize and connect capillaries
[nodes,conn,nv] = growCapillaries(nodes,conn,neur_ves,nv,vp,vres);

%% Add capillaries to rest of volume
cappidxs = find(cellfun(@isempty,{conn.locs}));
[neur_ves,~] = connToVol(nodes,conn,nv,cappidxs,neur_ves);

if isfield(vol_params,'vasc_sz')&&~isempty(vol_params.vasc_sz)&&nargout==3
  neur_ves_all = neur_ves;
  sz = [vol_params.vol_sz(1:2), vol_params.vol_depth+vol_params.vol_sz(3)]*vres;
  sz_diff  = ceil((vol_params.vasc_sz*vres-sz)/2);
  neur_ves = neur_ves((1:sz(1))+sz_diff(1),(1:sz(2))+sz_diff(2),(1:sz(3))+sz_diff(3));
elseif(nargout==3)
  neur_ves_all = neur_ves;
end
%%
if vol_params.verbose == 1
  fprintf('done.\n');
elseif vol_params.verbose >1
  Tdone = toc;
  fprintf('done (%f seconds).\n',Tdone);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
