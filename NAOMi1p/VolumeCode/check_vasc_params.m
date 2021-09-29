function vasc_params = check_vasc_params(vasc_params)

% vasc_params = check_vasc_params(vasc_params)
%  
% This function checks the elements of the struct vasc_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - vasc_params - Struct containing parameters for vasculature simulation
%       .flag            - On/off flag for  vasculature simulation
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
%                          um (default = [15 10 3] um) 
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
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(vasc_params)                                                    % Make sure that vasc_params is a struct
    clear vasc_params
    vasc_params = struct;
end

dParams.flag            = 1;                                               % amount of wobble allowed for blood vessels
dParams.ves_shift       = [5 15 5];                                        % amount of wobble allowed for blood vessels
dParams.depth_vasc      = 200;                                             % depth into tissue for which vasculature is simulated
dParams.depth_surf      = 15;                                              % depth into tissue of surface vasculature
dParams.distWeightScale = 2;                                               % scaling factor for weight of node distance scaling
dParams.randWeightScale = 0.1;                                             % scaling factor for weight of nodes (additional variability)
dParams.cappAmpScale    = 0.5;                                             % scaling factor for weights (capillaries, lateral)
dParams.cappAmpZscale   = 0.5;                                             % scaling factor for weights (capillaries, axial)
dParams.vesSize         = [15 9 2];                                        % vessel radius (surface, axial, capillaries) in um
dParams.vesFreq         = [125 200 50];                                    % vessel freqency in microns
dParams.sourceFreq      = 1000;                                            % 
dParams.vesNumScale     = 0.2;                                             % vessel number random scaling factor
dParams.sepweight       = 0.75;                                            %
dParams.distsc          = 4;                                               %

vasc_params = setParams(dParams,vasc_params);

if (~isfield(vasc_params,'node_params'))||isempty(vasc_params.node_params) % Check if the node parameter sub-struct exists. If not - make it
  vasc_params.node_params.maxit   = 25;                                                % Maximum iteration to place nodes
  vasc_params.node_params.lensc   = 50;                                                % 
  vasc_params.node_params.varsc   = 15;                                                % 
  vasc_params.node_params.mindist = 10;                                                % Minimum inter-node distance
  vasc_params.node_params.varpos  = 5;                                                 % 
  vasc_params.node_params.dirvar  = pi/8;                                              % 
  vasc_params.node_params.branchp = 0.02;                                              % 
  vasc_params.node_params.vesrad  = 25;                                                % 
else                                                                       % Check if the node parameter sub-struct exists
    if (~isfield(vasc_params.node_params,'maxit'))||isempty(vasc_params.node_params.maxit)
        vasc_params.node_params.maxit = 25;                                        
    end
    if (~isfield(vasc_params.node_params,'lensc'))||isempty(vasc_params.node_params.lensc)
        vasc_params.node_params.lensc = 50;                                        
    end
    if (~isfield(vasc_params.node_params,'varsc'))||isempty(vasc_params.node_params.varsc)
        vasc_params.node_params.varsc = 15;                                        
    end
    if (~isfield(vasc_params.node_params,'mindist'))||isempty(vasc_params.node_params.mindist)
        vasc_params.node_params.mindist = 10;                                        
    end
    if (~isfield(vasc_params.node_params,'varpos'))||isempty(vasc_params.node_params.varpos)
        vasc_params.node_params.varpos = 5;                                        
    end
    if (~isfield(vasc_params.node_params,'dirvar'))||isempty(vasc_params.node_params.dirvar)
        vasc_params.node_params.dirvar = pi/8;                                        
    end
    if (~isfield(vasc_params.node_params,'branchp'))||isempty(vasc_params.node_params.branchp)
        vasc_params.node_params.branchp = 0.02;                                        
    end
    if (~isfield(vasc_params.node_params,'vesrad'))||isempty(vasc_params.node_params.vesrad)
        vasc_params.node_params.vesrad = 25;                                        
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
