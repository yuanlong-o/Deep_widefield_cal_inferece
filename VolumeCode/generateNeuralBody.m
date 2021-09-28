function [Vcell,Vnuc,Tri,a] = generateNeuralBody(neur_params)

% [Vcell,Vnuc,Tri,a] = generateNeuralBody(neur_params)
%
% This function samples a neural shape using an Isotropic Gaussian process
% on a sphere. The function basically draws a set of points whos
% covariance matrix depends on the geodescic distance along the sphere.
% Also provided is a sphere embedded in the 
% 
% Inputs to this function are
%   - neur_params - Struct containing parameters for neuron generation
%       .n_samps     - Number of sphere samples to use in creating the mesh
%                      for generating soma and nucleus shapes (default =
%                      200) 
%       .l_scale     - length-scale for the isotropic GP of the soma
%                      shapes. This controls the shape `bumpiness' (default
%                      = 90) 
%       .p_scale     - Overall variance of the isotropic GP of the soma 
%                      shape. (default = 90) 
%       .avg_rad     - Average radius of each neuron in um (default =
%                      6.8 um) 
%       .nuc_fluorsc - Potential fluorescence in the nucleus (default =
%                      0)
%       .min_thic    - Minimum cytoplasmic thickness (default = 0.8)
%       .eccen       - Maximum eccentricity of neuron (default = 0.35)
%       .exts        - Parameters dictating the max/min of the soma radii
%                      (Default = [0.75,1.7])
%       .nexts       - Parameters dictating the extent to shrink and smooth
%                      the nucleus (Default = [60,20])
%       .neur_type   - Option for neuron type (Default 'pyr')
% 
% The outputs of this function are
%   - Vcell:       Nx3 matrix indicating points on the neural surface
%   - Vnuc         Nx3 matrix indicating points on the nucleus surface
%   - Tri          Triangulation matrix indicating the faces of the neural
%                  and nucleus surfaces (same triangulation for both)
%   - a            Vector with the rotation angle of the cell (Rx,Ry,Rz)
% 
% 2016 - Adam Charles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing

if isempty(neur_params)                                                    % Make sure the neur_params struct exists at all
    neur_params.TMP = [];
end
neur_params = check_neur_params(neur_params);                              % Check neuron parameters
pwr       = 1;                                                             % GP power: DON'T CHANGE - VERY SENSITIVE
nucoff    = 3;                                                             % Nucleus offset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get sampling on a sphere

if ((~isfield(neur_params,'S_samp'))||isempty(neur_params.S_samp))&&...
    ((~isfield(neur_params,'Tri'))||isempty(neur_params.Tri))              % Check if a sampling is already provided
    [V,Tri]=SpiralSampleSphere(neur_params.n_samps,false);                 % This function uses a spiral method to sample uniformly on a sphere
else
    V   = neur_params.S_samp;                                              % If provided, use the provided sampling and triangulation
    Tri = neur_params.Tri;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate covariance based on distances

if strcmp(neur_params.neur_type,'pyr')                                     % In the case of pyramidal neurons, add a tear-dropped mean to bias the GP
       Vtear = teardrop_poj(V,1);                                          % Change samples to samples on a teardrop
elseif strcmp(neur_params.neur_type,'peanut')                              % Take me out to the ball-game
       Vtear = teardrop_poj(V,2);                                          % Change samples to samples on a teardrop
else
    % Do nothing
end
if ((~isfield(neur_params,'dists'))||isempty(neur_params.dists))           % If 'dists' is not provided, calculate arc-lengths between points
    dists = bsxfun(@minus, reshape(V, [size(V,1),1,size(V,2)]), ...
                                     reshape(V, [1,size(V,1),size(V,2)])); % Calculate element-wise distance between points
    dists = sqrt(sum(dists.^2,3));                                         % Calculate full distance between points
    dists = 2*asin(dists./(2));                                            % The geodesic distance is the arc-length
    dists = neur_params.p_scale*exp(-(dists./(neur_params.l_scale)).^pwr); % The actual covariance elements are e^{-dist/l}
else
    dists = neur_params.dists;                                             % If provided, use the provided distances
    dists = neur_params.p_scale*exp(-(dists./(neur_params.l_scale)).^pwr); % The actual covariance elements are e^{-dist/l}
end

if ((~isfield(neur_params,'Rtear'))||isempty(neur_params.Rtear))
    if strcmp(neur_params.neur_type,'pyr')                                 % In the case of pyramidal neurons, add a tear-dropped mean to bias the GP
        Rtear = 1*sqrt(sum(Vtear.^2,2));                                   % Get radii of points on the teardrop
    elseif strcmp(neur_params.neur_type,'peanut')                          % Take me out to the ball-game
        Vtear = Vs;                                                        % Initialize points to the points on a sphere
        Rtear = 1*sqrt(sum(Vtear.^2,2));                                   % Get radii of points on the teardrop
    else 
        Rtear = 1;
    end
else
    Rtear = neur_params.Rtear;                                             % If provided, use the provided means
end
% eo.issym  = 1;                                                             % Use symmetric flag for eigenvalue finding 
eo.isreal = 1;                                                             % Use real flag for eigenvalue finding
min_eig   = 1.03*eigs(dists,1,'sa',eo);                                    % To combat ill conditioning from numerical errors, find the minimum eigenvalue
if min_eig < 0
    dists = dists + abs(min_eig)*eye(size(dists));                         % Find the value that makes sure the covariance is PSD
end
x_bnds = neur_params.exts*neur_params.avg_rad;                             % Set bounds for the maximum and minimum radii of the soma shape
x_base = abs(mvnrnd(0*Rtear,dists).');                                     % Sample from a Gaussian with the correct covariance matrix
x      = x_base - mean(x_base) + neur_params.avg_rad*Rtear;                % Re-center the samples 
xmin   = min(min(x),x_bnds(1));                                            % Calculate the minimum radius of the current shape
x      = (x_bnds(2) - x_bnds(1))*(x - xmin)...
                             /(max(max(x),x_bnds(2)) - xmin) + x_bnds(1);  % Renormalize the radii so that the cell isn't ever too big or too small in any direction
if strcmp(neur_params.neur_type,'pyr') 
    x2   = x_base - mean(x_base) + neur_params.avg_rad;                    % Re-center the samples (now for the nucleus)
    xmin = min(min(x2),x_bnds(1));                                         % Calculate the minimum radius of the current shape (now for the nucleus)
    x2   = (x_bnds(2) - x_bnds(1))*(x2 - xmin)...
                             /(max(max(x2),x_bnds(2)) - xmin) + x_bnds(1); % Renormalize the radii so that the nucleus isn't ever too big or too small in any direction (now for the nucleus)
else
    x2 = x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Elliptical Nucleus

% eccens = ones(1,3)+neur_params.eccen;                      % Get ellipse eccentricities
eccens = ones(1,3)+neur_params.eccen.*(rand(1,3)-[0.5 0.5 0]);                      % Get ellipse eccentricities
% eccens = ones(1,3)+neur_params.eccen*(rand(1,3)-0.5);                      % Get ellipse eccentricities
% eccens = eccens./prod(eccens);
eccens = eccens/(prod(eccens)^(1/3));

if strcmp(neur_params.neur_type,'pyr')  
    Vetear = bsxfun(@times,Vtear,eccens);                                  % Create an elliptical teardrop
    Vetear = Vetear/sqrt(mean(sum(Vetear.^2,2)));                          % Normalize the elliptical teardrop
else
    Vetear = bsxfun(@times,V,eccens);                                      % Create an elliptical shape
    Vetear = Vetear/sqrt(mean(sum(V.^2,2)));                               % Normalize the elliptical shape
end

Vcell   = bsxfun(@times, Vetear, x(:));                                    % Multiply through the sampled points to get the samples on the neural surface
Vcell   = bsxfun(@plus, Vcell, [0, 0, -nucoff]);                           % Shift the soma down (for easier modulation of the nucleus)
Vnorms  = sqrt(sum(Vcell.^2,2));                                           % Get soma surface point extents
% Vnuc    = bsxfun(@times, bsxfun(@times, V, [1,1,-1].*eccens), x2(:));      % Initialize the nucleus to an off-set cell shape
% Vnorms2 = sqrt(sum(Vnuc.^2,2));                                            % Get soma surface point extents
% Vnorms2 = neur_params.nexts(2)*(neur_params.nexts(1)*(Vnorms2 ...
%                 - min(Vnorms2)) + (1-neur_params.nexts(1))*max(Vnorms2));  % Shrink & smooth the nucleus
% Vnorms2 = Vnorms2 + min(Vnorms - Vnorms2) - neur_params.min_thic(1);          % Make sure the nucleus fits inside the soma with the minimum thickness
% Vnuc    = bsxfun(@times,Vnuc,Vnorms2./ sqrt(sum(Vnuc.^2,2)));              % Apply shrinkage values to soma surface points

Vnuc    = bsxfun(@times, bsxfun(@times, V, [1,1,-1]), x2(:));      % Initialize the nucleus to an off-set cell shape
Vnorms2 = sqrt(sum(Vnuc.^2,2));                                            % Get soma surface point extents
Vnorms2 = neur_params.nexts(2)*(neur_params.nexts(1)*(Vnorms2 ...
                - min(Vnorms2)) + (1-neur_params.nexts(1))*max(Vnorms2));  % Shrink & smooth the nucleus
Vnorms2 = Vnorms2 + min(Vnorms - Vnorms2) - neur_params.min_thic(1);          % Make sure the nucleus fits inside the soma with the minimum thickness
Vnuc    = bsxfun(@times, bsxfun(@times, Vnuc,eccens),Vnorms2./ sqrt(sum(Vnuc.^2,2)));              % Apply shrinkage values to soma surface points

lat_ang = rand(1)*2*pi;
lat_shft = (1-abs(rand(1)-rand(1)))*neur_params.min_thic(2)*[sin(lat_ang), cos(lat_ang)];
% lat_shft = rand(1)*neur_params.min_thic(2)*[sin(lat_ang), cos(lat_ang)];

Vcell   = bsxfun(@plus, Vcell, [0, 0, nucoff]);                            % Shift the soma up
Vnuc    = bsxfun(@plus, Vnuc, [lat_shft(1), lat_shft(2), nucoff]);                             % Shift the nucleus up

% Vcell   = bsxfun(@plus, Vcell, [0, 0, nucoff]);                            % Shift the soma up
% Vnuc    = bsxfun(@plus, Vnuc, [0, 0, nucoff]);                             % Shift the nucleus up

% Optional scaling of nucleus size to adjust to normalized radius
if(isfield(neur_params,'nuc_rad')&&~isempty(neur_params.nuc_rad))
  [~,VnucSz] = convhull(Vnuc(:,1),Vnuc(:,2),Vnuc(:,3));
  nucsz = (4/3)*pi*(neur_params.nuc_rad(1).^3);
  if(length(neur_params.nuc_rad)>1)
    Vnuc = Vnuc*(((nucsz/VnucSz)^(1/3))^(1/neur_params.nuc_rad(2)));
  else
    Vnuc = Vnuc*((nucsz/VnucSz)^(1/3));
  end
end


if ((~isfield(neur_params,'max_ang'))||isempty(neur_params.max_ang))
  max_ang = 20;
else
  max_ang = neur_params.max_ang;
end
a     = -abs(max_ang)+2*abs(max_ang)*rand(1,3);                            % Choose a random rotation angle
Rx    = [1 0 0; 0 cosd(a(1)) -sind(a(1)); 0 sind(a(1)) cosd(a(1))];        % ... Rotation around x
Ry    = [cosd(a(2)) 0 sind(a(2)); 0 1 0; -sind(a(2)) 0 cosd(a(2))];        % ... Rotation around y
Rz    = [cosd(a(3)) -sind(a(3)) 0; sind(a(3)) cosd(a(3)) 0; 0 0 1];        % ... Rotation around z
Vnuc  = Vnuc*Rx*Ry*Rz;                                                     % Apply the three rotation matrices to the nucleus shape 
Vcell = Vcell*Rx*Ry*Rz;                                                    % Apply the three rotation matrices to the soma shape
                              
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                         