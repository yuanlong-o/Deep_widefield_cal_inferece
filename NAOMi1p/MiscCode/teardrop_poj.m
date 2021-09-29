function Vtear = teardrop_poj(V, varargin)

% Vtear = teardrop_poj(V)
%
% Project a series of points on a sphere to points on a tear-drop shape.
% The projection for each point is given by
%  
%            [cos(\phi)sin(\theta)sin^m(\theta/2)]
%  Vtear =   [sin(\phi)sin(\theta)sin^m(\theta/2)]
%            [cos(\theta)                        ]
%
% where \phi and \theta are the azymouth and elevation angles,
% respectively. The inputs to this function are:
%   - V        - Nx3 unit norm points (sampled from a unit sphere)
%   - p        - Parameter of the tear-drop shape (default = 1)
%   - plot_opt - OPTIONAL parameter to produce a plot of the resulting
%                shape
%
% Output is
%   - Vtear - The set of points V projected on the unit sphere
%
% 2016 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 1                                                              % Check if the tear-drop parameter p was provided
    p = varargin{1};
else
    p = 1;
end

if nargin > 2                                                              % Check if the tear-drop parameter p was provided
    plot_opt = varargin{2};
else
    plot_opt = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate shape points

rr = sqrt(V(:,1).^2+V(:,2).^2);                                            % Get xy-plane radii
tt = pi - atan(rr./abs(V(:,3))) - (V(:,3)>0)*pi;                           % Create azymouth angles

% Vtear(:,1) = (V(:,1)./rr).*sin(tt).*(sin(0.1*tt).^p);                      % Set x-dimension
% Vtear(:,2) = (V(:,2)./rr).*sin(tt).*(sin(0.1*tt).^p);                      % Set y-dimension
Vtear(:,1) = (V(:,1)./rr).*sin(tt).*(sin(0.5*tt).^p);                      % Set x-dimension
Vtear(:,2) = (V(:,2)./rr).*sin(tt).*(sin(0.5*tt).^p);                      % Set y-dimension
Vtear(:,3) = -cos(tt);                                                     % Set z-dimension

Vtear(isnan(Vtear)) = 0;                                                   % Sometimes NaNs happen: but only when a value should be zero (here at least)

if plot_opt                                                                % Optional plotting for debugging purposes
    scatter3(Vtear(:,1), Vtear(:,2), Vtear(:,3))
    axis equal
    view(3)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%