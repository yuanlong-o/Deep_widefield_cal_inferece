function Uout = generateGaussianProfile(X,Y,rad,aper,k,fl,offset)

% Uout = generateGaussianProfile(X,Y,rad,aper,k,fl,offset)
%
% This function generates a Gaussian back aperture intensity profile with a
% fixed aperture size with a focused phase fl. The inputs are
%
% - X           - lateral position X [m]
% - Y           - lateral position X [m]
% - rad         - radius [m]
% - aper        - aperture distance [m]
% - k           - optical wavenumber [rad/m]
% - fl          - focal lenth of lens [m]
% - offset      - offset of Gaussian position [m]
%
% The output is
% - Uout        - scalar field for a gaussian beam with apodization 1
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<7)
  offset = [0 0];
end
rho2 = X.^2+Y.^2;
Uout = exp(-((X-offset(1)).^2+(Y-offset(2)).^2)/(rad^2));                  % Gaussian intensity profile
Uout = Uout.*(rho2<(aper^2));                                              % Fixed aperture
Uout = Uout.*exp(-1i*k/(2*fl)*rho2);                                       % Apply ideal phase
