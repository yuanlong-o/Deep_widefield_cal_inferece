function [psf,x,y,z,intensity] = gaussian_psf_na(na,lambda,sampling,matSize,theta,nidx)

% [psf,x,y,z] = gaussian_psf(psflen,lambda,sampling,matSize,theta)
% 
% This function provides a sampled point-spread function (PSF) for use in
% the scanning simulation. Specifically, this function generates Gaussian
% shaped PSFs given the following inputs:
%   psflen   - the length of the PSF given as the full-width half-max FWHM
%              (1/e2 over 1.7). 
%   lambda   - the wavelength of the PSF
%   sampling - the spacing between pixels
%   matSize  - the output psf size
% 
% The outputs for this function are:
%   psf - the output
%   x y and z - the coordinate positions for FWHM axial length, defined
%               as where the signal for a particular plane is 1/2 of the
%               central plane 
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing

if nargin<6
  nidx = 1.33;                                                             % index of refraction
end

if nargin<5
  theta = 0;                                                               % Default angle of the beam is zero (axially aligned)
end

if nargin<4
  matSize = [101 101 1001];                                                % Default matrix size is 101x101x101 (units determined by "sampling"
end

if nargin<3
  sampling = [0.1 0.1 0.1];                                                % Default sampling is 0.1um in each dimension
end

if length(sampling)<3
  sampling = [sampling(1) sampling(1) sampling(1)];
end

if nargin<2
  lambda = 0.92;                                                           % Default wavelength is 0.92
end

if nargin<1
  na = 0.6;                                                                % Default PSF NA
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the PSF

psflen = 0.626*lambda/(nidx-sqrt(nidx^2-na^2));
zr = psflen/2;                                                             % Cut the PSF length in two for calculation purposes

x = ((0:matSize(1)-1)-round(matSize(1)/2))*sampling(1);                      % Get x grid locations to evaluate on
y = ((0:matSize(2)-1)-round(matSize(2)/2))*sampling(2);                      % Get y grid locations to evaluate on
z = ((0:matSize(3)-1)-round(matSize(3)/2))*sampling(3);                      % Get z grid locations to evaluate on
[xg,yg,zg] = meshgrid(x,y,z);                                              % Generate the meshgrid to evaluate the PSF values on

R   = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
xg2 = R(1,1)*xg+R(1,2)*zg;
zg2 = R(2,1)*xg+R(2,2)*zg;
xg  = xg2;
zg  = zg2;
% intensity = (exp(-2*2*sqrt(log(2))*pi*nidx*(xg.^2+yg.^2)./(zr*lambda*(1+(zg/zr).^2)))...
%                                                   ./(1+(zg/zr).^2));    % Evaluate the PSF values at each grid point
intensity = (exp(-2*pi*nidx*(xg.^2+yg.^2)./(zr*lambda*(1+(zg/zr).^2)))...
                                                  ./(1+(zg/zr).^2));    % Evaluate the PSF values at each grid point
psf = intensity.^2;
psf = permute(psf,[2 1 3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%