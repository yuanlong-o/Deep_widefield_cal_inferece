function [psf,x,y,z] = gaussian_psf(psflen,lambda,sampling,matSize,theta)
% length is the FWHM (1/e2 over 1.7) of the PSF, lambda is the wavelength,
% sampling is the spacing between pixels, and matSize is the output psf
% size. psf is the output, x y and z are the coordinate positions
% for FWHM axial length, we describe this as where the signal for a
% particular plane is 1/2 of the central plane

if nargin<5
  theta = 0;
end

if nargin<4
  matSize = [101 101 1001];
end

if nargin<3
  sampling = [0.1 0.1 0.1];
end

if length(sampling)<3
  sampling = [sampling(1) sampling(1) sampling(1)];
end

if nargin<2
  lambda = 0.92;
end

if nargin<1
  psflen = 10;
end

n = 1;
% zr = 1/(sqrt(sqrt(2)-1))*psflen/2;
zr = psflen/2;

x = ((1:matSize(1))-round(matSize(1)/2))*sampling(1);
y = ((1:matSize(2))-round(matSize(2)/2))*sampling(2);
z = ((1:matSize(3))-round(matSize(3)/2))*sampling(3);
[xg,yg,zg] = meshgrid(x,y,z);

R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
xg2 = R(1,1)*xg+R(1,2)*zg;
zg2 = R(2,1)*xg+R(2,2)*zg;
xg = xg2;
zg = zg2;

psf = (exp(-2*pi*n*(xg.^2+yg.^2)./(zr*lambda*(1+(zg/zr).^2)))./(1+(zg/zr).^2)).^2;
