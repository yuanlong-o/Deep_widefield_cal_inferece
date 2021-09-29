function [Uout, UoutAll] = fresnel_propagation_multi(Uin, lambda, dx, z, phi, nidx)

% [Uout, UoutAll] = fresnel_propagation_multi(Uin, lambda, dx, z, phi, nidx)
%
% Simulates optical propagation through a medium with phase masks applied
% at specified points. The inputs are
%
%  - Uin          - incident field
%  - lambda       - wavelength [m]
%  - dx           - lateral pixel width at each step of propagatio [m]
%  - z            - position of each step of propagation [m]
%  - phi          - phaseshift at each step of propagation [rad]
%  - nidx         - index of refraction
%
% The outputs are
%  - Uout         - output field at end of propgation
%  - UoutAll      - output field at all steps of propagation
%
% Adapted from Numerical Simulation of Optical Wave Propagation by Jason D.
% Schmidt (2010)
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda   = lambda/nidx;
sz       = size(Uin);                                                           
[nx, ny] = meshgrid(-sz(1)/2:sz(1)/2-1,-sz(2)/2:sz(2)/2-1);
nx = nx';
ny = ny';
k   = 2*pi/lambda;                                                          
n   = length(z);                                                           
df  = 1./(sz(1)*dx);
dz  = z(2:n)-z(1:n-1);                                                   
sc  = dx(2:n)./dx(1:n-1);                                                
Uin = Uin.*exp(1i*k*((nx*dx(1)).^2 + (ny*dx(1)).^2)*(1-sc(1))/(2*dz(1))).*phi(:,:,1);

if nargout>1
  UoutAll        = zeros(sz(1),sz(2),n,'like',Uin);
  UoutAll(:,:,1) = Uin;
end
tol = 1e-12;
Q = exp(-1i*pi*lambda*dz(1)*((nx*df(1)).^2+(ny*df(1)).^2)/sc(1));
for i = 1:n-1
  if(i>1 && ~(abs(dz(i)-dz(i-1))<tol && abs(df(i)-df(i-1))<tol && all(abs(sc(i)-sc(i-1))<tol)))
    Q = exp(-1i*pi*lambda*dz(i)*((nx*df(i)).^2+(ny*df(i)).^2)/sc(i));
  end
    
  if(i==n-1 && size(phi,3)<n)
    Uin = ifftshift(ifft2(ifftshift(Q.*fftshift(fft2(fftshift(Uin/sc(i)))))));
  else
    Uin = phi(:,:,i+1).*ifftshift(ifft2(ifftshift(Q.*fftshift(fft2(fftshift(Uin/sc(i)))))));
  end
  if(nargout>1)
    UoutAll(:,:,i+1) = Uin;
  end
end
Uout = Uin.*exp(1i*k/2*(sc(n-1)-1)/(sc(n-1)*dz(i))*((nx*dx(n)).^2+(ny*dx(n)).^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
