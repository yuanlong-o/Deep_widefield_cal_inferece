function Z = zernike(i,x,y)

% function Z = zernike(i, x, y)
%
% Creates the Zernike polynomial with mode index i. The inputs are
%
% - i     - zernike index
% - x     - x position of evaluated points
% - y     - y position of evaluated points
%
% The output is
% - Z     - zernike polynomial profile
%
% Adapted from Numerical Simulation of Optical Wave Propagation by Jason D. Schmidt (2010)
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[theta, r] = cart2pol(x, y);
n = zidx(i,1);
m = zidx(i,2);
if m==0
  Z = sqrt(n+1)*zrf(n,0,r);
else
  if mod(i,2) == 0 
    Z = sqrt(2*(n+1))*zrf(n,m,r) .* cos(m*theta);
  else
    Z = sqrt(2*(n+1))*zrf(n,m,r) .* sin(m*theta);
  end
end
return

% Zernike radial function
function R = zrf(n, m, r)
R = 0;
for s = 0:(n-m)/2
  R = R+(-1)^s*gamma(n-s+1)/(gamma(s+1)*gamma((n+m)/2-s+1)* ... 
    gamma((n-m)/2-s+1))*r.^(n-2*s);
end

% Zernike index function
function out = zidx(i,j)
idx = ceil(sqrt(0.25+2*i)-1.5);
if j==1
  out = idx;
elseif j==2
  out = mod(idx+1,2).*2.*floor((i-(idx+1).*idx/2)/2) + ...
    mod(idx,2).*(2*ceil((i-(idx+1).*idx/2)/2)-1);
end