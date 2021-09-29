function [pos,pdf] = pseudoRandSample3D(sz, nsamps, width, weight, pdf)

% [pos,pdf] = pseudoRandSample3D(sz, nsamps, width, weight, pdf)
% 
% Function to sample pseudo randomly within a 3D array, with partial
% exclusion of previously sampled locations. The inputs are:
%
%  - sz         - Size of input matrix
%  - nsamps     - Number of locations to sample
%  - width      - Width of Gaussian exclusionary zone
%  - weight     - Weight of Gaussian (0-1, 1 fully excludes previously
%                 sampled locations)
%  - pdf        - Probability density function of allowable locations to
%                 sample
% 
% The outputs are:
%  - pos        - Sampled positions
%  - pdf        - Resulting probability density function of future
%                 positions that may be sampled
% 
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin<5)
  pdf = ones(sz,'single');
end

if(nargin<4)
  weight = 1;
end

if(nargin<3)
  width = 2;
end

[X,Y,Z] = meshgrid(-ceil(2*width):ceil(2*width));
gpdf = single(-weight*exp(-(X.^2+Y.^2+Z.^2)/(width^2)));
pos = zeros(nsamps,3);
i = 1;
while(i<=nsamps)
  rndpt = ceil(rand(1,3).*sz);
  if((pdf(rndpt(1),rndpt(2),rndpt(3))-rand)>0)
    xc = int32([max(0,2*width+1-rndpt(1)) min(0,sz(1)-rndpt(1)-2*width)+4*width]+1);
    yc = int32([max(0,2*width+1-rndpt(2)) min(0,sz(2)-rndpt(2)-2*width)+4*width]+1);
    zc = int32([max(0,2*width+1-rndpt(3)) min(0,sz(3)-rndpt(3)-2*width)+4*width]+1);
    
    xi = int32([max(1,rndpt(1)-2*width) min(sz(1),rndpt(1)+2*width)]);
    yi = int32([max(1,rndpt(2)-2*width) min(sz(2),rndpt(2)+2*width)]);
    zi = int32([max(1,rndpt(3)-2*width) min(sz(3),rndpt(3)+2*width)]);
    pdf(xi(1):xi(2),yi(1):yi(2),zi(1):zi(2)) = pdf(xi(1):xi(2),yi(1):yi(2),zi(1):zi(2))+gpdf(xc(1):xc(2),yc(1):yc(2),zc(1):zc(2));
    pos(i,:) = rndpt;
    i = i+1;
  end
end