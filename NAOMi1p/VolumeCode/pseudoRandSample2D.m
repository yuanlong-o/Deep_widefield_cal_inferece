function [pos,pdf] = pseudoRandSample2D(sz, nsamps, width, weight, pdf, maxit)

% [pos,pdf] = pseudoRandSample2D(sz, nsamps, width, weight, pdf, maxit)
% 
% Function to sample pseudo randomly within a 2D matrix, with partial
% exclusion of previously sampled locations. The inputs are:
%
%  - sz         - Size of input matrix
%  - nsamps     - Number of locations to sample
%  - width      - Width of Gaussian exclusionary zone
%  - weight     - Weight of Gaussian (0-1, 1 fully excludes previously
%                 sampled locations)
%  - pdf        - Probability density function of allowable locations to
%                 sample
%  - maxit      - Maximum number of iterations to look for new positions
%
% The outputs are:
%  - pos        - Sampled positions
%  - pdf        - Resulting probability density function of future
%                 positions that may be sampled
% 
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if(nargin<6)
  maxit = 1000;
end

if(nargin<5)
  pdf = ones(sz,'single');
end

if(nargin<4)
  weight = 1;
end

if(nargin<3)
  width = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main function

[X,Y] = meshgrid(-ceil(2*width):ceil(2*width));
gpdf = single(-weight*exp(-(X.^2+Y.^2)/(width^2)));
pos = zeros(nsamps,2);
i = 1;
numit = 0;
while(i<=nsamps && numit<maxit)
  numit = numit+1;
  rndpt = ceil(rand(1,2).*sz);
  if((pdf(rndpt(1),rndpt(2))-rand)>0)
    xc = int32([max(0,2*width+1-rndpt(1)) min(0,sz(1)-rndpt(1)-2*width)+4*width]+1);
    yc = int32([max(0,2*width+1-rndpt(2)) min(0,sz(2)-rndpt(2)-2*width)+4*width]+1);
    
    xi = int32([max(1,rndpt(1)-2*width) min(sz(1),rndpt(1)+2*width)]);
    yi = int32([max(1,rndpt(2)-2*width) min(sz(2),rndpt(2)+2*width)]);
    pdf(xi(1):xi(2),yi(1):yi(2)) = pdf(xi(1):xi(2),yi(1):yi(2))+gpdf(xc(1):xc(2),yc(1):yc(2));
    pos(i,:) = rndpt;
    i = i+1;
    numit = 0;
  end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
