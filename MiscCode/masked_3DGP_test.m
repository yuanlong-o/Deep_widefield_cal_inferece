function [gp_tmp,varargout] = masked_3DGP_test(grid_sz, l_scale, p_scale, mu, bin_mask, threshold, l_weights,type)

% [gp_vals] = masked_3DGP_v2(grid_sz, l_scale, p_scale, varargin)
% 
% This function draws from a 3-D GP with mean mu, length-scale l_scale and
% variance scale parameter p_scale. Basically it returns a volume X, where
%             X ~ N(mu, C), where C_{i,j} = p*e^{-(i-j)^2/(2*l^2)}
% The method used is to draw i.i.d. random variables in the frequency
% domain, apply the fft of C_{i,j} point-wise, and ise in ifft to return to
% the spatial domain. This function was modified for speedups and multiple
% spatial scales
% 
% It takes in parameters 
%     - grid_sz:  The dimensions of the sample to take (if grid_sz is a
%                 scalar, then the samples in each dimension will be that
%                 value).
%     - l_scale:  The length scale of the GP (i.e. l in the above
%                 definition of C_{i,j})
%     - p_scale:  The covariance scaling of the GP (i.e. p in the above
%                 definition of C_{i,j})
%     - mu:       The mean of the grid size. Can be either a scalar or a
%                 matrix of the same size as X.
%     - bin_mask: OPTIONAL Binary mask to set certain values of the output
%                 to zero
% 
% And returns the output
%     - X:        The sample from the GP, as defined above.
%
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Inputs

if nargin < 8 || isempty(type)
    type = 'single';
end
if nargin < 7 || isempty(l_weights)
    l_weights = 1;
end
if nargin < 6 || isempty(threshold)
    threshold = 1e-10;
end
if nargin < 5 || isempty(bin_mask)
    bin_mask = 1;
end
if numel(grid_sz) == 1
    grid_sz = grid_sz*[1,1,1];
end
if size(l_scale,2) == 1
    l_scale = l_scale*[1,1,1];
end
if numel(l_weights) == 1
    l_weights = ones(size(l_scale,1),1)*l_weights;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Kernel

wmx = pi/2;
if(strcmp(type,'single'))
  grid_x = reshape(single(linspace(-wmx,wmx,grid_sz(1))).^2,[],1,1);
  grid_y = reshape(single(linspace(-wmx,wmx,grid_sz(2))).^2,1,[],1);
  grid_z = reshape(single(linspace(-wmx,wmx,grid_sz(3))).^2,1,1,[]);
elseif(strcmp(type,'double'))
  grid_x = reshape(linspace(-wmx,wmx,grid_sz(1)).^2,[],1,1);
  grid_y = reshape(linspace(-wmx,wmx,grid_sz(2)).^2,1,[],1);
  grid_z = reshape(linspace(-wmx,wmx,grid_sz(3)).^2,1,1,[]);
end  
  
  % gp_vals = zeros(grid_sz,'single');
gp_tmp = zeros(grid_sz,type);

for i = 1:size(l_scale,1)
    ker_x = exp(-grid_x*l_scale(i,1).^2);
    ker_y = exp(-grid_y*l_scale(i,2).^2);
    ker_z = exp(-grid_z*l_scale(i,3).^2);
    gp_tmp = gp_tmp+prod(l_scale(i,:))*(bsxfun(@times,bsxfun(@times,ker_x,ker_y),ker_z).^2)*l_weights(i).^2;
end
gp_tmp = sqrt(gp_tmp);

if nargout > 1
    varargout{1} = gp_tmp;
end

sprev = rng;
rng('shuffle','simdTwister');
for i = 1:size(gp_tmp,3)
  gp_tmp(:,:,i) = (randn(grid_sz(1:2),type)+1j*randn(grid_sz(1:2),type)).*gp_tmp(:,:,i);
end
% gp_tmp = exp(1j*2*pi*rand(grid_sz,type)).*gp_tmp;
% gp_tmp = single(sqrt(chi2rnd(2,grid_sz)).*gp_tmp);
rng(sprev);

gp_tmp = 2*p_scale*bin_mask.*sqrt(numel(gp_tmp))*real(ifftshift(ifftn(ifftshift(gp_tmp))))+mu; % Move random process back to the spatial domain

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     ker_1 = bsxfun(@times,bsxfun(@times,ker_x,ker_y),ker_z);
%     all_ker = all_ker+prod(l_scale(i,:))*ker_1.^2*l_weights(i).^2;

% gp_vals = randn(grid_sz,type)+1j*randn(grid_sz,type);              % First make a random array
% gp_vals = gp_vals.*all_ker;

% gp_vals = p_scale*(2^4.5/pi^1.5)*bin_mask.*gp_vals/sqrt(length(l_weights))+mu;  % Apply a binary mask as needed