function [gp_vals] = masked_3DGP_v2(grid_sz, l_scale, p_scale, mu, bin_mask, threshold, l_weights)

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
grid_x = reshape(single(linspace(-wmx,wmx,grid_sz(1))).^2,[],1,1);
grid_y = reshape(single(linspace(-wmx,wmx,grid_sz(2))).^2,1,[],1);
grid_z = reshape(single(linspace(-wmx,wmx,grid_sz(3))).^2,1,1,[]);
gp_vals = zeros(grid_sz,'single');
for i = 1:size(l_scale,1)
    ker_x = exp(-grid_x*l_scale(i,1).^2);
    ker_y = exp(-grid_y*l_scale(i,2).^2);
    ker_z = exp(-grid_z*l_scale(i,3).^2);
    ker_1 = bsxfun(@times,bsxfun(@times,ker_x,ker_y),ker_z);
    ker_loc = (ker_1(:)>threshold);
    ker_len = sum(ker_loc);
    if(ker_len<prod(grid_sz)/2)
        sprev = rng;
        rng('shuffle','simdTwister');
        TMP = randn([ker_len 1],'single')+1j*randn([ker_len 1],'single');  % First make a random array
        TMP = (l_weights(i)*sqrt(prod(l_scale(i,:))))*TMP.*ker_1(ker_loc);
        rng(sprev);
        gp_vals(ker_loc) = gp_vals(ker_loc)+TMP;                           % Multiply with a kernel (Frequency-domain filtering) 
        clear TMP ker_loc
    else
        clear ker_loc
        sprev = rng;
        rng('shuffle','simdTwister');
        TMP = randn(grid_sz,'single')+1j*randn(grid_sz,'single');          % First make a random array
        TMP = (l_weights(i)*sqrt(prod(l_scale(i,:))))*TMP.*ker_1;
        rng(sprev);
        gp_vals = gp_vals+TMP;                                             % Multiply with a kernel (Frequency-domain filtering) 
        clear TMP
    end
end
gp_vals = sqrt(numel(gp_vals))*real(ifftshift(ifftn(ifftshift(gp_vals)))); % Move random process back to the spatial domain
gp_vals = p_scale*(2^4.5/pi^1.5)*bin_mask.* ...
                                       gp_vals/sqrt(length(l_weights))+mu; % Apply a binary mask as needed

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
