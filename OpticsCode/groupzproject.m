function out = groupzproject(images,groupsize,type)

% out = groupzproject(images,groupsize,type)
%
% Function to project a 3D matrix along the third dimension in groups of
% the specified size. The inputs are
%
%   - images      - 3D matrix to project in groups
%   - groupsize   - Size of projection group. Needs to be a factor of
%                   size(images,3)
%   - type        - One of 'sum' 'prod' 'mean' 'median' 'mode' 'max' 'min' 
%                   for projection type. Default is 'mean'
%
% The output is
%   - out         - group projected 3D matrix
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3; type = 'mean'; end                                          % If the projection type is unspecified, use the mean

if(rem(size(images,3),groupsize)==0)
    switch type
        case 'sum'
          out = squeeze(sum(reshape(images,size(images,1),size(images,2), ...
            groupsize,size(images,3)/groupsize),3));
        case 'prod'
          out = squeeze(prod(reshape(images,size(images,1),size(images,2), ...
            groupsize,size(images,3)/groupsize),3));
        case 'mean'
          out = squeeze(mean(reshape(images,size(images,1),size(images,2), ...
            groupsize,size(images,3)/groupsize),3));
        case 'median'
          out = squeeze(median(reshape(images,size(images,1),size(images,2), ...
            groupsize,size(images,3)/groupsize),3));
        case 'mode'
          out = squeeze(mode(reshape(images,size(images,1),size(images,2), ...
            groupsize,size(images,3)/groupsize),3));
        case 'max'
          out = squeeze(max(reshape(images,size(images,1),size(images,2), ...
            groupsize,size(images,3)/groupsize),[],3));
        case 'min'
            out = squeeze(min(reshape(images,size(images,1),size(images,2), ...
                                groupsize,size(images,3)/groupsize),[],3));
        otherwise
            error('need a valid type');
    end
else
    switch type
        case 'sum'
          out = squeeze(sum(reshape(images(:,:,1:groupsize*floor(size(images,3)/groupsize)), ...
            size(images,1),size(images,2), groupsize,floor(size(images,3)/groupsize)),3));
          out = cat(3,out,sum(images(:,:,size(images,3)-rem(size(images,3),groupsize)+1:end),3));
        case 'prod'
          out = squeeze(prod(reshape(images(:,:,1:groupsize*floor(size(images,3)/groupsize)), ...
            size(images,1),size(images,2), groupsize,floor(size(images,3)/groupsize)),3));
          out = cat(3,out,prod(images(:,:,size(images,3)-rem(size(images,3),groupsize)+1:end),3));
        case 'mean'
          out = squeeze(mean(reshape(images(:,:,1:groupsize*floor(size(images,3)/groupsize)), ...
            size(images,1),size(images,2), groupsize,floor(size(images,3)/groupsize)),3));
          out = cat(3,out,mean(images(:,:,size(images,3)-rem(size(images,3),groupsize)+1:end),3));
        case 'median'
          out = squeeze(median(reshape(images(:,:,1:groupsize*floor(size(images,3)/groupsize)), ...
            size(images,1),size(images,2), groupsize,floor(size(images,3)/groupsize)),3));
          out = cat(3,out,median(images(:,:,size(images,3)-rem(size(images,3),groupsize)+1:end),3));
        case 'mode'
          out = squeeze(mode(reshape(images(:,:,1:groupsize*floor(size(images,3)/groupsize)), ...
            size(images,1),size(images,2), groupsize,floor(size(images,3)/groupsize)),3));
          out = cat(3,out,mode(images(:,:,size(images,3)-rem(size(images,3),groupsize)+1:end),3));
        case 'max'
          out = squeeze(max(reshape(images(:,:,1:groupsize*floor(size(images,3)/groupsize)), ...
            size(images,1),size(images,2), groupsize,floor(size(images,3)/groupsize)),[],3));
          out = cat(3,out,max(images(:,:,size(images,3)-rem(size(images,3),groupsize)+1:end),[],3));
        case 'min'
          out = squeeze(min(reshape(images(:,:,1:groupsize*floor(size(images,3)/groupsize)), ...
            size(images,1),size(images,2), groupsize,floor(size(images,3)/groupsize)),[],3));
          out = cat(3,out,min(images(:,:,size(images,3)-rem(size(images,3),groupsize)+1:end),[],3));
      otherwise
          error('need a valid type');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%