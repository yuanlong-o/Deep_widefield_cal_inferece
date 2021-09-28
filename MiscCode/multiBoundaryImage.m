function I = multiBoundaryImage(X,varargin)

% I = multiBoundaryImage(X)
% 
% Function to create a series of bounding curves for multiple images in a
% 3D stack. 
% 
% Input:  X - a 3D array of images, each slice in the 3rd dimension will be
%             processed independently.
%         t - a threshold parameter for what is "inside" the object.
%             (defaults to 0.01)
%
% Output: I - A cell array the same length as the number of slices in X.
%             Each index contains the bounding curve information related to
%             the corresponding image in X.
% 
% 2018 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 1
    t = varargin{1};
else
    t = 0.01;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get all the boundaries

I = cell(numel(X,3),1);

for kk = 1:size(X,3)
    TMP       = X(:,:,kk);
    TMP       = TMP > t*max(TMP(:));
    if sum(TMP(:)) > 0
        [row,col] = find(TMP == 1,1,'first');
        I{kk}     = bwtraceboundary(TMP,[row, col],'N');
    else
        I{kk} = [];
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%