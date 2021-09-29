function [path,pathM] = getDendritePath2(M,node,root)

% [path,pathM] = getDendritePath2(M,node,root)
%
% Function to retrieve the path of a dendrite from the full paths matrix until
% it hits the root node. Inputs are:
%   - M    - pathfrom matrix (matrix containing obsticals to avoid)
%   - node - End node location in the volume
%   - root - Starting location for the dendrite path
%
% Outputs are:
%   - path  - retrieved path (list of positions)
%   - pathM - full path matrix (binary matrix of path)
% 
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = node;                                                               % Initialize the path to the starting location node

if length(node)==2
    while(~isequal(node,root))
        node = squeeze(M(node(1),node(2),:))';
        path = [path;node];
    end
  
    if nargout > 1
        pathM = false(size(M,1),size(M,2));
        pathM(sub2ind([size(M,1),size(M,2)],path(:,1),path(:,2)))=1;
    end
elseif length(node)==3
    path      = zeros(sum(size(M)),3);                                     % Initialize path to zeros
    path(1,:) = node;                                                      % Initialize path start to the starting node
    i         = 1;
    while(~isequal(node,root))
        try
            i         = i+1;
            node      = reshape(M(node(1),node(2),node(3),:),1,[]);
            path(i,:) = node;
        catch
            path = [];
            break;
        end
    end
  if ~isempty(path)
    path = path(1:i,:);
  end
  if nargout > 1
    pathM = false(size(M,1),size(M,2),size(M,3));
    if(~isempty(path))
      pathM(sub2ind([size(M,1),size(M,2),size(M,3)],path(:,1),path(:,2),path(:,3)))=1;
    end
  end
else
  error('number of dimension of node is not 2 or 3')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%