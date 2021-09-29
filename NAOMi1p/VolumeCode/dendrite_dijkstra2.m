function [distance,pathfrom] = dendrite_dijkstra2(M,dims,root)

% function [distance,pathfrom2] = dendrite_dijkstra2(M,dims,root)
% 
% Function to run Dijkstra's algorithm for growing out dendrites. Inputs to
% the function are
%    - M    - dims-size matrix that indicates blockages in the volume (i.e.
%             other dendrites, blood vessels, etc.)
%    - dims - 3x1 Size of the volume to grow dendrites in
%    - root - 3x1 Starting point of the path
%
% Outputs are:
%    - distance - Distance the path has traveled (dendrite length)
%    - pathfrom - Path through the volume that the dendrite takes
% 
% 2017 - Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    root2 = sub2ind(dims,root(1),root(2),root(3));                         % Get index set location of the root node from the subscript location provided
catch
    disp(root);
    error('root');
end

if(~isa(M,'single'))
    warning('M must be a single, casting as a single')
    M = single(M);                                                         % Make sure M is a single
end
e  = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];                        % Set the adjacent edges to loop through: R L U D T B (right/left/up/down/towards/back)
pe = e*[1;dims(1);dims(1)*dims(2)];                                        % Set distance based on dimensions of the volumes

[distance,pathfrom_tmp] = dendrite_dijkstra_cpp(M,int32(pe),int32(root2)); % Run the mexed dijkstra algorithm for speed

distance  = reshape(distance,dims);                                        % Reshape outputs to the size of the volume
pathfrom  = zeros(prod(dims),3);                                           % 
pidxs     = find(pathfrom_tmp>0);                                              % 
[pathfrom(pidxs,1),pathfrom(pidxs,2),pathfrom(pidxs,3)] = ...
                                        ind2sub(dims,pathfrom_tmp(pidxs)); % Change the path locations from an index set to subscript indexing
pathfrom = reshape(pathfrom,[dims 3]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%