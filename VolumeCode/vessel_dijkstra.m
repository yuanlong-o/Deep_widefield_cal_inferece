function [distance,pathfrom] = vessel_dijkstra(distMat,proot)

% [distance,pathfrom] = vessel_dijkstra(distMat,proot)
%
% Function to apply Dijkstra's algorithm to growing out vasculature. 
% The inputs to this function are:
%
%  - distMat  - distance matrix between different nodes
%  - proot    - starting root node
%
% The ouputs of this function are:
%  - distance - minimum path distance from all nodes to root node
%  - pathfrom - path from root to each node
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dims             = size(distMat,1);
tovisit          = true(dims,1);
unvisited        = zeros(dims,1);
unvisited(proot) = 1;
distance         = inf(dims,1);
distance(proot)  = 0;
pathfrom         = nan(dims,1);
cn               = proot;                                                  % current node being examined

while(nnz(unvisited))
    tovisit(cn)   = 0;
    nextidx       = find(unvisited);
    [~,idx]       = min(nonzeros(unvisited));
    cn            = nextidx(idx);
    unvisited(cn) = 0;
    for nn = 1:dims
        if tovisit(nn)
            ndist = distance(cn)+distMat(cn,nn);
            if ndist < distance(nn)
                unvisited(nn) = max(eps,ndist);
                distance(nn)  = ndist;
                pathfrom(nn)  = cn;
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%