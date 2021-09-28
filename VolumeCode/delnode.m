function nodes = delnode(nodes,num)

% nodes = delnode(nodes,num)
%
% Function to delete a single node (branch point of blood vessels). Nodes 
% are arranged as a structure array with fields:
%
%   node.num      - identifying number
%   node.root     - root node, if it exists (0 is no root)
%   node.conn     - connecting node(s)
%   node.pos      - position of node
%   node.type     - string, type of node (edge, surf, vert, sfvt (surf/vert), capp)
%   node.misc     - used for misc. parameters
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = nodes(num).conn                                                    % Iterate through all the nodes
    nodes(i).conn(nodes(i).conn==num) = [];
    if(nodes(i).root==num)
        nodes(i).root = [];
    end
end
nodes(num) = gennode;                                                

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%