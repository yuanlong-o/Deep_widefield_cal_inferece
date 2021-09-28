function conn = nodesToConn(nodes)

% conn = nodesToConn(nodes)
%
% Function that generates the connection structure from the node structure
% for blood vessels. The input is:
%
% - nodes       - struct array containing blood vessel node information
%     .num      - identifying number
%     .root     - root node, if it exists (0 is no root)
%     .conn     - connecting node(s)
%     .pos      - position of node
%     .type     - string, type of node (edge, surf, vert, sfvt (surf/vert), capp) 
%     .misc     - used for misc. parameters
%
% The output is:
% - conn        - struct array containing blood vessel edge information
%     .start    - start node
%     .ends     - end node
%     .weight   - node weight
%     .locs     - position of locations
%     .misc     - used for misc. parameters
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ends = find(cellfun(@length,{nodes.conn})==1);
connmat = spalloc(length(nodes),length(nodes),ceil(length([nodes.conn])/2));

for i = 1:length(ends)
  curr_node = ends(i);
  while(nodes(curr_node).root>0)
    connmat(nodes(curr_node).num,nodes(curr_node).root) = ...
      sqrt(connmat(nodes(curr_node).num,nodes(curr_node).root)^2+nodes(ends(i)).misc^2);
    curr_node = nodes(curr_node).root;
  end
end
[x,y,w] = find(connmat);
conn = genconn();
for i = 1:length(x)
  conn(i) = genconn(x(i),y(i),w(i));
end
