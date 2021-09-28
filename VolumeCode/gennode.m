function node = gennode(num,root,conn,pos,type,misc)

% node = gennode(num,root,conn,pos,type,misc)
%
% Function to create a node struct (branch point of blood vessels).  Nodes 
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if(nargin<1);num  = [];end
if(nargin<2);root = [];end
if(nargin<3);conn = [];end
if(nargin<4);pos  = [];end
if(nargin<5);type = [];end
if(nargin<6);misc = [];end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set node values

node.num  = num;
node.root = root;
node.conn = conn;
node.pos  = pos;
node.type = type;
node.misc = misc;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%