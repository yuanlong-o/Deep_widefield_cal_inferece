function conn = genconn(start,ends,weight,locs,misc)

% conn = genconn(start,ends,weight,locs,misc)
%
% Function to generate a connection structure (piece of blood vessel; 
% edge between nodes) Connections are a structure array with fields: 
%
%   conn.start      - start node
%   conn.ends       - end node
%   conn.weight     - node weight
%   conn.locs       - position of locations
%   conn.misc       - used for misc. parameters
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<1);start  = [];end
if(nargin<2);ends   = [];end
if(nargin<3);weight = [];end
if(nargin<4);locs   = [];end
if(nargin<5);misc   = [];end

conn.start  = start;
conn.ends   = ends;
conn.weight = weight;
conn.locs   = locs;
conn.misc   = misc;
