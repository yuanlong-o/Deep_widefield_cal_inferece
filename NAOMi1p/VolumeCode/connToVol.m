function [neur_ves,conn] = connToVol(nodes,conn,nv,idxs,neur_ves)

% [neur_ves,conn] = connToVol(nodes,conn,nv,idxs,neur_ves)
% 
% Function to add vasculature connectivity to a volume. Inputs to this function
% are:
% 
% - nodes       - struct array containing blood vessel node information
%     .num      - identifying number
%     .root     - root node, if it exists (0 is no root)
%     .conn     - connecting node(s)
%     .pos      - position of node
%     .type     - string, type of node (edge, surf, vert, sfvt (surf/vert), capp) 
%     .misc     - used for misc. parameters
% - conn        - struct array containing blood vessel edge information
%     .start    - start node
%     .ends     - end node
%     .weight   - node weight
%     .locs     - position of locations
%     .misc     - used for misc. parameters
% - nv.size     - volume size (pixels)
% - idxs        - indexes of vessels to generate
% - neur_ves    - simulated blood vessel volume
% 
% The ouputs of this function are
% - neur_ves    - Updated simulated blood vessel volume
% - conn        - Updated struct array containing blood vessel edge information
%     .start    - start node
%     .ends     - end node
%     .weight   - node weight
%     .locs     - position of locations
% 
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<5)
  neur_ves = false(nv.size);
end
if(nargin<4)
  idxs = 1:length(conn);
end

for j = 1:length(idxs)
  i = idxs(j);
  TMPL = setxor(nodes(conn(i).start).conn,conn(i).ends);
  if(~isempty(TMPL))
    TMPL = TMPL(randi(length(TMPL)));
  end
  TMPU = setxor(nodes(conn(i).ends).conn,conn(i).start);
  if(~isempty(TMPU))
    TMPU = TMPU(randi(length(TMPU)));
  end
  numsamp = 2*norm(nodes(conn(i).start).pos'-nodes(conn(i).ends).pos');
  if(~isempty(TMPL) && ~isempty(TMPU))
    spts = cscvn([nodes(TMPL).pos' nodes(conn(i).start).pos' nodes(conn(i).ends).pos' nodes(TMPU).pos']);
    ves_loc = ceil(fnval(spts,linspace(spts.breaks(2),spts.breaks(3),numsamp))');
  elseif(~isempty(TMPL))
    spts = cscvn([nodes(TMPL).pos' nodes(conn(i).start).pos' nodes(conn(i).ends).pos']);
    ves_loc = ceil(fnval(spts,linspace(spts.breaks(2),spts.breaks(3),numsamp))');
  elseif(~isempty(TMPU))
    spts = cscvn([nodes(conn(i).start).pos' nodes(conn(i).ends).pos' nodes(TMPU).pos']);
    ves_loc = ceil(fnval(spts,linspace(spts.breaks(1),spts.breaks(2),numsamp))');
  else
    spts = cscvn([nodes(conn(i).start).pos' nodes(conn(i).ends).pos']);
    ves_loc = ceil(fnval(spts,linspace(spts.breaks(1),spts.breaks(2),numsamp))');
  end
  ves_loc = bsxfun(@max,ves_loc,[1 1 1]);
  ves_loc = bsxfun(@min,ves_loc,nv.size);
  ves_loc = unique(ves_loc,'rows');
  conn(i).locs = ves_loc;
  
  min_idx = max(min(ves_loc)-ceil(conn(i).weight),[1 1 1]);
  max_idx = min(max(ves_loc)+ceil(conn(i).weight),nv.size);
  ves_loc = bsxfun(@minus,ves_loc,min_idx)+1;
  
  TMP = false(max_idx-min_idx+1);
  TMP(sub2ind(max_idx-min_idx+1,ves_loc(:,1),ves_loc(:,2),ves_loc(:,3))) = 1;
  [x,y,z] = ndgrid(-ceil(conn(i).weight):ceil(conn(i).weight));
  se      = strel(sqrt(x.^2 + y.^2 + z.^2) <=conn(i).weight);
  TMP     = imdilate(TMP,se);
  neur_ves(min_idx(1):max_idx(1),min_idx(2):max_idx(2),min_idx(3):max_idx(3)) = ...
    neur_ves(min_idx(1):max_idx(1),min_idx(2):max_idx(2),min_idx(3):max_idx(3))+TMP;
end