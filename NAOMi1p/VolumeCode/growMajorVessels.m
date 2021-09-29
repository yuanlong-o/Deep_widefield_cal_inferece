function [nodes,nv] = growMajorVessels(nv,np,vp)

% [nodes,nv] = growMajorVessels(nv,np,vp)
%
% Function to sample vasculature locations and setup nodes and connections of
% vasculature (including connections to penetrating vessels)
%
% - nv           - struct for number of vessel parameters
%   .size        - size of volume [px]
%   .szum        - size of volume [um]
%   .nsource     - number of source vessels from volume surface
%   .ncapp       - number of capilliaries
%   .nvert       - number of penetrating vessels
%   .nnodes      - number of nodes
%   .nconn       - number of edges
%   .nvert_sum   - number of penetrating vessel connections to capilliaries
% - np           - node placement parameters
%   .maxit       - Maximum iteration to place nodes
%   .lensc       - average distance between vasculature branch points 
%   .varsc       - standard deviation of distances between vascualure 
%                  branch points  
%   .mindist     - minimum inter-node distance 
%   .varpos      - standard deviation of vasculature placement
%   .dirvar      - maximum branching angle
%   .branchp     - probability of branching surface vasculature
%   .vesrad      - radius of surface vasculature
% - vp           - struct for vasculature parameters
%   .depth_surf  - Depth into tissue of surface vasculature
%   .sepweight   - Set weight that guides how far, on average, the
%                  vasculature nodes are placed
%   .mindists    - set minimum distance between nodes
%   .distWeightScale - scaling factor for weight of node distance scaling 
%   .randWeightScale - scaling factor for weight of nodes (variability)
%   .maxcappdist - maximum capilliary distance
%
% - nodes        - struct array containing blood vessel node information
%   .num         - identifying number
%   .root        - root node, if it exists (0 is no root)
%   .conn        - connecting node(s)
%   .pos         - position of node
%   .type        - string, type of node (edge, surf, vert, sfvt, capp) 
%   .misc        - used for misc. parameters
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize a few points for vertical vessels, Initialize some points in surface for surface vessels
%% Create node struct type and grow initial (large) vessels from edges
nodes = gennode();                               
for i = 1:nv.nsource
  TMPIDX = (rand(1,2)>=[(nv.size(2)/(nv.size(1)+nv.size(2))) 0.5]);
  if      TMPIDX(1) &&  TMPIDX(2)
    TMPPOS = [ceil(nv.size(1)*rand) 1 vp.depth_surf];
  elseif  TMPIDX(1) && ~TMPIDX(2)
    TMPPOS = [ceil(nv.size(1)*rand) nv.size(2) vp.depth_surf];
  elseif ~TMPIDX(1) &&  TMPIDX(2)
    TMPPOS = [1 ceil(nv.size(2)*rand) vp.depth_surf];
  elseif ~TMPIDX(1) && ~TMPIDX(2)
    TMPPOS = [nv.size(1) ceil(nv.size(2)*rand) vp.depth_surf];
  end
  nodes(i) = gennode(i,0,[],TMPPOS,'edge',TMPIDX);
end

neur_surf = false(nv.size(1:2));
for i = 1:nv.nsource
  if      nodes(i).misc(1) &&  nodes(i).misc(2)
    randDir = 0.5*pi+randn*np.dirvar;
  elseif  nodes(i).misc(1) && ~nodes(i).misc(2)
    randDir = 1.5*pi+randn*np.dirvar;
  elseif ~nodes(i).misc(1) &&  nodes(i).misc(2)
    randDir = 0.0*pi+randn*np.dirvar;
  elseif ~nodes(i).misc(1) && ~nodes(i).misc(2)
    randDir = 1.0*pi+randn*np.dirvar;
  end
  [nodes,neur_surf] = branchGrowNodes(nodes,neur_surf,np,i,randDir);  
end
nv.nlinks = length(nodes)-nv.nsource;

for i = nv.nsource+1:(nv.nlinks+nv.nsource)
  nodes(i).pos = round([nodes(i).pos vp.depth_surf]);
end

%% Sample additional locations such that diving vessels can semi-uniformly cover the surface
neur_surf = imdilate(neur_surf,strel('disk',np.vesrad*2));

surfpos = pseudoRandSample2D(nv.size(1:2),nv.nsurf,vp.mindists(1),vp.sepweight,single(1-neur_surf)); % 
surfpos = cat(2,surfpos,vp.depth_surf*ones(nv.nsurf,1));   
surfpos = cat(1,reshape([nodes(1:nv.nlinks+nv.nsource).pos],[],nv.nlinks+nv.nsource)',surfpos);

% Forms a block matrix of connections:
% M = [A B]
%     [C D]
surfmat = pos2dists(surfpos);
surfmat(1:nv.nlinks+nv.nsource,1:nv.nlinks+nv.nsource) = inf;
surfmat(1:nv.nsource,1:nv.nsource) = 0;
for i = 1:nv.nlinks+nv.nsource
  if(nodes(i).root>0)
    surfmat(nodes(i).root,i) = norm(nodes(i).pos-nodes(nodes(i).root).pos);
    surfmat(i,nodes(i).root) = surfmat(nodes(i).root,i);
  end
end
surfmat(nv.nlinks+nv.nsource+1:end,1:nv.nlinks+nv.nsource) = inf;

TMPsurfmat = (surfmat.^vp.distWeightScale).*(1+...
                      vp.randWeightScale*randn(size(surfmat)));
[~,surfpath] = vessel_dijkstra(TMPsurfmat,1);
surfpath(1:nv.nsource)  = 1:nv.nsource;
for i = 1:nv.nsurf+nv.nsource
    if surfpath(i) == 1
        [~,surfpath(i)] = min(surfmat(i,1:nv.nsource));
    end
end

for i = (nv.nlinks+nv.nsource+(1:nv.nsurf))
  if(sum(surfpos(i,1)==[1 nv.size(1)]) || sum(surfpos(i,2)==[1 nv.size(2)]))
    nodes(i) = gennode(i,surfpath(i),surfpath(i),surfpos(i,:),'edge');    
  else
    nodes(i) = gennode(i,surfpath(i),surfpath(i),surfpos(i,:),'surf');
  end
end

nv.nnodes = nv.nlinks+nv.nsource+nv.nsurf;
for i = 1:nv.nnodes
  if(nodes(i).root>0)
    nodes(nodes(i).root).conn = union(nodes(nodes(i).root).conn,i);
  end
end

%% Prune some surface vasculature and choose diving vessels
se = strel('disk',round(vp.mindists(1)*2));
neur_vert = false(nv.size(1:2));
for i = 1:nv.nnodes
  if(strcmp(nodes(i).type,'surf') && (length(nodes(i).conn)==1))
    if((neur_vert(nodes(i).pos(1),nodes(i).pos(2))==0))
      nodes(i).type = 'sfvt';
      TMP = false(nv.size(1:2));
      TMP(nodes(i).pos(1),nodes(i).pos(2)) = 1;
      neur_vert = neur_vert+imdilate(TMP,se);
    else
      nodes = delnode(nodes,i);
    end
  end
end

nodes = nodes(1:nv.nnodes);
surfidx = find(strcmp({nodes.type},'surf'));
surfpos = reshape([nodes((strcmp({nodes.type},'surf'))).pos],3,[])';
surfpos = sub2ind(nv.size(1:2),surfpos(:,1),surfpos(:,2));
while(sum(strcmp({nodes.type},'sfvt'))<nv.nvert && sum(neur_vert(surfpos)==0)>0)
  TMPIDX = find(neur_vert(surfpos)==0);
  TMPIDX = surfidx(TMPIDX(ceil(rand*length(TMPIDX))));
  nodes(TMPIDX).type = 'sfvt';
  TMP = false(nv.size(1:2));
  TMP(nodes(TMPIDX).pos(1),nodes(TMPIDX).pos(2)) = 1;
  neur_vert = neur_vert+imdilate(TMP,se);
end

%% Grow diving vessels to bottom of volume
vertidx = find(strcmp({nodes.type},'sfvt'));
TMPIDX = nv.nnodes;
for i = 1:length(vertidx)
  curr_node = vertidx(i);
  while(nodes(curr_node).pos(3)<nv.size(3))
    TMPIDX = TMPIDX + 1;
    node_pos = nodes(curr_node).pos+ceil([randn(1,2) 1].*[np.varpos np.varpos  max(np.varsc*randn+np.lensc,np.mindist)]);
    node_pos = min(max(node_pos,[1 1 1]),nv.size);
    nodes(TMPIDX) = gennode(TMPIDX,curr_node,curr_node,node_pos,'vert');
    nodes(nodes(TMPIDX).root).conn = union(nodes(nodes(TMPIDX).root).conn,TMPIDX);
    curr_node = TMPIDX;
  end
end
nv.nvert = sum(strcmp({nodes.type},'sfvt'));
nv.nvertconn = curr_node-nv.nnodes;
nv.nnodes = curr_node;

%% End nodes are initialized to a size
ends = find(cellfun(@length,{nodes.conn})==1);
for i = 1:length(ends)
  nodes(ends(i)).misc = vp.vesSize(3)+gamrnd(3,(vp.vesSize(2)-vp.vesSize(3))/3); % gamma distribution with shape parameter 3 to set the vessel distribution
end