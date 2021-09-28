function [nodes,conn,nv] = growCapillaries(nodes,conn,neur_ves,nv,vp,vres)

% [nodes,conn,nv] = growCapillaries(nodes,conn,neur_ves,nv,vp,vres)
%
% Function to sample capillary locations and setup nodes and connections of
% capillaries (including connections to penetrating vessels)
%
% - nodes        - struct array containing blood vessel node information
%   .num         - identifying number
%   .root        - root node, if it exists (0 is no root)
%   .conn        - connecting node(s)
%   .pos         - position of node
%   .type        - string, type of node (edge, surf, vert, sfvt, capp) 
%   .misc        - used for misc. parameters
% - conn         - struct array containing blood vessel edge information
%   .start       - start node
%   .ends        - end node
%   .weight      - node weight
%   .locs        - position of locations
%   .misc        - used for misc. parameters
% - nv           - struct for number vessel parameters
%   .szum        - size of volume [um]
%   .ncapp       - number of capilliaries
%   .nvert       - number of penetrating vessels
%   .nnodes      - number of nodes
%   .nconn       - number of edges
%   .nvert_sum   - number of penetrating vessel connections to capilliaries
% - vp           - struct for vasculature parameters
%   .maxcappdist - maximum capilliary distance
%   .mindists    - set minimum distance between nodes
%   .vesFreq     - blood vessel freqency in px
%   .vesSize     - vessel radius (surface, axial, capillaries) in px
%   .sepweight   - Set weight that guides how far, on average, the
%                  vasculature nodes are placed 
%   .ves_shift   - 3-vector of the amount of wobble allowed for blood
%                  vessels
%   .distsc      - How strongly local capillary connections are. Higher 
%                  numbers indicate more local connections
%   .vesNumScale - blood vessel number random scaling factor
% - vres         - volume resolution [1/um]
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample capillary locations and setup nodes and connections

% Initialize a pseudouniform sampling of 3D positions for capillaries
dilrad = ceil(vp.mindists(1)/vres);
[x,y,z] = ndgrid(-dilrad:dilrad);
se = exp(-2*(x.^2 + y.^2 + z.^2)/dilrad^2);
TMPvol = zeros(nv.szum,'single');
for i = 1:length(conn)
  TMPpos = ceil(conn(i).locs/vres);
  TMPpos = TMPpos(1:max(1,floor(dilrad/3)):end,:);
  for j = 1:size(TMPpos,1)
    TMPL = TMPpos(j,:)-dilrad;
    TMPL = (TMPL<1).*(1-TMPL);
    TMPU = TMPpos(j,:)+dilrad;
    TMPU = (TMPU>nv.szum).*(TMPU-nv.szum);
    TMP = se(TMPL(1)+1:end-TMPU(1),TMPL(2)+1:end-TMPU(2),TMPL(3)+1:end-TMPU(3));
    TMPL = TMPpos(j,:)-dilrad+TMPL;
    TMPU = TMPpos(j,:)+dilrad-TMPU;
    TMPvol(TMPL(1):TMPU(1),TMPL(2):TMPU(2),TMPL(3):TMPU(3)) = max(TMPvol(TMPL(1):TMPU(1),TMPL(2):TMPU(2),TMPL(3):TMPU(3)),TMP);
  end
end

capppos  = pseudoRandSample3D(nv.size/vres,nv.ncapp,vp.mindists(3)/vres,vp.sepweight,single(1-TMPvol));
capppos  = capppos*vres+1-randi(vres,size(capppos));

% Setup connections from diving vessels to capillaries
nv_vert_conn = randi(ceil(nv.szum(3)/vp.vesFreq(3)),[1 nv.nvert]);
nv.nvert_sum = sum(nv_vert_conn);
nodeIdx = nv.nnodes;
connIdx = nv.nconn;
vertidxs = find(strcmp({nodes.type},'sfvt'));
if(length(vp.vesSize)<4)
  vp.vesSize(4) = vp.vesSize(3)*0;
%   vp.vesSize(4) = vp.vesSize(3)*0.25;
end
for i = 1:length(vertidxs)
  vesidx = vertidxs(i);
  flag = 1;
  while(flag)
    TMPidx = nodes(vesidx(end)).conn;
    TMPidx = TMPidx(strcmp({nodes(TMPidx).type},'vert'));
    TMPidx = intersect(setxor(vesidx,TMPidx),TMPidx);
    if(isempty(TMPidx))
      flag = 0;
    else
      vesidx = [vesidx;TMPidx];
    end
  end
  vesidx = vesidx(2:end);
  for j = 1:nv_vert_conn(i)
    TMPidx = vesidx(randi(length(vesidx)));
    nodeIdx = nodeIdx+1;
    connIdx = connIdx+1;
    [~,TMP] = min(sum(bsxfun(@minus,capppos,nodes(TMPidx).pos).^2,2));
    nodes(nodeIdx) = gennode(nodeIdx,TMPidx,TMPidx,capppos(TMP,:),'capp');
    nodes(TMPidx).conn = union(nodes(TMPidx).conn,nodeIdx);
    conn(connIdx) = genconn(nodeIdx,TMPidx,max(1,normrnd(vp.vesSize(3),vp.vesSize(4))),[],'vtcp');
    capppos(TMP,:) = nan(1,3);
  end
end
nv.nnodes = nodeIdx;
nv.nconn = connIdx;

nodeIdx = nv.nnodes;
for i = 1:size(capppos,1)
  if(~isnan(capppos(i,:)))
    nodeIdx = nodeIdx+1;
    nodes(nodeIdx) = gennode(nodeIdx,[],[],capppos(i,:),'capp');
  end
end
nv.nnodes = nodeIdx;

% Setup connections from capillaries to each other
vertconnidxs = find(~(cellfun(@isempty,{nodes.root})).*strcmp({nodes.type},'capp'));
cappconnidxs = find((cellfun(@isempty,{nodes.root})).*strcmp({nodes.type},'capp'));
connidxs = [vertconnidxs cappconnidxs];
capppos = reshape([nodes(vertconnidxs).pos nodes(cappconnidxs).pos],3,[])';

cappmat  = pos2dists(capppos);                                             % scale by vessels and blocking in way
cappmat(logical(eye(nv.ncapp))) = inf;
cappmat(1:nv.nvert_sum,1:nv.nvert_sum) = inf;


cappconnmat = zeros(nv.ncapp);
[~,mincapp] = min(cappmat);

cappconnmat(sub2ind([nv.ncapp nv.ncapp],1:nv.ncapp,mincapp)) = 1;
cappconnmat(sub2ind([nv.ncapp nv.ncapp],mincapp,1:nv.ncapp)) = 1;
cappmat(sub2ind([nv.ncapp nv.ncapp],1:nv.ncapp,mincapp)) = inf;
cappmat(sub2ind([nv.ncapp nv.ncapp],mincapp,1:nv.ncapp)) = inf;

cappmat(cappmat>vp.maxcappdist) = inf;
for i=1:nv.ncapp
  if(sum(cappconnmat(i,:))>=3)
    cappmat(i,:) = inf;
    cappmat(:,i) = inf;
  end
  cappmat(i,logical(cappconnmat(i,:))) = inf;
  cappmat(logical(cappconnmat(i,:)),i) = inf;
end

for i = (nv.nvert_sum+1):nv.ncapp
  for j = (i+1):nv.ncapp
    if(cappmat(i,j)<inf)
      xpix = ceil([linspace(capppos(i,1),capppos(j,1),2*cappmat(i,j))' ...
        linspace(capppos(i,2),capppos(j,2),2*cappmat(i,j))' ...
        linspace(capppos(i,3),capppos(j,3),2*cappmat(i,j))']);
      if(sum(neur_ves(sub2ind(nv.size,xpix(:,1),xpix(:,2),xpix(:,3)))))
        cappmat(i,j) = inf;
        cappmat(j,i) = inf;
      end
    end
  end
end
lflag = 1;
while(lflag)
  cappsum = sum(cappconnmat);
  if(min(cappsum((nv.nvert_sum+1):end))>1)
    lflag = 0;
  else
    [~,idxs] = find(cappsum==1);
    if(min(min(cappmat(:,idxs)))==inf)
      lflag = 0;
    end
    rndidx = idxs(ceil(length(idxs)*rand));
    if(min(cappmat(rndidx,:))<inf)
      capdistinv = 1./(cappmat(rndidx,:).^vp.distsc);                      % distance needs to avoid certain obstacles scale in cappmat
      capcdf = [0 cumsum(capdistinv)/sum(capdistinv)];
      lnkidx = find(diff(capcdf>rand));
      
      cappconnmat(rndidx,lnkidx) = 1;
      cappconnmat(lnkidx,rndidx) = 1;
      
      cappmat(logical(cappconnmat(rndidx,:)),lnkidx) = inf;
      cappmat(logical(cappconnmat(lnkidx,:)),rndidx) = inf;
      cappmat(lnkidx,logical(cappconnmat(rndidx,:))) = inf;
      cappmat(rndidx,logical(cappconnmat(lnkidx,:))) = inf;
      cappmat(rndidx,lnkidx) = inf;
      cappmat(lnkidx,rndidx) = inf;
      
      if(sum(cappconnmat(lnkidx,:))>=3)
        cappmat(lnkidx,:) = inf;
        cappmat(:,lnkidx) = inf;
      end
    end
  end
end
[connS,connF] = ind2sub([nv.ncapp nv.ncapp],find(triu(cappconnmat)));

% Setup connection structure from capillary connections
connIdx = nv.nconn;
connMat = sparse(length(nodes),length(nodes));
for i = 1:length(connS)
  nodes(connidxs(connS(i))).conn = union(nodes(connidxs(connS(i))).conn,connidxs(connF(i)));
  nodes(connidxs(connF(i))).conn = union(nodes(connidxs(connF(i))).conn,connidxs(connS(i)));
  connIdx = connIdx+1;
  conn(connIdx) = genconn(connidxs(connS(i)),connidxs(connF(i)),NaN,[],'capp');
  connMat(connidxs(connS(i)),connidxs(connF(i))) = connIdx;
end
nv.nconn = connIdx;
nodesToConnect = [conn(strcmp({conn.misc},'vtcp')).start];
toConnect = [];
for i = 1:length(nodesToConnect)
  TMPconn = nodes(nodesToConnect(i)).conn;
  for j = 1:length(TMPconn)
    if(connMat(nodesToConnect(i),TMPconn(j)))
      toConnect = [toConnect; connMat(nodesToConnect(i),TMPconn(j))];
    end
  end
end
toConnect = full(toConnect);
connToConnect = find(strcmp({conn.misc},'vtcp'));
for i = 1:length(connToConnect)
  connMat(conn(connToConnect(i)).ends,conn(connToConnect(i)).start) = connToConnect(i);
end
connMat = connMat+connMat';

while(~isempty(toConnect))
  currConn = toConnect(1);
  if(isnan(conn(currConn).weight))
    connStart = conn(currConn).start;
    connEnd = conn(currConn).ends;
    
    startConns = nodes(conn(currConn).start).conn;
    endConns = nodes(conn(currConn).ends).conn;
    startConns = full(connMat(connStart,startConns));
    endConns = full(connMat(connEnd,endConns));
    startConns(startConns==currConn) = [];
    startConns(startConns==0) = [];
    endConns(endConns==currConn) = [];
    endConns(endConns==0) = [];
    
    startWeights = [conn(startConns).weight];
    endWeights = [conn(endConns).weight];
    endflag = 0;
    startflag = 0;
    if(sum(isnan(startWeights)))
      weight1 = NaN;
    else
      if(length(startWeights)==1)
        startflag = 1;
        weight1 = startWeights;
      else
        TMP1 = (max(startWeights)^2)-(min(startWeights)^2);
        TMP2 = (max(startWeights)^2)+(min(startWeights)^2);
        weight1 = sqrt((rand*(TMP2-TMP1))+TMP1);
      end
    end
    
    if(sum(isnan(endWeights)))
      weight2 = NaN;
    else
      if(length(endWeights)==1)
        endflag = 1;
        weight2 = endWeights;
      else
        TMP1 = (max(endWeights)^2)-(min(endWeights)^2);
        TMP2 = (max(endWeights)^2)+(min(endWeights)^2);
        weight2 = sqrt((rand*(TMP2-TMP1))+TMP1);
      end
    end
    if(isnan(weight1))
      if(isnan(weight2))
        connweight = max(1,normrnd(vp.vesSize(3),vp.vesSize(4)));
      else
        connweight = weight2;
      end
    else
      if(isnan(weight2))
        connweight = weight1;
      else
        if(startflag)
          connweight = weight1;
        elseif(endflag)
          connweight = weight2;
        else
          connweight = (weight1+weight2)/2;
        end
      end
    end
    conn(currConn).weight = connweight;
    toConnect = [toConnect;endConns(isnan(endWeights))'; startConns(isnan(startWeights))'];
  end
  toConnect = toConnect(2:end);
end

for i = 1:nv.nconn
  if(isempty(conn(i).weight)||isnan(conn(i).weight))
    conn(i).weight = max(1,normrnd(vp.vesSize(3),vp.vesSize(4)));
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
