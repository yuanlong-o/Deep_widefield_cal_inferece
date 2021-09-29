function [neur_num,dendnumAD,dend_params,gp_soma] = growNeuronDendrites(vol_params, dend_params, neur_soma, neur_ves, neur_locs, gp_nuc, gp_soma, rotAng)

% [neur_num,cellVolumeAD,dend_params] = growNeuronDendrites(vol_params, dend_params, neur_soma, neur_ves, neur_locs, gp_nuc, gp_soma, rotAng)
% 
% Function to grow dendrites in a neural volume. The inputs to this
% function are:
%   - vol_params  - Struct containing parameters for the volume generation
%       .vol_sz   - 3-element vector with the size (in um) of the volume to
%                   generate (default = 100x100x30um)
%       .min_dist - Minimum distance between neurons (default = 15um)
%       .N_neur   - Number of neurons to generate (default = 50)
%       .vres     - resolution to simulate volume at (default = 2
%                   samples/um)
%       .N_den    - Width of apical dendrites (default = 10)
%       .N_bg     - Number of background/neuropil components to simulate
%                   (default = 50)
%       .vol_depth- Depth of the volume under the brain surface (default =
%                   100um)
%       .verbose  - Level of verbosity in the output during the volume
%                   generation. Can be 0,1,2. 0 = no text updates, 1 = some
%                   some text outputs. 2 = detailed text outputs (default =
%                   1)
%   - dend_params - Struct containing parameters for dendrite simulation
%       .dtParams        - dendritic tree number,radius in um of branches
%                          (uniform across circle),radius in um of branches
%                          (uniform in z), width scaling, variation in
%                          dendritic tree number (default = [40 150 50 1 10])
%       .atParams        - Apical dendrite number,radius in um of branches
%                          (uniform across circle),radius in um of branches
%                          (uniform in z),offset from center in um (default
%                          = [1 5 2 2 4]) 
%       .dweight         - Weight for path planning randomness in the
%                          dendrites (default = 10) 
%       .bweight         - Weight for obstruction (default = 50)
%       .thicknessScale  - Scaling for dendrite thickness (int). Should be
%                          1 for 1um sampling,(4 for 0.5um sampling)
%                          (default = 0.75) 
%       .dims            - dims set at 10 um per space (default = [20 20 20])
%       .dimsSS          - dims subsampling factor (10 samples per dim
%                          grid) (default = [10 10 10]) 
%   - neur_soma - Array where the locations (actual points) of the k^th
%                 neural soma are represented by the value 'k'
%   - neur_ves  - Array deliniating where the vasculature in the volume is
%                 located
%   - neur_locs - Nx3 location of neurons in the volume
%   - gp_nuc    - Weights of the nucleus
%   - rotAng    - Nx3 Vector with the rotation angle of the cell (Rx,Ry,Rz)
%
%  The outputs for this function are:
%   - neur_num     - An array where the k^th neuron's locations in the
%                    volume (both the soma and dendrites) are deliniated by
%                    the value 'k'
%   - cellVolumeAD - Volume array containing the number for each apical
%                    dendrite at each location it occupies
%   - dend_params  - The updated struct of dendrite parameters (in case the
%                    default values were added)
% 
% 2017 - Alex Song and Adam Charles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing

if(nargin<8)
  rotAng = [];
end
vol_params  = check_vol_params(vol_params);                                % Check volume parameters
dend_params = check_dend_params(dend_params);                              % Check dendrite parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% General Parameter Setup



if vol_params.verbose == 1
    fprintf('Growing out dendrites')
elseif vol_params.verbose >1
    fprintf('Growing out dendrites...\n')
end

dtParams       = dend_params.dtParams;                                     %
atParams       = dend_params.atParams;                                     %
dweight        = dend_params.dweight;                                      %
bweight        = dend_params.bweight;                                      %
thicknessScale = dend_params.thicknessScale;                               %
dims           = dend_params.dims;                                         %
dimsSS         = dend_params.dimsSS;                                       %
rallexp        = dend_params.rallexp;                                      %
vres           = vol_params.vres;                                          % Extract the volume resolution variable (points per micron)
N_neur         = vol_params.N_neur;                                        % Extract the [total] number of neurons
vol_sz         = vol_params.vol_sz;                                        % Extract volume size
dims           = min(dims, [vol_sz(1), vol_sz(2), vol_sz(3)]./dimsSS);     % dims set at dimsSS um per space
fulldims       = vol_sz*vres;                                              % Get dimension of the volume to work with (sampling times size)
dims           = dims*vres;                                                % Calculate dimensions of the sampled array
dtParams(2:3)  = dtParams(2:3)*vres;
atParams(2:4)  = atParams(2:4)*vres;
thicknessScale = thicknessScale*vres*vres;                                 % Calculate thicknes of dendrites

vol_depth      = vol_params.vol_depth*vol_params.vres;                     % Calculate the depth of the volume
cellVolume     = single(neur_soma)+(vol_params.N_den+vol_params.N_neur+...
                vol_params.N_bg+1)*neur_ves(:,:,(vol_depth+1):vol_depth...
                                             +vol_sz(3)*vol_params.vres);% Create a cell volume to store all the values of somas and dendrites in
for kk = 1:N_neur
    cellVolume(gp_nuc{kk,1}) = kk;
end

if nargin < 7
    gp_soma = cell(N_neur,1);
    for kk = 1:N_neur
        gp_soma{kk,1} = int32(vec(find(neur_soma==kk)));
    end
end

neur_num = uint16(cellVolume);
  
cellVolumeIdx  = zeros(size(neur_soma), 'single');
cellVolumeVal  = zeros(size(neur_soma), 'single');                         % Initialize the volume of apical dendrites
cellVolumeAD   = false(size(neur_soma));                                   % Initialize the volume of apical dendrites
allroots       = ceil(max(vres*neur_locs,1e-4));                           % Initialize an array for the roots of the dendrites
fdims          = dims.*dimsSS;                                             % 
fdims          = min(fdims, fulldims);                                     % 
ML             = inf([fdims 6],'single');                                           % Make sure that fdims is always greater than 6
if (~isfield(dend_params,'dendVar'))||isempty(dend_params.dendVar)
  dendVar = 0.25;
else
  dendVar = dend_params.dendVar;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grow Dendrites

smallZ = fulldims(3)<=fdims(3);
for j = 1:N_neur
    if vol_params.verbose > 1
        tic
    end
    [aproot(1),aproot(2),aproot(3)] = ind2sub(fulldims,min(gp_soma{j}));
    
    numdt = max(1,dtParams(1)+round(dtParams(5)*randn));

    borderflag = 0;
    try
      if(smallZ)
        rootL = [(fdims(1)/2+1) (fdims(2)/2+1) allroots(j,3)];               % Center origin point (in um)
        obstruction = cellVolume(allroots(j,1)+(1:fdims(1))-fdims(1)/2-1,...
                                 allroots(j,2)+(1:fdims(2))-fdims(2)/2-1,:);
      else
        rootL = [(fdims(1)/2+1) (fdims(2)/2+1) (fdims(3)/2+1)];              % Center origin point (in um)
        obstruction = cellVolume(allroots(j,1)+(1:fdims(1))-fdims(1)/2-1,...
                                 allroots(j,2)+(1:fdims(2))-fdims(2)/2-1,...
                                 allroots(j,3)+(1:fdims(3))-fdims(3)/2-1);
      end
    catch
      if(smallZ)
        rootL = [(fdims(1)/2+1) (fdims(2)/2+1) allroots(j,3)];               % Center origin point (in um)
        obstruction = zeros(fdims,'single');                               % Set up an array of obstructions (somas, vasculature, other dendrites)
        Xlims = [(allroots(j,1)-fdims(1)/2) (allroots(j,1)+fdims(1)/2-1)];
        Ylims = [(allroots(j,2)-fdims(2)/2) (allroots(j,2)+fdims(2)/2-1)];
        obstruction( ...
         max([1 -Xlims(1)+2]):(fdims(1)+min([0 -Xlims(2)+fulldims(1)])),...
         max([1 -Ylims(1)+2]):(fdims(2)+min([0 -Ylims(2)+fulldims(2)])),:) = ...
         cellVolume( ...
         max([1 Xlims(1)]):min([Xlims(2) fulldims(1)]), ...
         max([1 Ylims(1)]):min([Ylims(2) fulldims(2)]),:);
         borderflag = 1;
      else
        rootL = [(fdims(1)/2+1) (fdims(2)/2+1) (fdims(3)/2+1)];              % Center origin point (in um)
        obstruction = zeros(fdims,'single');                               % Set up an array of obstructions (somas, vasculature, other dendrites)
        Xlims = [(allroots(j,1)-fdims(1)/2) (allroots(j,1)+fdims(1)/2-1)];
        Ylims = [(allroots(j,2)-fdims(2)/2) (allroots(j,2)+fdims(2)/2-1)];
        Zlims = [(allroots(j,3)-fdims(3)/2) (allroots(j,3)+fdims(3)/2-1)];
        obstruction( ...
         max([1 -Xlims(1)+2]):(fdims(1)+min([0 -Xlims(2)+fulldims(1)])),...
         max([1 -Ylims(1)+2]):(fdims(2)+min([0 -Ylims(2)+fulldims(2)])),...
         max([1 -Zlims(1)+2]):(fdims(3)+min([0 -Zlims(2)+fulldims(3)]))) = ...
         cellVolume( ...
         max([1 Xlims(1)]):min([Xlims(2) fulldims(1)]), ...
         max([1 Ylims(1)]):min([Ylims(2) fulldims(2)]), ...
         max([1 Zlims(1)]):min([Zlims(2) fulldims(3)]));% Figure out where all the obstructions are in the volume
        borderflag = 1;
      end
    end
    cellBody = find(obstruction==j);
    obstruction(cellBody) = 0;                                             % Set obstructions

    root = ceil(rootL./dimsSS);
    root = min(root,dims);
    M = 1+dweight*rand(dims(1),dims(2),dims(3),6);                         % Lattice of values, in order of R,L,U,D,T,B

    aprootS = root+round((aproot-allroots(j,:))/dimsSS(3));
    aprootS = max(aprootS,[1 1 1]);
    aprootS = min(aprootS,dims);    
    if(aprootS(1)>root(1))
      M(root(1):aprootS(1),root(2),root(3),1) = 0;
    elseif(aprootS(1)<root(1))
      M(aprootS(1):root(1),root(2),root(3),2) = 0;
    end
    if(aprootS(2)>root(2))
      M(aprootS(1),root(2):aprootS(2),root(3),3) = 0;
    elseif(aprootS(2)<root(2))
      M(aprootS(1),aprootS(2):root(2),root(3),4) = 0;
    end
    if(aprootS(3)>root(3))
      M(aprootS(1),aprootS(2),root(3):aprootS(3),5) = 0;
    elseif(aprootS(3)<root(3))
      M(aprootS(1),aprootS(2),aprootS(3):root(3),6) = 0;
    end

    M(1,:,:,1)=inf;M(end,:,:,2)=inf;M(:,1,:,3)=inf;M(:,end,:,4)=inf;
    M(:,:,1,5)=inf;M(:,:,end,6)=inf;
    
    fillfrac = reshape(single(obstruction>0),dimsSS(1),dims(1),dimsSS(2),dims(2),...
                                                        dimsSS(3),dims(3));
    fillfrac = squeeze(sum(sum(sum(fillfrac,5),3),1))/prod(dimsSS);
    M = bsxfun(@plus,M,-bweight*log(1-(2*max(0,fillfrac-0.5))));
    
    
    M = reshape(M,prod(dims),6);
    [~,pathfrom] = dendrite_dijkstra2(M,dims,root);
  
    % Find endpoints
    endsT = zeros(numdt,3);
    for i = 1:numdt
        flag   = 1;
        distSC = 1;
        numit =  0;
        while(flag && numit<100)
            theta = rand*2*pi;                                             % Get a random orientation
            r = sqrt(rand)*dtParams(2)*distSC;                             % Get a random distance
            endT = floor([r.*cos(theta)+rootL(1);r.*sin(theta)+...
                             rootL(2);2*dtParams(3)*(rand-0.5)+rootL(3)])';% Propose an end-point based on the orientation and distance drawn
            if(endT(3)>fdims(3)); endT(3)=fdims(3); end                    %---
            if(endT(2)>fdims(2)); endT(2)=fdims(2); end                    % |
            if(endT(1)>fdims(1)); endT(1)=fdims(1); end                    % These lines make sure that the end-point is actually within the simulated volume
            if(endT(3)<1); endT(3)=1; end                                  % |
            if(endT(2)<1); endT(2)=1; end                                  % |
            if(endT(1)<1); endT(1)=1; end                                  %---
            if(obstruction(endT(1),endT(2),endT(3))==0)
                endsT(i,:) = endT;                                         % If a feasible end-point was found, save the end-point
                flag       = 0;                                            % Clear the flag to get out of the loop
            end
            distSC     = distSC*1.01;                                      % In the case of failure, increase the range that an end-point is looked for in
            numit      = numit+1;                                          % Incriment the number of attempts to get a feasible end-point
        end
    end
    endsTC = ceil(bsxfun(@rdivide,endsT,dimsSS));                           % Store the end points
    endsA  = zeros(atParams(1),3);                                         % Initialize the apical dendrite ends
    rootA  = floor([rootL(1)+2*atParams(4)*(rand-0.5),rootL(2)+2*...
                                    atParams(4)*(rand-0.5),1+atParams(3)]);
    if(size(rotAng,1)>=j)
        rootA2 = rootA+rootL(3)*sind([rotAng(j,2) -rotAng(j,1) 0]);
    else
        rootA2 = rootA;
    end
                                  
    for i = 1:atParams(1)
        flag   = 1;
        distSC = 1;
        numit =  0;
        while(flag && numit<100)
            theta = rand*2*pi;                                         % Pick a random direction to grow the dendrite in
            r     = sqrt(rand)*atParams(2)*distSC;                     % Pick on off-center position for the apical dendrite
            endA  = floor([r.*cos(theta)+rootA2(1);r.*sin(theta)+rootA2(2)...
                                      ;2*atParams(3)*(rand-0.5)+rootA2(3)])';% Pick the end-point of the apical dendrite
            if(endA(3)>fdims(3)); endA(3)=fdims(3); end                    %---
            if(endA(2)>fdims(2)); endA(2)=fdims(2); end                    % |
            if(endA(1)>fdims(1)); endA(1)=fdims(1); end                    % These lines make sure that the end-point is actually within the simulated volume
            if(endA(3)<1); endA(3)=1; end                                  % |
            if(endA(2)<1); endA(2)=1; end                                  % |
            if(endA(1)<1); endA(1)=1; end                                  %---
            if(obstruction(endA(1),endA(2),endA(3))==0)
                endsA(i,:) = endA;                                         % If a feasible end-point was found, save the end-point
                flag       = 0;                                            % Clear the flag to get out of the loop
            end
            distSC     = distSC*1.01;                                      % In the case of failure, increase the range that an end-point is looked for in
            numit      = numit+1;                                          % Incriment the number of attempts to get a feasible end-point     
        end
    end
    endsAC = ceil(bsxfun(@rdivide,endsA,dimsSS));
    ends   = [endsTC;endsAC];
    endsL  = [endsT;endsA];
    nends  = numdt + atParams(1);
  
    % Retrive paths
    paths = false(dims);
    for i = 1:nends
        path = getDendritePath2(pathfrom,ends(i,:),root);
        if(~isempty(path))
            paths(sub2ind(dims,path(:,1),path(:,2),path(:,3))) = true;
        end
    end    
    
    % Refine paths
    rootL   = round((root-[0.5 0.5 0.5]).*dimsSS);
    denLocs = find(paths);
    for i = 1:length(denLocs)
        temp       = 1+dweight*rand(dimsSS(1),dimsSS(2),dimsSS(3),6);      % Lattice of values, in order of R,L,U,D,T,B
        [lx,ly,lz] = ind2sub(dims,denLocs(i));
        filled = obstruction((lx-1)*dimsSS(1)+(1:dimsSS(1)),(ly-1)*dimsSS(2)+...
                   (1:dimsSS(2)),(lz-1)*dimsSS(3)+(1:dimsSS(3)))*inf;
        filled(isnan(filled))=0;
        temp       = bsxfun(@plus,temp,filled);
        ML((lx-1)*dimsSS(1)+(1:dimsSS(1)),(ly-1)*dimsSS(2)+...
                   (1:dimsSS(2)),(lz-1)*dimsSS(3)+(1:dimsSS(3)),:) = temp; %
    end
    aprootL = rootL+aproot-allroots(j,:);
    aprootL = max(aprootL,[1 1 1]);
    aprootL = min(aprootL,fdims);  
    if(aprootL(1)>rootL(1))
      ML(rootL(1):aprootL(1),rootL(2),rootL(3),1) = 0;
    elseif(aprootL(1)<rootL(1))
      ML(aprootL(1):rootL(1),rootL(2),rootL(3),2) = 0;
    end
    if(aprootL(2)>rootL(2))
      ML(aprootL(1),rootL(2):aprootL(2),rootL(3),3) = 0;
    elseif(aprootL(2)<rootL(2))
      ML(aprootL(1),aprootL(2):rootL(2),rootL(3),4) = 0;
    end
    if(aprootL(3)>rootL(3))
      ML(aprootL(1),aprootL(2),rootL(3):aprootL(3),5) = 0;
    elseif(aprootL(3)<rootL(3))
      ML(aprootL(1),aprootL(2),aprootL(3):rootL(3),6) = 0;
    end

    ML(1,:,:,1)=inf;ML(end,:,:,2)=inf;ML(:,1,:,3)=inf;ML(:,end,:,4)=inf;
    ML(:,:,1,5)=inf;ML(:,:,end,6)=inf;

    [~,pathfromL] = dendrite_dijkstra2(reshape(ML,prod(fdims),6),...
                                                             fdims,rootL); % Run dijksra's algorithm to get the path
    
    for i = 1:length(denLocs)
        [lx,ly,lz] = ind2sub(dims,denLocs(i));
        ML((lx-1)*dimsSS(1)+(1:dimsSS(1)),(ly-1)*dimsSS(2)+...
                   (1:dimsSS(2)),(lz-1)*dimsSS(3)+(1:dimsSS(3)),:) = inf; %
    end
    
    finepathsIdx = zeros(fdims,'single');
    fineIdxs = [];
    allpaths = cell(numdt+atParams(1),1);
    for i = 1:numdt
        path = getDendritePath2(pathfromL,endsL(i,:),rootL);
        allpaths{i} = path;
        if ~isempty(path)
            dendSz = max(0,normrnd(1,dendVar))^2;
            pathW = dendSz*single(1-(1-1/sqrt(2))*[0;sum(abs(diff(abs(diff(path)))),2)/2;0]);
            finepathsIdx(sub2ind(fdims,path(:,1),path(:,2),path(:,3))) ...
              = finepathsIdx(sub2ind(fdims,path(:,1),path(:,2),path(:,3)))+pathW; % 
            fineIdxs =[fineIdxs; sub2ind(fdims,path(:,1),path(:,2),path(:,3))];
        end
    end
    fineIdxs = unique(fineIdxs);
    finepathsIdx(fineIdxs) = thicknessScale*dtParams(4)*(finepathsIdx(fineIdxs).^(1/rallexp));
    
    finepathsVal = zeros(fdims,'single');
    fineIdxs2 = [];
    for i = 1:atParams(1)
        path = getDendritePath2(pathfromL,endsL(i+numdt,:),rootL);         % 
        allpaths{i+numdt} = path;
        if ~isempty(path)
            if(size(path,1)>2)
                dendSz = max(0,normrnd(1,dendVar));
                pathW  = dendSz*single(1-(1-1/sqrt(2))*[0;...
                                  sum(abs(diff(abs(diff(path)))),2)/2;0]); % Get second derivatives of dendrite path length
            else
                pathW = 1;                                                 % In the special case of only 1 or 2 nodes, make pathW a default
%                 pathW = [1; 1];                                            % In the special case of only 2 nodes, make pathW a default
            end
            try
            finepathsVal(sub2ind(fdims,path(:,1),path(:,2),path(:,3))) ...
            = finepathsVal(sub2ind(fdims,path(:,1),path(:,2),path(:,3)))+pathW; %
            catch ME
                disp(size(pathW))
                disp(size(finepathsVal(sub2ind(fdims,path(:,1),path(:,2),path(:,3)))))
                rethrow(ME)
            end
            fineIdxs2 = [fineIdxs2;sub2ind(fdims,path(:,1),path(:,2),path(:,3))];
        end
    end
    fineIdxs2 = unique(fineIdxs2);
    finepathsAD = false(fdims);
    finepathsAD(fineIdxs2) = true;
    finepathsVal(fineIdxs2) = thicknessScale*atParams(5)*(finepathsVal(fineIdxs2).^(1/rallexp));

    finepathsVal(fineIdxs) = finepathsIdx(fineIdxs)+finepathsVal(fineIdxs);% used to save on space
    
    finepathsIdx = zeros(fdims,'single');

    cellBody = cat(1,cellBody,sub2ind(fdims,rootL(1),rootL(2),rootL(3)));
    cellBodySmoothed = smoothCellBody(allpaths,cellBody,fdims);
    cellBodySmoothed = setxor(cellBodySmoothed,cellBody);
    fineIdxs3 = unique([fineIdxs;fineIdxs2;cellBodySmoothed]);
    finepathsIdx(fineIdxs3) = j;
    finepathsIdx(cellBodySmoothed)   = j;    %for some reason this is not working?
    finepathsVal(cellBodySmoothed) = finepathsVal(cellBodySmoothed)+1;
    
    finepathsIdx(cellBody)   = 0;    %for some reason this is not working?
    finepathsVal(cellBody) = 0;
    finepathsAD(cellBody) = 0;
    
    
    % Set matrix values    
    if ~borderflag
      if(smallZ)
        [xi,yi,zi] = ind2sub(fdims,fineIdxs3);
        xi = xi+allroots(j,1)-fdims(1)/2-1;
        yi = yi+allroots(j,2)-fdims(2)/2-1;
        Didxs = sub2ind(fulldims,xi,yi,zi);
        Fidxs = fineIdxs3;        
      else
        [xi,yi,zi] = ind2sub(fdims,fineIdxs3);
        xi = xi+allroots(j,1)-fdims(1)/2-1;
        yi = yi+allroots(j,2)-fdims(2)/2-1;
        zi = zi+allroots(j,3)-fdims(3)/2-1;
        Didxs = sub2ind(fulldims,xi,yi,zi);
        Fidxs = fineIdxs3;
      end
    else
      if(smallZ)
        [xi,yi,zi] = ind2sub(fdims,fineIdxs3);
        xi = xi+max([1 Xlims(1)])-max([1 -Xlims(1)+2]);
        yi = yi+max([1 Ylims(1)])-max([1 -Ylims(1)+2]);
        gidxs = (xi<=fulldims(1))&(yi<=fulldims(2))&(zi<=fulldims(3))& ...
                 (xi>=1)&(yi>=1)&(zi>=1);
        xi = xi(gidxs);
        yi = yi(gidxs);
        zi = zi(gidxs);
        Didxs = sub2ind(fulldims,xi,yi,zi);
        Fidxs = fineIdxs3(gidxs);
      else
        [xi,yi,zi] = ind2sub(fdims,fineIdxs3);
        xi = xi+max([1 Xlims(1)])-max([1 -Xlims(1)+2]);
        yi = yi+max([1 Ylims(1)])-max([1 -Ylims(1)+2]);
        zi = zi+max([1 Zlims(1)])-max([1 -Zlims(1)+2]);
        gidxs = (xi<=fulldims(1))&(yi<=fulldims(2))&(zi<=fulldims(3))& ...
                 (xi>=1)&(yi>=1)&(zi>=1);
        xi = xi(gidxs);
        yi = yi(gidxs);
        zi = zi(gidxs);
        Didxs = sub2ind(fulldims,xi,yi,zi);
        Fidxs = fineIdxs3(gidxs);
      end
    end
    
    cellVolume(Didxs) = cellVolume(Didxs) + finepathsIdx(Fidxs);
    cellVolumeIdx(Didxs) = cellVolumeIdx(Didxs) + finepathsIdx(Fidxs);
    cellVolumeVal(Didxs) = cellVolumeVal(Didxs) + finepathsVal(Fidxs);
    cellVolumeAD(Didxs) = cellVolumeAD(Didxs) + finepathsAD(Fidxs);
    
    if(nargout>3)
      gp_soma{j,2} = cellBodySmoothed;
    end
    if vol_params.verbose == 1
        fprintf('.');
    elseif vol_params.verbose >1
        Tdone = toc;
        fprintf('%d (%f seconds).\n',j,Tdone);
    end
end

%%
cellVolumeVal = uint16(floor(cellVolumeVal)+(mod(cellVolumeVal,1)>rand(size(cellVolumeVal))));
% cellVolumeVal = uint16(ceil(cellVolumeVal));
cellVolumeIdx = uint16(cellVolumeIdx);
cellVolumeAD = uint16(cellVolumeAD);
cellVolumeBD = uint16(~cellVolumeAD);
% before dilating, turns in dendrites should be weighted down to reduce the
% importance of edge effects.
% [~,dendnum] = dilateDendritePathAll(cellVolumeVal,cellVolumeIdx,neur_num);

[~,dendnumAD] = dilateDendritePathAll(cellVolumeVal.*cellVolumeAD,cellVolumeIdx.*cellVolumeAD,neur_num);
[~,dendnumBD] = dilateDendritePathAll(cellVolumeVal.*cellVolumeBD,cellVolumeIdx.*cellVolumeBD,neur_num);
for kk = 1:N_neur       
    dendnumAD(gp_nuc{kk,1})     = uint16(0);                                % Remove nuceus growths from neur_num
    dendnumBD(gp_nuc{kk,1})     = uint16(0);                                % Remove nuceus growths from neur_num
    dendnumAD(gp_soma{kk,1})    = uint16(0);                               % Make sure somas are STILL where they are supposed to be!!
    dendnumBD(gp_soma{kk,1})    = uint16(0);                               % Make sure somas are STILL where they are supposed to be!!
end

dendnumBD(dendnumAD>0) = dendnumAD(dendnumAD>0);
neur_num(dendnumBD>0)  = dendnumBD(dendnumBD>0);                               

% Remove any dendrites that grew into the nucleus
for kk = 1:N_neur       
    neur_num(gp_nuc{kk,1})     = uint16(0);                                % Remove nuceus growths from neur_num
    neur_num(gp_soma{kk,1})    = uint16(kk);                               % Make sure somas are STILL where they are supposed to be!!
end 

if vol_params.verbose >= 1
    fprintf('done.\n')
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%