function [neur_num,neur_num_AD,dend_params] = growApicalDendrites(vol_params, dend_params, neur_num, cellVolumeAD, gp_nuc, gp_soma)

% [neur_num,neur_num_AD,dend_params] = growApicalDendrites(vol_params, dend_params, neur_num, cellVolumeAD, gp_nuc, neur_soma)
% 
% Function to grow dendrites in a neural volume. The inputs to this
% function are:
%   - vol_params      - Struct with parameters for the volume generation
%       .vol_sz       - 3-element vector with the size (in um) of the 
%                       volume to generate (default = 100x100x30um)
%       .min_dist     - Minimum distance between neurons (default = 15um)
%       .N_neur       - Number of neurons to generate (default = 50)
%       .vres         - resolution to simulate volume at (default = 2
%                       samples/um)
%       .N_den        - Width of apical dendrites (default = 10)
%       .N_bg         - Number of background/neuropil components to 
%                       simulate (default = 50)
%       .vol_depth    - Depth of the volume under the brain surface
%                       (default = 100um)
%       .dendrite_tau - Dendrite decay strength exponential distance
%                       (default = 5)
%       .verbose      - Level of verbosity in the output during the volume
%                       generation. Can be 0,1,2. 0 = no text updates, 1 =
%                       some text outputs. 2 = detailed text outputs
%                       (default = 1)
%   - dend_params - Struct containing parameters for dendrite simulation
%       .dtParams        - dendritic tree number,radius in um of branches
%                          (uniform across circle),radius in um of branches
%                          (uniform in z) (default = [40 150 50 1 10])
%       .atParams        - Apical dendrite number,radius in um of branches
%                          (uniform across circle),radius in um of branches
%                          (uniform in z),offset from center in um (default
%                          = [1 5 2 2 4]) 
%       .atParams2       - Through-volume apical dendrite number,radius in
%                          um of branches (uniform across circle),radius in
%                          um of branches (uniform in z),offset from center
%                          in um (default = = [1 5 2 2 6])
%       .dweight         - Weight for path planning randomness in the
%                          dendrites (default = 10) 
%       .bweight         - Weight for obstruction (default = 50)
%       .thicknessScale  - Scaling for dendrite thickness (int). Should be
%                          1 for 1um sampling,(4 for 0.5um sampling)
%                          (default = 0.75) 
%       .dims            - dims set at 10 um/space (default = [30 30 30])
%       .dimsSS          - dims subsampling factor (10 samples per dim
%                          grid) (default = [10 10 10]) 
%       .rallexp         - Rall exponent that controls the cahnge in size
%                          of dendrites at branching locations (default =
%                          1.5) 
%   - neur_num     - An array where the k^th neuron's locations in the
%                    volume (both the soma and dendrites) are deliniated by
%                    the value 'k'
%   - cellVolumeAD - 
%   - gp_nuc       - 
%   - neur_soma    - vol_params.vol_sz-sized array where the integers at 
%                    each location indicate which neuron exists in that
%                    voxel 
% 
%  The outputs for this function are:
%   - neur_num     - Updated array where the k^th neuron's locations in the
%                    volume (both the soma and dendrites) are deliniated by
%                    the value 'k'
%   - neur_num_AD  - Volume array containing the number for each apical
%                    dendrite at each location it occupies
%   - dend_params  - The updated struct of dendrite parameters (in case the
%                    default values were added)
% 
% 2017 - Alex Song and Adam Charles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing

vol_params  = check_vol_params(vol_params);                                % Check volume parameters
dend_params = check_dend_params(dend_params);                              % Check dendrite parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grow Apical Dendrites

if vol_params.verbose == 1
    fprintf('Growing out apical dendrites')
elseif vol_params.verbose >1
    fprintf('Growing out apical dendrites...\n')
end

atParams       = dend_params.atParams2;
dweight        = dend_params.dweight;
bweight        = dend_params.bweight;
thicknessScale = dend_params.thicknessScale;
dims           = dend_params.dims;
dimsSS         = dend_params.dimsSS;
rallexp        = dend_params.rallexp;
vres           = vol_params.vres;                                          % Extract the volume resolution variable (points per micron)
N_neur         = vol_params.N_neur;                                        % Extract the [total] number of neurons
N_den          = vol_params.N_den;                                         % Extract number of dendrites
vol_sz         = vol_params.vol_sz;                                        % Extract volume size
dims           = min(dims, [vol_sz(1), vol_sz(2), vol_sz(3)]./dimsSS);     % dims set at 10 um per space
fulldims       = vol_sz*vres;                                              % Get dimension of the volume to work with (sampling times size)
dims           = dims*vres;                                                % Calculate dimensions of the sampled array
atParams(2:4)  = atParams(2:4)*vres;
thicknessScale = thicknessScale*vres*vres;                                 % Calculate thicknes of dendrites
fdims          = dims.*dimsSS;
fdims          = min(fdims, fulldims);

cellVolume     = single(neur_num);                                         % Make sure everying is single precision for memory's sake
for kk = 1:N_neur
    cellVolume(gp_nuc{kk,1}) = kk;
end
root_den       = zeros(N_den,3);                                           % Initialize an array for the roots of the dendrites

for j = 1:N_den                                                             % set up start and end points for apical dendrites
    while(~root_den(j,3))
        root = ceil(fulldims(1:2).*rand(1,2));
        if(cellVolume(root(1),root(2),1)==0)
            root_den(j,:) = [root fulldims(3)];
        end
    end
end

cellVolumeIdx   = zeros(size(neur_num), 'single');
cellVolumeVal   = zeros(size(neur_num), 'single');                         % Initialize the volume of apical dendrites
ML = inf([fdims 6]);

if (~isfield(dend_params,'apicalVar'))||isempty(dend_params.apicalVar)
  if (~isfield(dend_params,'dendVar'))||isempty(dend_params.dendVar)  
    dendVar = 0.35;
  else
    dendVar = dend_params.dendVar;
  end
else
  dendVar = dend_params.apicalVar;
end


for j = 1:N_den
    if vol_params.verbose >1
        tic
    end
    
    rootL = [(fdims(1)/2+1) (fdims(2)/2+1) root_den(j,3)]; % Center origin point (in um)
    % Create paths
    borderflag = 0;
    try
        obstruction = cellVolume(root_den(j,1)+(1:fdims(1))-fdims(1)/2-1,root_den(j,2)+(1:fdims(2))-fdims(2)/2-1,:);
    catch
        obstruction = zeros(fdims,'single');
        Xlims = [(root_den(j,1)-fdims(1)/2) (root_den(j,1)+fdims(1)/2-1)];
        Ylims = [(root_den(j,2)-fdims(2)/2) (root_den(j,2)+fdims(2)/2-1)];
        obstruction(max([1 -Xlims(1)+2]):(fdims(1)+min([0 -Xlims(2)+fulldims(1)])),max([1 -Ylims(1)+2]):(fdims(2)+min([0 -Ylims(2)+fulldims(2)])),:) ...
          = cellVolume(max([1 Xlims(1)]):min([Xlims(2) fulldims(1)]),max([1 Ylims(1)]):min([Ylims(2) fulldims(2)]),:);
        borderflag = 1;
    end
    cellBody = (obstruction==j+N_neur);
    obstruction(cellBody) = 0;

    root = ceil(rootL./dimsSS);
    M = 1+dweight*rand(dims(1),dims(2),dims(3),6);                         % Lattice of values, in order of R,L,U,D,T,B
    M(1,:,:,1)=inf;M(end,:,:,2)=inf;M(:,1,:,3)=inf;M(:,end,:,4)=inf;M(:,:,1,5)=inf;M(:,:,end,6)=inf;
  
    fillfrac = reshape(single(obstruction>0),dimsSS(1),dims(1),dimsSS(2),dims(2),dimsSS(3),dims(3));
    fillfrac = squeeze(sum(sum(sum(fillfrac,5),3),1))/prod(dimsSS);
    M = bsxfun(@plus,M,-bweight*log(1-(2*max(0,fillfrac-0.5))));
    
    M = reshape(M,prod(dims),6);
    [~,pathfrom] = dendrite_dijkstra2(M,dims,root);
  
    % Find endpoints
    endsA = zeros(atParams(1),3);
    rootA = floor([rootL(1)+2*atParams(4)*(rand-0.5),rootL(2)+2*atParams(4)*(rand-0.5),fdims(3)]);
    for i = 1:atParams(1)
        flag = 1;
        distSC = 1;
        numit =  0;
        while(flag && numit<100)
            theta = rand*2*pi;                                             % Get a random orientation
            r = sqrt(rand)*atParams(2)*distSC;                             % Get a random distance
            endA = floor([r.*cos(theta)+rootA(1);(r.*sin(theta)...
                                                         +rootA(2));1])';  % Propose an end-point based on the orientation and distance drawn
            if(endA(3)>fdims(3));endA(3)=fdims(3);end                      %---
            if(endA(2)>fdims(2));endA(2)=fdims(2);end                      % |
            if(endA(1)>fdims(1));endA(1)=fdims(1);end                      % These lines make sure that the end-point is actually within the simulated volume
            if(endA(3)<1);endA(3)=1;end                                    % |
            if(endA(2)<1);endA(2)=1;end                                    % |
            if(endA(1)<1);endA(1)=1;end                                    %---
            if(obstruction(endA(1),endA(2),endA(3))==0)
                endsA(i,:) = endA;                                         % If a feasible end-point was found, save the end-point
                flag       = 0;                                            % Clear the flag to get out of the loop
            end
            distSC     = distSC*1.01;                                      % In the case of failure, increase the range that an end-point is looked for in
            numit      = numit+1;                                          % Incriment the number of attempts to get a feasible end-point
        end
    end
    endsAC = ceil(bsxfun(@rdivide,endsA,dimsSS));

    % Retrive paths
    paths = zeros(dims);
    for i = 1:atParams(1)
        [~,pathM] = getDendritePath2(pathfrom,endsAC(i,:),root);
        paths = paths+pathM;
    end
  
    % Refine paths
    ML = inf*ML;
    rootL = round((root-[0.5 0.5 0.5]).*dimsSS);
    denLocs = find(paths);
    for i = 1:length(denLocs)
        temp = 1+dweight*rand(dimsSS(1),dimsSS(2),dimsSS(3),6);            % Lattice of values, in order of R,L,U,D,T,B
        [lx,ly,lz] = ind2sub(dims,denLocs(i));
        ML((lx-1)*dimsSS(1)+(1:dimsSS(1)),(ly-1)*dimsSS(2)+...
                    (1:dimsSS(2)),(lz-1)*dimsSS(3)+(1:dimsSS(3)),:) = temp;
    end
    ML(1,:,:,1)=inf;ML(end,:,:,2)=inf;ML(:,1,:,3)=inf;ML(:,end,:,4)=inf;
    ML(:,:,1,5)=inf;ML(:,:,end,6)=inf;
    filled = obstruction*inf;
    filled(isnan(filled))=0;
    ML = bsxfun(@plus,ML,filled);
    [~,pathfromL] = dendrite_dijkstra2(reshape(ML,prod(fdims),6),fdims,...
                                                                    rootL);% Run dijksra's algorithm to get the path

    finepathsVal = zeros(fdims,'single');
    for i = 1:atParams(1)
        path = getDendritePath2(pathfromL,endsA(i,:),rootL);
        if ~isempty(path)
            dendSz = max(0,normrnd(1,dendVar));
            pathW = dendSz*single(1-(1-1/sqrt(2))*[0;sum(abs(diff(abs(diff(path)))),2)/2;0]);
            finepathsVal(sub2ind(fdims,path(:,1),path(:,2),path(:,3))) ...
              = finepathsVal(sub2ind(fdims,path(:,1),path(:,2),path(:,3)))+pathW; % 
        end
    end
    finepathsVal(finepathsVal>0) = thicknessScale*atParams(5)*(finepathsVal(finepathsVal>0).^(1/rallexp));
    finepathsVal(cellBody) = 0;
    finepathsIdx = (j+N_neur)*single(finepathsVal>0);
    
    
    % Set matrix values
    if ~borderflag
        cellVolume(root_den(j,1)+(1:fdims(1))-fdims(1)/2-1,root_den(j,2)...
            +(1:fdims(2))-fdims(2)/2-1,:) = cellVolume(root_den(j,1)+...
            (1:fdims(1))-fdims(1)/2-1,root_den(j,2)+(1:fdims(2))-...
                                                 fdims(2)/2-1,:)+finepathsIdx;
        cellVolumeVal(root_den(j,1)+(1:fdims(1))-fdims(1)/2-1,root_den(j,2)...
            +(1:fdims(2))-fdims(2)/2-1,:) = cellVolumeVal(root_den(j,1)+...
            (1:fdims(1))-fdims(1)/2-1,root_den(j,2)+(1:fdims(2))-...
                                                 fdims(2)/2-1,:)+finepathsVal;
        cellVolumeIdx(root_den(j,1)+(1:fdims(1))-fdims(1)/2-1,root_den(j,2)...
            +(1:fdims(2))-fdims(2)/2-1,:) = cellVolumeIdx(root_den(j,1)+...
            (1:fdims(1))-fdims(1)/2-1,root_den(j,2)+(1:fdims(2))-...
                                                 fdims(2)/2-1,:)+finepathsIdx;
    else
        cellVolume(max([1 Xlims(1)]):min([Xlims(2) fulldims(1)]),...
            max([1 Ylims(1)]):min([Ylims(2) fulldims(2)]),:) = ...
            cellVolume(max([1 Xlims(1)]):min([Xlims(2) fulldims(1)]),...
            max([1 Ylims(1)]):min([Ylims(2) fulldims(2)]),:)+ ...
            finepathsIdx(max([1 -Xlims(1)+2]):(fdims(1)+min([0 -Xlims(2)+...
            fulldims(1)])),max([1 -Ylims(1)+2]):(fdims(2)+...
                                        min([0 -Ylims(2)+fulldims(2)])),:);
        cellVolumeVal(max([1 Xlims(1)]):min([Xlims(2) fulldims(1)]),...
            max([1 Ylims(1)]):min([Ylims(2) fulldims(2)]),:) = ...
            cellVolumeVal(max([1 Xlims(1)]):min([Xlims(2) fulldims(1)]),...
            max([1 Ylims(1)]):min([Ylims(2) fulldims(2)]),:)+ ...
            finepathsVal(max([1 -Xlims(1)+2]):(fdims(1)+min([0 -Xlims(2)+...
            fulldims(1)])),max([1 -Ylims(1)+2]):(fdims(2)+...
                                        min([0 -Ylims(2)+fulldims(2)])),:);
        cellVolumeIdx(max([1 Xlims(1)]):min([Xlims(2) fulldims(1)]),...
            max([1 Ylims(1)]):min([Ylims(2) fulldims(2)]),:) = ...
            cellVolumeIdx(max([1 Xlims(1)]):min([Xlims(2) fulldims(1)]),...
            max([1 Ylims(1)]):min([Ylims(2) fulldims(2)]),:)+ ...
            finepathsIdx(max([1 -Xlims(1)+2]):(fdims(1)+min([0 -Xlims(2)+...
            fulldims(1)])),max([1 -Ylims(1)+2]):(fdims(2)+...
                                        min([0 -Ylims(2)+fulldims(2)])),:);
    end
    if vol_params.verbose == 1
        fprintf('.')
    elseif vol_params.verbose >1
        Tdone = toc;
        fprintf('%d (%f seconds).\n',j,Tdone);
    end
end

% Dilate paths
cellVolumeVal = single(ceil(cellVolumeVal));
cellVolumeIdx = single(cellVolumeIdx);
[~,dendnum] = dilateDendritePathAll(cellVolumeVal,cellVolumeIdx,neur_num);


cellVolumeAD = uint16(cellVolumeAD)+uint16(dendnum);
neur_num     = uint16(neur_num)+uint16(dendnum);

for kk = 1:N_neur       
    neur_num(gp_nuc{kk,1})     = uint16(0);                                % Remove nucleus growths from neur_num
end

for kk = 1:N_neur       
    neur_num(gp_soma{kk,1})    = uint16(kk);                               % Make sure somas are STILL where they are supposed to be!!
end 

neur_num_AD  = uint16(cellVolumeAD);

for kk = 1:N_neur       
    neur_num_AD(gp_nuc{kk,1})     = uint16(0);                             % Remove nucleus growths from neur_num_AD
    neur_num_AD(gp_soma{kk,1})     = uint16(0);                            
end
 neur_num_AD((neur_num_AD-neur_num)>0) = uint16(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TESTING CODE: REMOVE IN FINAL VERSION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for kk = 1:max(neur_soma(:))
%     TMP = sum((neur_soma(:)==kk)&(neur_num(:)~=kk)); 
%     if TMP~=0 
%         disp(kk)
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if vol_params.verbose >= 1
    fprintf('done.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%