function ideals = separateIdealByCorrelation(ideals,est,cnmf,pcaica,suite2p,cuttoff)

if isempty(cuttoff)
    cuttoff = 0.5;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get pairings

% Get pairing infor for CNMF
cnmf.allpairs.strongpairs    = cnmf.pairs(cnmf.corrvals>=cuttoff,:);
cnmf.allpairs.weakpairs      = cnmf.pairs(cnmf.corrvals<cuttoff,:);
cnmf.allpairs.unpaired       = setdiff(1:size(cnmf.compSpatial,3),cnmf.pairs(:,2));
% Get pairing infor for PCA/ICA
pcaica.allpairs.strongpairs  = pcaica.pairs(pcaica.corrvals>=cuttoff,:);
pcaica.allpairs.weakpairs    = pcaica.pairs(pcaica.corrvals<cuttoff,:);
pcaica.allpairs.unpaired     = setdiff(1:size(pcaica.compSpatial,3),pcaica.pairs(:,2));
% Get pairing infor for Suite2p
suite2p.allpairs.strongpairs = suite2p.pairs(suite2p.corrvals>=cuttoff,:);
suite2p.allpairs.weakpairs   = suite2p.pairs(suite2p.corrvals<cuttoff,:);
suite2p.allpairs.unpaired    = setdiff(1:size(suite2p.compSpatial,3),suite2p.pairs(:,2));
% Get pairing infor for Ideals
est.allpairs.strongpairs = est.pairs(est.corrvals>=cuttoff,:);
est.allpairs.weakpairs   = est.pairs(est.corrvals<cuttoff,:);
est.allpairs.unpaired    = setdiff(1:size(est.compsIdeal,3),est.pairs(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract the ideal data into the appropriate format

ideals.strongpair  = unique([pcaica.allpairs.strongpairs(:,1);cnmf.allpairs.strongpairs(:,1);suite2p.allpairs.strongpairs(:,1)]);
ideals.weakpair    = setdiff(unique([pcaica.allpairs.weakpairs(:,1);cnmf.allpairs.weakpairs(:,1);suite2p.allpairs.weakpairs(:,1)]),ideals.strongpair);
ideals.nopair      = setdiff(1:sum(est.estactIdxs(1:size(est.compsIdealAll,3))),unique([ideals.strongpair;ideals.weakpair]));

% Strongly paried ideal components
ideals.strongcomps = est.compsIdealAll(:,:,ideals.strongpair);
ideals.strongcomps = bsxfun(@rdivide,ideals.strongcomps,max(max(ideals.strongcomps,[],1),[],2));
ideals.strongbound = multiBoundaryImage(ideals.strongcomps);
ideals.strongcomps = max(ideals.strongcomps,[],3);

% Weakly paried ideal components
ideals.weakcomps   = est.compsIdealAll(:,:,ideals.weakpair);
ideals.weakcomps   = bsxfun(@rdivide,ideals.weakcomps,max(max(ideals.weakcomps,[],1),[],2));
ideals.weakbound   = multiBoundaryImage(ideals.weakcomps);
ideals.weakcomps   = max(ideals.weakcomps,[],3);

% unparied ideal components
ideals.npcomps     = est.compsIdealAll(:,:,ideals.nopair);
ideals.npcomps     = bsxfun(@rdivide,ideals.npcomps,max(max(ideals.npcomps,[],1),[],2));
ideals.npbound     = multiBoundaryImage(ideals.npcomps);
ideals.npcomps     = max(ideals.npcomps,[],3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Also get images for ideal components not found by the algorithms

ideals.strongpairI = setdiff(est.allpairs.strongpairs(:,1), ideals.strongpair);
ideals.weakpairI   = setdiff(est.allpairs.weakpairs, union(setdiff(est.allpairs.weakpairs(:,1), ideals.weakpair),ideals.strongpairI));
ideals.strongpairX = setdiff(ideals.strongpair, est.allpairs.strongpairs(:,1));
ideals.weakpairX   = setdiff(ideals.weakpair, union(setdiff(ideals.weakpair, est.allpairs.weakpairs(:,1)),ideals.strongpairI));

% Strongly paried ideal components
ideals.strongcompsI = est.compsIdealAll(:,:,ideals.strongpairI);
ideals.strongcompsI = bsxfun(@rdivide,ideals.strongcompsI,max(max(ideals.strongcompsI,[],1),[],2));
ideals.strongboundI = multiBoundaryImage(ideals.strongcompsI);
ideals.strongcompsI = max(ideals.strongcompsI,[],3);

% Weakly paried ideal components
ideals.weakcompsI   = est.compsIdealAll(:,:,ideals.weakpairI);
ideals.weakcompsI   = bsxfun(@rdivide,ideals.weakcompsI,max(max(ideals.weakcompsI,[],1),[],2));
ideals.weakboundI   = multiBoundaryImage(ideals.weakcompsI);
ideals.weakcompsI   = max(ideals.weakcompsI,[],3);

% Strongly paried ideal components
ideals.strongcompsX = est.compsIdealAll(:,:,ideals.strongpairX);
ideals.strongcompsX = bsxfun(@rdivide,ideals.strongcompsX,max(max(ideals.strongcompsX,[],1),[],2));
ideals.strongboundX = multiBoundaryImage(ideals.strongcompsX);
ideals.strongcompsX = max(ideals.strongcompsX,[],3);

% Weakly paried ideal components
ideals.weakcompsX   = est.compsIdealAll(:,:,ideals.weakpairX);
ideals.weakcompsX   = bsxfun(@rdivide,ideals.weakcompsX,max(max(ideals.weakcompsX,[],1),[],2));
ideals.weakboundX   = multiBoundaryImage(ideals.weakcompsX);
ideals.weakcompsX   = max(ideals.weakcompsX,[],3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% OLDER CODE

% ideals.strongcomps = ideals.strongcomps(:,:,est.estactIdxs(1:size(ideals.strongcomps,3)));
% ideals.weakcomps   = ideals.weakcomps(:,:,est.estactIdxs(1:size(ideals.weakcomps,3)));
% ideals.npcomps     = ideals.npcomps(:,:,est.estactIdxs(1:size(ideals.npcomps,3)));

