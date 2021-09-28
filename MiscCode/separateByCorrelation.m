function [cnmf,suite2p,pcaica,est] = separateByCorrelation(cnmf,suite2p,pcaica,est,cuttoff)

% [cnmf,suite2p,pcaica,est] = separateByCorrelation(cnmf,suite2p,pcaica,est,cuttoff)
%
%
%
% 2018 - Adam Charles & Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing/initializations

if isempty(cuttoff)
    cuttoff = 0.5;
end

TMP_sz = [500,500,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get pairing info

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
%% Isolate the different components
% CNMF - strongly paired components
cnmf.allpairs.strongcomp = cnmf.compSpatial(:,:,cnmf.allpairs.strongpairs(:,2));
cnmf.strongbound         = multiBoundaryImage(cnmf.allpairs.strongcomp);
cnmf.allpairs.strongcomp = mean(bsxfun(@rdivide,cnmf.allpairs.strongcomp,max(max(cnmf.allpairs.strongcomp,[],1),[],2)),3);
cnmf.allpairs.strongcomp = cnmf.allpairs.strongcomp./mean(cnmf.allpairs.strongcomp(cnmf.allpairs.strongcomp~=0));
cnmf.allpairs.strongcomp = cat(3, cnmf.allpairs.strongcomp, zeros(TMP_sz), zeros(TMP_sz));
% CNMF - weakly paired components
cnmf.allpairs.weakcomp = cnmf.compSpatial(:,:,cnmf.allpairs.weakpairs(:,2));
cnmf.weakbound         = multiBoundaryImage(cnmf.allpairs.weakcomp);
cnmf.allpairs.weakcomp = mean(bsxfun(@rdivide,cnmf.allpairs.weakcomp,max(max(cnmf.allpairs.weakcomp,[],1),[],2)),3);
cnmf.allpairs.weakcomp = cnmf.allpairs.weakcomp./mean(cnmf.allpairs.weakcomp(cnmf.allpairs.weakcomp~=0));
cnmf.allpairs.weakcomp = cat(3, cnmf.allpairs.weakcomp, zeros(TMP_sz), zeros(TMP_sz));
% CNMF - unpaired components
cnmf.allpairs.upcomp = cnmf.compSpatial(:,:,cnmf.allpairs.unpaired(:));
cnmf.upbound         = multiBoundaryImage(cnmf.allpairs.upcomp);
cnmf.allpairs.upcomp = mean(bsxfun(@rdivide,cnmf.allpairs.upcomp,max(max(cnmf.allpairs.upcomp,[],1),[],2)),3);
cnmf.allpairs.upcomp = cnmf.allpairs.upcomp./mean(cnmf.allpairs.upcomp(cnmf.allpairs.upcomp~=0));
cnmf.allpairs.upcomp = cat(3, cnmf.allpairs.upcomp, zeros(TMP_sz), zeros(TMP_sz));

% PCA/ICA - strongly paired components
pcaica.allpairs.strongcomp = pcaica.compSpatialSc(:,:,pcaica.allpairs.strongpairs(:,2));
pcaica.strongbound         = multiBoundaryImage(pcaica.allpairs.strongcomp);
pcaica.allpairs.strongcomp = mean(bsxfun(@rdivide,pcaica.allpairs.strongcomp,max(max(pcaica.allpairs.strongcomp,[],1),[],2)),3);
pcaica.allpairs.strongcomp = pcaica.allpairs.strongcomp./mean(pcaica.allpairs.strongcomp(pcaica.allpairs.strongcomp~=0));
pcaica.allpairs.strongcomp = cat(3, zeros(TMP_sz), pcaica.allpairs.strongcomp, zeros(TMP_sz));
% PCA/ICA - weakly paired components
pcaica.allpairs.weakcomp = pcaica.compSpatialSc(:,:,pcaica.allpairs.weakpairs(:,2));
pcaica.weakbound         = multiBoundaryImage(pcaica.allpairs.weakcomp);
pcaica.allpairs.weakcomp = mean(bsxfun(@rdivide,pcaica.allpairs.weakcomp,max(max(pcaica.allpairs.weakcomp,[],1),[],2)),3);
pcaica.allpairs.weakcomp = pcaica.allpairs.weakcomp./mean(pcaica.allpairs.weakcomp(pcaica.allpairs.weakcomp~=0));
pcaica.allpairs.weakcomp = cat(3, zeros(TMP_sz), pcaica.allpairs.weakcomp, zeros(TMP_sz));
% PCA/ICA - unpaired components
pcaica.allpairs.upcomp = pcaica.compSpatialSc(:,:,pcaica.allpairs.unpaired(:));
pcaica.upbound         = multiBoundaryImage(pcaica.allpairs.upcomp);
pcaica.allpairs.upcomp = mean(bsxfun(@rdivide,pcaica.allpairs.upcomp,max(max(pcaica.allpairs.upcomp,[],1),[],2)),3);
pcaica.allpairs.upcomp = pcaica.allpairs.upcomp./mean(pcaica.allpairs.upcomp(pcaica.allpairs.upcomp~=0));
pcaica.allpairs.upcomp = cat(3, zeros(TMP_sz), pcaica.allpairs.upcomp, zeros(TMP_sz));

% Suite2p - strongly paired components
suite2p.allpairs.strongcomp = suite2p.compSpatial(:,:,suite2p.allpairs.strongpairs(:,2));
suite2p.strongbound         = multiBoundaryImage(suite2p.allpairs.strongcomp);
suite2p.allpairs.strongcomp = mean(bsxfun(@rdivide,suite2p.allpairs.strongcomp,max(max(suite2p.allpairs.strongcomp,[],1),[],2)),3);
suite2p.allpairs.strongcomp = suite2p.allpairs.strongcomp./mean(suite2p.allpairs.strongcomp(suite2p.allpairs.strongcomp~=0));
suite2p.allpairs.strongcomp = cat(3, zeros(TMP_sz), zeros(TMP_sz), suite2p.allpairs.strongcomp);
% Suite2p - weakly paired components
suite2p.allpairs.weakcomp = suite2p.compSpatial(:,:,suite2p.allpairs.weakpairs(:,2));
suite2p.weakbound         = multiBoundaryImage(suite2p.allpairs.weakcomp);
suite2p.allpairs.weakcomp = mean(bsxfun(@rdivide,suite2p.allpairs.weakcomp,max(max(suite2p.allpairs.weakcomp,[],1),[],2)),3);
suite2p.allpairs.weakcomp = suite2p.allpairs.weakcomp./mean(suite2p.allpairs.weakcomp(suite2p.allpairs.weakcomp~=0));
suite2p.allpairs.weakcomp = cat(3, zeros(TMP_sz), zeros(TMP_sz), suite2p.allpairs.weakcomp);
% Suite2p - unpaired components
suite2p.allpairs.upcomp = suite2p.compSpatial(:,:,suite2p.allpairs.unpaired(:));
suite2p.upbound         = multiBoundaryImage(suite2p.allpairs.upcomp);
suite2p.allpairs.upcomp = mean(bsxfun(@rdivide,suite2p.allpairs.upcomp,max(max(suite2p.allpairs.upcomp,[],1),[],2)),3);
suite2p.allpairs.upcomp = suite2p.allpairs.upcomp./mean(suite2p.allpairs.upcomp(suite2p.allpairs.upcomp~=0));
suite2p.allpairs.upcomp = cat(3, zeros(TMP_sz), zeros(TMP_sz), suite2p.allpairs.upcomp);

% Ideal - strongly paired components
est.allpairs.strongcomp = est.compsIdeal(:,:,est.allpairs.strongpairs(:,2));
est.strongbound         = multiBoundaryImage(est.allpairs.strongcomp);
est.allpairs.strongcomp = mean(bsxfun(@rdivide,est.allpairs.strongcomp,eps+max(max(est.allpairs.strongcomp,[],1),[],2)),3);
est.allpairs.strongcomp = est.allpairs.strongcomp./(eps+mean(est.allpairs.strongcomp(est.allpairs.strongcomp~=0)));
% Ideal - weakly paired components
est.allpairs.weakcomp   = est.compsIdeal(:,:,est.allpairs.weakpairs(:,2));
est.weakbound           = multiBoundaryImage(est.allpairs.weakcomp);
est.allpairs.weakcomp   = mean(bsxfun(@rdivide,est.allpairs.weakcomp,eps+max(max(est.allpairs.weakcomp,[],1),[],2)),3);
est.allpairs.weakcomp   = est.allpairs.weakcomp./(eps+mean(est.allpairs.weakcomp(est.allpairs.weakcomp~=0)));
% Ideal - unpaired components
est.allpairs.upcomp     = est.compsIdeal(:,:,est.allpairs.unpaired(:));
est.upbound             = multiBoundaryImage(est.allpairs.upcomp);
est.allpairs.upcomp     = mean(bsxfun(@rdivide,est.allpairs.upcomp,eps+max(max(est.allpairs.upcomp,[],1),[],2)),3);
est.allpairs.upcomp     = est.allpairs.upcomp./(eps+mean(est.allpairs.upcomp(est.allpairs.upcomp~=0)));


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%