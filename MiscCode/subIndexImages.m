function subIm = subIndexImages(cnmf,suite2p, pcaica, est, Idxs)

if isempty(Idxs)
    Idxs = 1;
end

TMP_sz     = [500,500,1];
subIm.Idxs = Idxs(:);

subIm.cnmfIdx   = cnmf.pairs(any(bsxfun(@eq, cnmf.pairs(:,1), subIm.Idxs.'),2),2);
% CNMF - strongly paired components
subIm.cnmf      = cnmf.compSpatial(:,:,subIm.cnmfIdx);
subIm.cnmfbound = multiBoundaryImage(subIm.cnmf);
subIm.cnmf      = mean(bsxfun(@rdivide,subIm.cnmf,max(max(subIm.cnmf,[],1),[],2)),3);
subIm.cnmf      = subIm.cnmf./mean(subIm.cnmf(subIm.cnmf~=0));
subIm.cnmf      = cat(3, subIm.cnmf, zeros(TMP_sz), zeros(TMP_sz));

% PCA/ICA - strongly paired components
subIm.pcaicaIdx   = pcaica.pairs(any(bsxfun(@eq, pcaica.pairs(:,1), subIm.Idxs.'),2),2);
subIm.pcaica      = pcaica.compSpatialSc(:,:,subIm.pcaicaIdx);
subIm.pcaicabound = multiBoundaryImage(subIm.pcaica);
subIm.pcaica      = mean(bsxfun(@rdivide,subIm.pcaica,max(max(subIm.pcaica,[],1),[],2)),3);
subIm.pcaica      = subIm.pcaica./mean(subIm.pcaica(subIm.pcaica~=0));
subIm.pcaica      = cat(3, zeros(TMP_sz), subIm.pcaica, zeros(TMP_sz));

% Suite2p - strongly paired components
subIm.suite2pIdx   = suite2p.pairs(any(bsxfun(@eq, suite2p.pairs(:,1), subIm.Idxs.'),2),2);
subIm.suite2p      = suite2p.compSpatial(:,:,subIm.suite2pIdx);
subIm.suite2pbound = multiBoundaryImage(subIm.suite2p);
subIm.suite2p      = mean(bsxfun(@rdivide,subIm.suite2p,max(max(subIm.suite2p,[],1),[],2)),3);
subIm.suite2p      = subIm.suite2p./mean(subIm.suite2p(subIm.suite2p~=0));
subIm.suite2p      = cat(3, zeros(TMP_sz), zeros(TMP_sz), subIm.suite2p);

subIm.ideals     = bsxfun(@rdivide,est.compsIdealAll,max(max(est.compsIdealAll,[],1),[],2));
subIm.ideals     = subIm.ideals(:,:,subIm.Idxs);
subIm.idealbound = multiBoundaryImage(subIm.ideals);
subIm.ideals     = max(subIm.ideals,[],3);

%     ideals.strongcomps = ideals.strongcomps(:,:,ideals.strongpair);
%     ideals.strongbound = multiBoundaryImage(ideals.strongcomps);
%     ideals.strongcomps = max(ideals.strongcomps,[],3);

end

