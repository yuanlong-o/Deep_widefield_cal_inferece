function roc = genAlgorithmPairROC(cnmf,pcaica,suite2p,est,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 4
    N = varargin{1};
else
    N = 0.5;
end

if nargin > 5
    cutoff = varargin{2};
else
    cutoff = 0.5;
end

if isempty(cutoff)
    cutoff = 0.5;
end

if isempty(N)
    N = 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

ord_thresh = linspace(0,1,N);

cnmfTT    = bsxfun(@times, vec(sum(sum(cnmf.compSpatial))),   cnmf.compTimecourse);
pcaicaTT  = bsxfun(@times, vec(sum(sum(pcaica.compSpatial))), pcaica.compTimecourse);
estTT     = bsxfun(@times, vec(sum(sum(est.compsIdeal))),     est.estact(1:end-1,:));
suite2pTT = bsxfun(@times, vec(sum(sum(suite2p.compSpatial>0))), bsxfun(@plus,  -median(suite2p.compTimecourse,2), suite2p.compTimecourse));
% suite2pTT = bsxfun(@times, vec(sum(sum(suite2p.compSpatial>0))), suite2p.compTimecourse);
% suite2pTT = bsxfun(@plus,  -median(suite2pTT,2), suite2pTT);

roc.max.cnmf      = max(max(cnmfTT));
roc.max.suite2p   = max(max(suite2pTT));
roc.max.pcaica    = max(max(pcaicaTT));
roc.max.est       = max(max(estTT));

roc.min.cnmf      = min(max(cnmfTT,[],2));
roc.min.suite2p   = min(max(suite2pTT,[],2));
roc.min.pcaica    = min(max(pcaicaTT,[],2));
roc.min.est       = min(max(estTT,[],2));


roc.Ntrue.cnmf    = sum(cnmf.corrvals>cutoff);
roc.Ntrue.suite2p = sum(suite2p.corrvals>cutoff);
roc.Ntrue.pcaica  = sum(pcaica.corrvals>cutoff);
roc.Ntrue.est     = sum(est.corrvals>cutoff);

roc.Nflse.cnmf    = size(cnmfTT,1) - sum(cnmf.corrvals>cutoff);
roc.Nflse.suite2p = size(suite2pTT,1) - sum(suite2p.corrvals>cutoff);
roc.Nflse.pcaica  = size(pcaicaTT,1) - sum(pcaica.corrvals>cutoff);
roc.Nflse.est     = size(estTT,1) - sum(est.corrvals>cutoff);

roc.FA.cnmf    = zeros(numel(ord_thresh),1);
roc.FA.suite2p = zeros(numel(ord_thresh),1);
roc.FA.pcaica  = zeros(numel(ord_thresh),1);
roc.FA.est     = zeros(numel(ord_thresh),1);

roc.TP.cnmf    = zeros(numel(ord_thresh),1);
roc.TP.suite2p = zeros(numel(ord_thresh),1);
roc.TP.pcaica  = zeros(numel(ord_thresh),1);
roc.TP.est     = zeros(numel(ord_thresh),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate ROC values

for ll = 1:numel(ord_thresh)
    % Pick time traces above a certain threshold
    roc.pick.cnmf    = find(max(cnmfTT,[],2)    >= roc.min.cnmf    + ord_thresh(ll)*(roc.max.cnmf-roc.min.cnmf));
    roc.pick.suite2p = find(max(suite2pTT,[],2) >= roc.min.suite2p + ord_thresh(ll)*(roc.max.suite2p-roc.min.suite2p));
    roc.pick.pcaica  = find(max(pcaicaTT,[],2)  >= roc.min.pcaica  + ord_thresh(ll)*(roc.max.pcaica-roc.min.pcaica));
    roc.pick.est     = find(max(estTT,[],2)     >= roc.min.est     + ord_thresh(ll)*(roc.max.est-roc.min.est));
    
    cnmf_HIT    = intersect(roc.pick.cnmf,    cnmf.pairs(cnmf.corrvals>cutoff,2));
    suite2p_HIT = intersect(roc.pick.suite2p, suite2p.pairs(suite2p.corrvals>cutoff,2));
    pcaica_HIT  = intersect(roc.pick.pcaica,  pcaica.pairs(pcaica.corrvals>cutoff,2));
    est_HIT     = intersect(roc.pick.est,     est.pairs(est.corrvals>cutoff,2));
    
    roc.TP.cnmf(ll)    = numel(cnmf_HIT)/roc.Ntrue.cnmf;
    roc.TP.suite2p(ll) = numel(suite2p_HIT)/roc.Ntrue.suite2p;
    roc.TP.pcaica(ll)  = numel(pcaica_HIT)/roc.Ntrue.pcaica;
    roc.TP.est(ll)     = numel(est_HIT)/roc.Ntrue.est;
    
    roc.FA.cnmf(ll)    = numel(setdiff(roc.pick.cnmf,cnmf_HIT))/roc.Nflse.cnmf;
    roc.FA.suite2p(ll) = numel(setdiff(roc.pick.suite2p,suite2p_HIT))/roc.Nflse.suite2p;
    roc.FA.pcaica(ll)  = numel(setdiff(roc.pick.pcaica,pcaica_HIT))/roc.Nflse.pcaica;
    roc.FA.est(ll)     = numel(setdiff(roc.pick.est,est_HIT))/roc.Nflse.est;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%