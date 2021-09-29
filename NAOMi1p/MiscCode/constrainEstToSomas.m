function est2 = constrainEstToSomas(est,sL)

% est2 = constrainEstToSomas(est,sL)
% 
% Function to constrain the est data struct to a subset of profiles
% that represent, for example, somatic components.
% 
% 2020 - Adam Charles

tmpSL     = false(size(est.estactIdxs));
tmpSL(sL) = true;

%% First: indexing hell

est2.estactIdxs      = tmpSL&est.estactIdxs;                          % Full logical list of global components to keep
est2.estactidealIdxs = tmpSL&est.estactidealIdxs;                     % |- Same for ideal 

est2.Idxs            = nan(size(est2.estactIdxs));                    % Create nan list
est2.idealIdxs       = nan(size(est2.estactidealIdxs));               % |- Same for ideal
est2.Idxs(est2.estactIdxs)           = 1:sum(est2.estactIdxs);        % Populate with index of component in new list (for pairing purposes
est2.idealIdxs(est2.estactidealIdxs) = 1:sum(est2.estactidealIdxs);   % |- Same for ideal

%% Need indexing to subselect new components from previously wittled list

keepIdxs       = est2.Idxs(est.estactIdxs);
keepIdealIdxs  = est2.idealIdxs(est.estactidealIdxs);

est2.estact        = [est.estact(~isnan(keepIdxs),:); est.estact(end,:)];
est2.estactideal   = [est.estactideal(~isnan(keepIdealIdxs),:); est.estactideal(end,:)];
est2.compsIdealAll = est.compsIdealAll;
est2.compsIdeal    = est.compsIdeal(:,:,~isnan(keepIdxs));

est2.corrIdxs      = est.corrIdxs;
est2.corrvals      = est.corrvals;
est2.corrvalsIdeal = est.corrvalsIdeal;
est2.corrvalsIdeal(~isnan(est2.idealIdxs)) = NaN;
[est2.corrSortIdeal,IXC] = sort(est2.corrvalsIdeal,'descend');
est2.corrIdxsIdeal       = est2.idealIdxs(IXC);

% x - estact
% x - estactIdxs
% x - estactideal
% x - estactidealIdxs
% x - Idxs
% x - idealIdxs
% x - corrvals
% x - corrvalsIdeal
% x - corrSort
% x - corrIdxs
% x - corrSortIdeal
% x - corrIdxsIdeal
% x - compsIdealAll
% x - compsIdeal
% o - pairs 
% o - simTraces
% o - pairedTraces
% o - allpairs
% o - strongbound
% o - weakbound
% o - upbound


end
