function displayPairingStatistics(cnmf,pcaica,suite2p,est)



clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Single algorithm statistics

fprintf('Ideal:')
fprintf('Number componants: %d\n', size(est.estact,1))
fprintf('Number paired componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', sum(est.corrvals>0.1),sum(est.corrvals>0.3),sum(est.corrvals>0.5))
fprintf('Number unique componants found (>0.1,>0.3,>0.5): (%d,%d,%d)\n', ...
                         numel(unique(est.pairs(est.corrvals>0.1,1))), ...
                         numel(unique(est.pairs(est.corrvals>0.3,1))), ...
                         numel(unique(est.pairs(est.corrvals>0.5,1)))); 
fprintf('Number doubled componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', ...
    sum((hist(est.pairs(est.corrvals>0.1,1),unique(est.pairs(est.corrvals>0.1,1)))>1)),...
    sum((hist(est.pairs(est.corrvals>0.3,1),unique(est.pairs(est.corrvals>0.3,1)))>1)),...
    sum((hist(est.pairs(est.corrvals>0.5,1),unique(est.pairs(est.corrvals>0.5,1)))>1)))

fprintf('CNMF:')
fprintf('Number componants: %d\n', size(cnmf.compFluoresence,1))
fprintf('Number paired componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', sum(cnmf.corrvals>0.1),sum(cnmf.corrvals>0.3),sum(cnmf.corrvals>0.5))
fprintf('Number unique componants found (>0.1,>0.3,>0.5): (%d,%d,%d)\n', ...
                         numel(unique(cnmf.pairs(cnmf.corrvals>0.1,1))), ...
                         numel(unique(cnmf.pairs(cnmf.corrvals>0.3,1))), ...
                         numel(unique(cnmf.pairs(cnmf.corrvals>0.5,1)))); 
fprintf('Number doubled componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', ...
    sum((hist(cnmf.pairs(cnmf.corrvals>0.1,1),unique(cnmf.pairs(cnmf.corrvals>0.1,1)))>1)),...
    sum((hist(cnmf.pairs(cnmf.corrvals>0.3,1),unique(cnmf.pairs(cnmf.corrvals>0.3,1)))>1)),...
    sum((hist(cnmf.pairs(cnmf.corrvals>0.5,1),unique(cnmf.pairs(cnmf.corrvals>0.5,1)))>1)))


fprintf('Suite2p:')
fprintf('Number componants: %d\n', size(suite2p.compTimecourse,1))
fprintf('Number paired componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', sum(suite2p.corrvals>0.1),sum(suite2p.corrvals>0.3),sum(suite2p.corrvals>0.5))
fprintf('Number unique componants found (>0.1,>0.3,>0.5): (%d,%d,%d)\n', ...
                         numel(unique(suite2p.pairs(suite2p.corrvals>0.1,1))), ...
                         numel(unique(suite2p.pairs(suite2p.corrvals>0.3,1))), ...
                         numel(unique(suite2p.pairs(suite2p.corrvals>0.5,1)))); 
fprintf('Number doubled componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', ...
    sum((hist(suite2p.pairs(suite2p.corrvals>0.1,1),unique(suite2p.pairs(suite2p.corrvals>0.1,1)))>1)),...
    sum((hist(suite2p.pairs(suite2p.corrvals>0.3,1),unique(suite2p.pairs(suite2p.corrvals>0.3,1)))>1)),...
    sum((hist(suite2p.pairs(suite2p.corrvals>0.5,1),unique(suite2p.pairs(suite2p.corrvals>0.5,1)))>1)))

fprintf('PCA/ICA:')
fprintf('Number componants: %d\n', size(pcaica.compTimecourse,1))
fprintf('Number paired componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', sum(pcaica.corrvals>0.1),sum(pcaica.corrvals>0.3),sum(pcaica.corrvals>0.5))
fprintf('Number unique componants found (>0.1,>0.3,>0.5): (%d,%d,%d)\n', ...
                         numel(unique(pcaica.pairs(pcaica.corrvals>0.1,1))), ...
                         numel(unique(pcaica.pairs(pcaica.corrvals>0.3,1))), ...
                         numel(unique(pcaica.pairs(pcaica.corrvals>0.5,1)))); 
fprintf('Number doubled componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', ...
    sum((hist(pcaica.pairs(pcaica.corrvals>0.1,1),unique(pcaica.pairs(pcaica.corrvals>0.1,1)))>1)),...
    sum((hist(pcaica.pairs(pcaica.corrvals>0.3,1),unique(pcaica.pairs(pcaica.corrvals>0.3,1)))>1)),...
    sum((hist(pcaica.pairs(pcaica.corrvals>0.5,1),unique(pcaica.pairs(pcaica.corrvals>0.5,1)))>1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison statistics

TMPcs1 = intersect(unique(cnmf.pairs(cnmf.corrvals>0.1,1)), unique(suite2p.pairs(suite2p.corrvals>0.1,1)));
TMPcp1 = intersect(unique(cnmf.pairs(cnmf.corrvals>0.1,1)), unique(pcaica.pairs(pcaica.corrvals>0.1,1)));
TMPsp1 = intersect(unique(pcaica.pairs(pcaica.corrvals>0.1,1)), unique(suite2p.pairs(suite2p.corrvals>0.1,1)));
TMPxx1 = intersect(TMPcs1,intersect(TMPcp1,TMPsp1));

TMPcs3 = intersect(unique(cnmf.pairs(cnmf.corrvals>0.3,1)), unique(suite2p.pairs(suite2p.corrvals>0.3,1)));
TMPcp3 = intersect(unique(cnmf.pairs(cnmf.corrvals>0.3,1)), unique(pcaica.pairs(pcaica.corrvals>0.3,1)));
TMPsp3 = intersect(unique(pcaica.pairs(pcaica.corrvals>0.3,1)), unique(suite2p.pairs(suite2p.corrvals>0.3,1)));
TMPxx3 = intersect(TMPcs3,intersect(TMPcp3,TMPsp3));

TMPcs5 = intersect(unique(cnmf.pairs(cnmf.corrvals>0.5,1)), unique(suite2p.pairs(suite2p.corrvals>0.5,1)));
TMPcp5 = intersect(unique(cnmf.pairs(cnmf.corrvals>0.5,1)), unique(pcaica.pairs(pcaica.corrvals>0.5,1)));
TMPsp5 = intersect(unique(pcaica.pairs(pcaica.corrvals>0.5,1)), unique(suite2p.pairs(suite2p.corrvals>0.5,1)));
TMPxx5 = intersect(TMPcs5,intersect(TMPcp5,TMPsp5));

fprintf('Number of unique cells found by:\n')
fprintf('\t - Only CNMF (>0.1,>0.3,>0.5): %d, %d, %d\n', ...
    numel(setdiff(unique(cnmf.pairs(cnmf.corrvals>0.1,1)),union(TMPcs1,TMPcp1))), ...
    numel(setdiff(unique(cnmf.pairs(cnmf.corrvals>0.3,1)),union(TMPcs3,TMPcp3))), ...
    numel(setdiff(unique(cnmf.pairs(cnmf.corrvals>0.5,1)),union(TMPcs5,TMPcp5))))
fprintf('\t - Only Suite2p (>0.1,>0.3,>0.5): %d, %d, %d\n', ...
    numel(setdiff(unique(suite2p.pairs(suite2p.corrvals>0.1,1)),union(TMPcs1,TMPsp1))), ...
    numel(setdiff(unique(suite2p.pairs(suite2p.corrvals>0.3,1)),union(TMPcs3,TMPsp3))), ...
    numel(setdiff(unique(suite2p.pairs(suite2p.corrvals>0.5,1)),union(TMPcs5,TMPsp5))))
fprintf('\t - Only PCAICA (>0.1,>0.3,>0.5): %d, %d, %d\n', ...
    numel(setdiff(unique(pcaica.pairs(pcaica.corrvals>0.1,1)),union(TMPcs1,TMPsp1))), ...
    numel(setdiff(unique(pcaica.pairs(pcaica.corrvals>0.3,1)),union(TMPcs3,TMPsp3))), ...
    numel(setdiff(unique(pcaica.pairs(pcaica.corrvals>0.5,1)),union(TMPcs5,TMPsp5))))
fprintf('\t - CNMF and Suite2p (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPcs1), numel(TMPcs3), numel(TMPcs5))
fprintf('\t - CNMF and PCA/ICA (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPcp1), numel(TMPcp3), numel(TMPcp5))
fprintf('\t - PCA/ICA and Suite2p (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPsp1), numel(TMPsp3), numel(TMPsp5))
fprintf('\t - All algorithms (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPxx1), numel(TMPxx3), numel(TMPxx5))

%% Comparison to ideal

% TMPic1 = intersect(unique(est.pairs(est.corrvals>0.1,1)), unique(suite2p.pairs(suite2p.corrvals>0.1,1)));
% TMPis1 = intersect(unique(est.pairs(est.corrvals>0.1,1)), unique(suite2p.pairs(suite2p.corrvals>0.1,1)));
% TMPip1 = intersect(unique(est.pairs(est.corrvals>0.1,1)), unique(pcaica.pairs(pcaica.corrvals  >0.1,1)));
% TMPic3 = intersect(unique(est.pairs(est.corrvals>0.3,1)), unique(suite2p.pairs(suite2p.corrvals>0.3,1)));
% TMPis3 = intersect(unique(est.pairs(est.corrvals>0.3,1)), unique(suite2p.pairs(suite2p.corrvals>0.3,1)));
% TMPip3 = intersect(unique(est.pairs(est.corrvals>0.3,1)), unique(pcaica.pairs(pcaica.corrvals  >0.3,1)));
% TMPic5 = intersect(unique(est.pairs(est.corrvals>0.5,1)), unique(suite2p.pairs(suite2p.corrvals>0.5,1)));
% TMPis5 = intersect(unique(est.pairs(est.corrvals>0.5,1)), unique(suite2p.pairs(suite2p.corrvals>0.5,1)));
% TMPip5 = intersect(unique(est.pairs(est.corrvals>0.5,1)), unique(pcaica.pairs(pcaica.corrvals  >0.5,1)));

TMPinc1 = setdiff(unique(est.pairs(est.corrvals>0.1,1)), unique(cnmf.pairs(cnmf.corrvals      >0.1,1)));
TMPins1 = setdiff(unique(est.pairs(est.corrvals>0.1,1)), unique(suite2p.pairs(suite2p.corrvals>0.1,1)));
TMPinp1 = setdiff(unique(est.pairs(est.corrvals>0.1,1)), unique(pcaica.pairs(pcaica.corrvals  >0.1,1)));
TMPinc3 = setdiff(unique(est.pairs(est.corrvals>0.3,1)), unique(cnmf.pairs(cnmf.corrvals      >0.3,1)));
TMPins3 = setdiff(unique(est.pairs(est.corrvals>0.3,1)), unique(suite2p.pairs(suite2p.corrvals>0.3,1)));
TMPinp3 = setdiff(unique(est.pairs(est.corrvals>0.3,1)), unique(pcaica.pairs(pcaica.corrvals  >0.3,1)));
TMPinc5 = setdiff(unique(est.pairs(est.corrvals>0.5,1)), unique(cnmf.pairs(cnmf.corrvals      >0.5,1)));
TMPins5 = setdiff(unique(est.pairs(est.corrvals>0.5,1)), unique(suite2p.pairs(suite2p.corrvals>0.5,1)));
TMPinp5 = setdiff(unique(est.pairs(est.corrvals>0.5,1)), unique(pcaica.pairs(pcaica.corrvals  >0.5,1)));

fprintf('Number of unique cells found by the Ideal components:\n')
fprintf('\t - and NOT CNMF: (>0.1,>0.3,>0.5): %d, %d, %d\n',    numel(TMPinc1), numel(TMPinc3), numel(TMPinc5))
fprintf('\t - and NOT Suite2p: (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPins1), numel(TMPins3), numel(TMPins5))
fprintf('\t - and NOT PCA/ICA: (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPinp1), numel(TMPinp3), numel(TMPinp5))

TMPcni1 = setdiff(unique(cnmf.pairs(cnmf.corrvals      >0.1,1)), unique(est.pairs(est.corrvals>0.1,1)));
TMPsni1 = setdiff(unique(suite2p.pairs(suite2p.corrvals>0.1,1)), unique(est.pairs(est.corrvals>0.1,1)));
TMPpni1 = setdiff(unique(pcaica.pairs(pcaica.corrvals  >0.1,1)), unique(est.pairs(est.corrvals>0.1,1)));
TMPcni3 = setdiff(unique(cnmf.pairs(cnmf.corrvals      >0.3,1)), unique(est.pairs(est.corrvals>0.3,1)));
TMPsni3 = setdiff(unique(suite2p.pairs(suite2p.corrvals>0.3,1)), unique(est.pairs(est.corrvals>0.3,1)));
TMPpni3 = setdiff(unique(pcaica.pairs(pcaica.corrvals  >0.3,1)), unique(est.pairs(est.corrvals>0.3,1)));
TMPcni5 = setdiff(unique(cnmf.pairs(cnmf.corrvals      >0.5,1)), unique(est.pairs(est.corrvals>0.5,1)));
TMPsni5 = setdiff(unique(suite2p.pairs(suite2p.corrvals>0.5,1)), unique(est.pairs(est.corrvals>0.5,1)));
TMPpni5 = setdiff(unique(pcaica.pairs(pcaica.corrvals  >0.5,1)), unique(est.pairs(est.corrvals>0.5,1)));

fprintf('Number of unique cells not found by the Ideal components:\n')
fprintf('\t - but found by CNMF: (>0.1,>0.3,>0.5): %d, %d, %d\n',    numel(TMPcni1), numel(TMPcni3), numel(TMPcni5))
fprintf('\t - but found by Suite2p: (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPsni1), numel(TMPsni3), numel(TMPsni5))
fprintf('\t - but found by PCA/ICA: (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPpni1), numel(TMPpni3), numel(TMPpni5))


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
