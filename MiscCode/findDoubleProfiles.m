function ideals = findDoubleProfiles(ideals, cnmf, suite2p, pcaica)

TMP = unique(cnmf.allpairs.strongpairs(:,1));
ideals.doubling.cnmf = find(hist(cnmf.allpairs.strongpairs(:,1),TMP)>1);
ideals.doubling.cnmf = TMP(ideals.doubling.cnmf);
TMP = unique(cnmf.allpairs.strongpairs(:,1));
ideals.doubling.suite2p = find(hist(suite2p.allpairs.strongpairs(:,1),TMP)>1);
ideals.doubling.suite2p = TMP(ideals.doubling.suite2p);
TMP = unique(cnmf.allpairs.strongpairs(:,1));
ideals.doubling.pcaica  = find(hist(pcaica.allpairs.strongpairs(:,1),TMP)>1);
ideals.doubling.pcaica = TMP(ideals.doubling.pcaica);

end

