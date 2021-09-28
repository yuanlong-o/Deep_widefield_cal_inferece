function [paths,pathnums] = dilateDendritePathAll(paths,pathnums,obstruction)

% [paths,pathnums] = dilateDendritePathAll(paths,pathnums,obstruction)
%
% Function to simultaneously dilate all the dendrite paths in the volume.
% This function looks for empty space around the established dendrite
% positions to allocate to the expanded dendrite size, within taxicab 3
% The inputs of this function are:
%
% - paths       - Full set of simulated paths
% - pathnums    - Corresponding cell number
% - obstruction - Occupied space in volume
%
% The ouputs of this function are
%
% - paths       - Updated set of simulated paths
% - pathnums    - Updated corresponding cell number
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxDist = 20; % maximum number of voxels that an object can be dilated

[x,y,z] = meshgrid(-maxDist:maxDist);
dists = (x.^2+y.^2+z.^2);
dsz = size(dists);
[dval,didx]= sort(dists(:));
dpos = find(diff(dval));
clear x y z dists dval


paths = single(paths);
paths(logical(obstruction)) = nan;
dims = size(paths);
pdims = prod(dims);
dshifts = [-dims(1)*dims(2);dims(1)*dims(2);-dims(1);dims(1);-1;1];

idxs = find(paths>1);
i = 1;
while(i<maxDist^2 && ~isempty(idxs))
  [dx,dy,dz] = ind2sub(dsz,didx(dpos(i)+1:dpos(i+1)));                     % shifts from center representing 
  dx = dx-maxDist-1;
  dy = dy-maxDist-1;
  dz = dz-maxDist-1;
  jidxs = (dz)*dims(1)*dims(2)-(dy)*dims(1)-dx;                            % shifts for tested indexes
  for j = 1:length(idxs)
    pidxs = idxs(j)+jidxs;                                                 % make sure indexes are not out of bounds
    pidxs = pidxs(pidxs>0);
    pidxs = pidxs(pidxs<pdims);
    pidxs = pidxs(paths(pidxs)==0);
    didxt = false(length(pidxs),1);
    numval = pathnums(idxs(j));
    for k = 1:length(pidxs)                                                % make sure that indexes are connected to rest of expanded process
      didxs = pidxs(k)+dshifts;
      didxs = didxs(didxs>0);
      didxs = didxs(didxs<pdims);
      if(sum(pathnums(didxs)==numval))
        didxt(k) = true;
      end
    end
    pidxs = pidxs(didxt);
    if(~isempty(pidxs))
      while(paths(idxs(j))>1 && ~isempty(pidxs))
        ridx = randi(length(pidxs));
        pidx = pidxs(ridx);
        pidxs(ridx) = [];
        paths(idxs(j)) = paths(idxs(j))-1;
        paths(pidx) = 1;
        pathnums(pidx) = numval;
      end
    end
  end
  idxs = find(paths>1);
  i = i+1;
end