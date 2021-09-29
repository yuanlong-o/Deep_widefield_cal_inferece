function output = smoothCellBody(allpaths,cellBody,fdims)

connIdxRoot = zeros(length(allpaths),3);
emptyIdxs = false(length(allpaths),1);
for i = 1:length(allpaths)
  path = allpaths{i};
  if(~isempty(path))
    pathInd = sub2ind(fdims,path(:,1),path(:,2),path(:,3));
    pathIntersect = ismember(pathInd,cellBody);
    try
      connIdxRoot(i,:) = path(find(pathIntersect,1),:);
    catch
      emptyIdxs(i) = 1;
    end
  else
    emptyIdxs(i) = 1;
  end
end

distMat = sqrt(bsxfun(@minus,connIdxRoot(:,1),connIdxRoot(:,1)').^2+...
  bsxfun(@minus,connIdxRoot(:,2),connIdxRoot(:,2)').^2+...
  bsxfun(@minus,connIdxRoot(:,3),connIdxRoot(:,3)').^2);

distMat = double(distMat==0);
if(sum(emptyIdxs))
  distMat(emptyIdxs,:) = NaN;
end

j = 1;
dendGroups = cell(1);
for i = 1:size(distMat,1)
  if(~isnan(distMat(i,i)))
    dendGroups{j} = find(distMat(i,:));
    distMat(dendGroups{j},:) = NaN;
    j = j+1;
  end
end

%%
offset = 2;
connIdx = zeros(length(dendGroups),3);
connRoots = zeros(length(dendGroups),3);
for i = 1:length(dendGroups)
  path = allpaths{dendGroups{i}(1)};
  pathInd = sub2ind(fdims,path(:,1),path(:,2),path(:,3));
  pathIntersect = ismember(pathInd,cellBody);
  try
    connIdx(i,:) = path(find(pathIntersect,1)-round(offset*sqrt(length(dendGroups{i}))),:);
  catch
    connIdx(i,:) = path(1,:);
  end
  connRoots(i,:) = path(find(pathIntersect,1),:);
end

%%
[xi,yi,zi] = ind2sub(fdims,cellBody);
cellInd = [xi yi zi];
clear xi yi zi
cellMin = min(cellInd);
cellMax = max(cellInd);

cellMat = false(fdims);
cellMat(cellBody) = true;

cellCrop = cellMat(cellMin(1):cellMax(1),cellMin(2):cellMax(2),cellMin(3):cellMax(3));

cellDiff = cellCrop(1:end-2,2:end-1,2:end-1)+cellCrop(3:end,2:end-1,2:end-1) + ...
  cellCrop(2:end-1,1:end-2,2:end-1)+cellCrop(2:end-1,3:end,2:end-1) + ...
  cellCrop(2:end-1,2:end-1,1:end-2)+cellCrop(2:end-1,2:end-1,3:end);

cellBorders = cellCrop;
cellBorders(2:end-1,2:end-1,2:end-1) = (cellDiff>0)&(cellDiff<6)&cellBorders(2:end-1,2:end-1,2:end-1);

cellBorders2 = false(fdims);
cellBorders2(cellMin(1):cellMax(1),cellMin(2):cellMax(2),cellMin(3):cellMax(3)) = cellBorders;

[xi,yi,zi] = ind2sub(fdims,find(cellBorders2));
bordersSub = [xi yi zi];
% bordersInd = find(cellBorders2(:));
clear xi yi zi

%%
% borderDist = sqrt(sum(bsxfun(@minus,connIdx2,bordersInd).^2,2));
cellProcessed = false(fdims);
testDist = [0 4 10];
numsamp = 20;
for j = 1:size(connRoots,1)
  distOff = min(max(testDist(2),round(offset*sqrt(length(dendGroups{j})))),testDist(3));
  borderDist = bsxfun(@minus,connRoots(j,:),bordersSub);
  % do ellipsoid scaling by doing a coordinate transform to the vector from
  % connIdx2 to connIdx2B
  borderDist = sqrt(sum(borderDist.^2,2));
  testIdx = find((borderDist<distOff).*(borderDist>testDist(1)));
  
  testSub = [];
  for i = 1:length(testIdx)
    spts = cscvn([connRoots(j,:)' connIdx(j,:)' bordersSub(testIdx(i),:)']);
    if(length(spts.breaks)>2)
      dpts = round(fnval(spts,linspace(spts.breaks(1),spts.breaks(3),numsamp))');
      testSub = cat(1,testSub,dpts);
    else
      dpts = round(fnval(spts,linspace(spts.breaks(1),spts.breaks(2),numsamp))');
      testSub = cat(1,testSub,dpts);      
    end
  end
  if(~isempty(testSub))
    testSub = bsxfun(@max,testSub,[1 1 1]);
    testSub = bsxfun(@min,testSub,fdims);
    testInd = sub2ind(fdims,testSub(:,1),testSub(:,2),testSub(:,3));
    
    cellBump = cellBorders2;
    cellBump(cellBody) = 1;
    cellBump(testInd) = 1;
    
    [xi,yi,zi] = ind2sub(fdims,find(cellBump));
    cellInd2 = [xi yi zi];
    clear xi yi zi
    cellMin2 = min(cellInd2);
    cellMax2 = max(cellInd2);
    
    numdiff = sum(cellBump(:));
    while(numdiff>0)
      cellCrop2 = cellBump(cellMin2(1):cellMax2(1),cellMin2(2):cellMax2(2),cellMin2(3):cellMax2(3));
      cellDiff2 = cellCrop2(1:end-2,2:end-1,2:end-1)+cellCrop2(3:end,2:end-1,2:end-1) + ...
        cellCrop2(2:end-1,1:end-2,2:end-1)+cellCrop2(2:end-1,3:end,2:end-1) + ...
        cellCrop2(2:end-1,2:end-1,1:end-2)+cellCrop2(2:end-1,2:end-1,3:end);
      cellCrop2(2:end-1,2:end-1,2:end-1) = (cellDiff2>=4)+cellCrop2(2:end-1,2:end-1,2:end-1);
      numdiff = sum(cellCrop2(:))-sum(cellBump(:));
      cellBump(cellMin2(1):cellMax2(1),cellMin2(2):cellMax2(2),cellMin2(3):cellMax2(3)) = cellCrop2;
    end
    cellProcessed = cellProcessed+cellBump;
  end
end

output = find(cellProcessed>0);