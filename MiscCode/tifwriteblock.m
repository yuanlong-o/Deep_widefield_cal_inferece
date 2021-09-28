function tifwriteblock(imData,imFilename,imHeader,dataType,blocksize)


if length(size(imData)) ~= 3
  error('image is not 2D or 3D')
end

if nargin<5 || isempty(blocksize)
  blocksize = 1000;
end

if nargin<4 || isempty(dataType)
  dataType = class(imData);
else
  imData   = eval([dataType '(imData)']);
end

if nargin<3
  imHeader = 'empty';
end

[path,name,ext] = fileparts(imFilename);
for i = 1:ceil(size(imData,3)/blocksize)
    fprintf('writing file %d of %d...\n', i, ceil(size(imData,3)/blocksize))
    TMPname = fullfile(path,sprintf([name '_%05d' ext],i));
    if(blocksize+((i-1)*blocksize)<=size(imData,3))
        tifwrite(imData(:,:,(1:blocksize)+((i-1)*blocksize)),TMPname,...
                                                     imHeader,dataType);
    else
        tifwrite(imData(:,:,(1+((i-1)*blocksize)):end),TMPname,...
                                                     imHeader,dataType);        
    end
end
fprintf('finished.\n')
