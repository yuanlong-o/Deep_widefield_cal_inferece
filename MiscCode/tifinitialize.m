function [tiflink,tagStruct] = tifinitialize(filename,imsize)

% [tiflink,tagStruct] = tifinitialize(filename,imsize)
%
% Function to write initialize saving appended tif to disk
%
% 2019 - Alex Song

tagStruct.ImageLength = imsize(1);
tagStruct.ImageWidth = imsize(2);
tagStruct.Photometric = Tiff.Photometric.MinIsBlack;
tagStruct.SamplesPerPixel = 1;
tagStruct.RowsPerStrip = 8;
tagStruct.Compression = Tiff.Compression.None;
tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagStruct.Software = 'MATLAB';
tagStruct.BitsPerSample = 32;
tagStruct.SampleFormat = Tiff.SampleFormat.IEEEFP;

[path,name,~] = fileparts(filename);
if(~isempty(path))
  mkdir(path);
end
filename2 = fullfile(path,sprintf([name '_%05d.tif'],1));
tiflink = Tiff(filename2,'w');
tiflink.setTag(tagStruct);
