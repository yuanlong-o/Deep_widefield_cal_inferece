function tifwrite(imData,imFilename,imHeader,dataType)
% writes a tif file with data, filename, header, and type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin<4
  dataType = class(imData);
else
  imData = eval([dataType '(imData)']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(size(imData)) > 3 || length(size(imData)) < 2
  error('image is not 2D or 3D')
end

tifLink = Tiff(imFilename,'w');
if nargin >= 3
  tagStruct.ImageDescription = imHeader;
end

tagStruct.ImageLength         = size(imData,1);
tagStruct.ImageWidth          = size(imData,2);
tagStruct.Photometric         = Tiff.Photometric.MinIsBlack;
tagStruct.SamplesPerPixel     = 1;
tagStruct.RowsPerStrip        = 8;
tagStruct.Compression         = Tiff.Compression.None;
tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagStruct.Software            = 'MATLAB';


switch dataType
  case 'logical'
    tagStruct.BitsPerSample = 8;
    tagStruct.SampleFormat  = Tiff.SampleFormat.UInt;
    imData                  = uint8(imData);
  case 'char'
    tagStruct.BitsPerSample = 8;
    tagStruct.SampleFormat  = Tiff.SampleFormat.UInt;
    imData                  = uint8(imData);
  case 'uint8'
    tagStruct.BitsPerSample = 8;
    tagStruct.SampleFormat  = Tiff.SampleFormat.UInt;
  case 'int8'
    tagStruct.BitsPerSample = 8;
    tagStruct.SampleFormat  = Tiff.SampleFormat.Int;
  case 'uint16'
    tagStruct.BitsPerSample = 16;
    tagStruct.SampleFormat  = Tiff.SampleFormat.UInt;
  case 'int16'
    tagStruct.BitsPerSample = 16;
    tagStruct.SampleFormat  = Tiff.SampleFormat.Int;
  case 'uint32'
    tagStruct.BitsPerSample = 32;
    tagStruct.SampleFormat  = Tiff.SampleFormat.UInt;
  case 'int32'
    tagStruct.BitsPerSample = 32;
    tagStruct.SampleFormat  = Tiff.SampleFormat.Int;
  case 'single'
    tagStruct.BitsPerSample = 32;
    tagStruct.SampleFormat  = Tiff.SampleFormat.IEEEFP;
  case 'double'
    tagStruct.BitsPerSample = 64;
    tagStruct.SampleFormat  = Tiff.SampleFormat.IEEEFP;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% White to file

tifLink.setTag(tagStruct);
if length(size(imData))==2
    tifLink.setTag(tagStruct)
    tifLink.write(imData);
else
    tifLink.setTag(tagStruct)
    tifLink.write(imData(:,:,1));
    for i = 2:size(imData,3)
        tifLink.writeDirectory();
        tifLink.setTag(tagStruct)
        tifLink.write(imData(:,:,i));
    end
end
tifLink.close();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
