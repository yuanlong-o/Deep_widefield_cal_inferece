function [finalImage, infoImage] = tifread(fileTif,subImages)
%Tifread - FileTif is a string specifying file location,
%subimages is len2 vector specifying start and end images (default all)

infoImage=imfinfo(fileTif);
mImage=infoImage(1).Width;
nImage=infoImage(1).Height;

if nargin < 2
  subImages = [1 length(infoImage)];
end

tifLink = Tiff(fileTif, 'r');
tifLink.setDirectory(subImages(1));
im = tifLink.read();
finalImage=zeros(nImage,mImage,subImages(2)-subImages(1)+1,'like',im);
finalImage(:,:,1)=im;
for i=1:subImages(2)-subImages(1)
  tifLink.setDirectory(i+subImages(1));
  finalImage(:,:,i+1)=tifLink.read();
end
tifLink.close();
