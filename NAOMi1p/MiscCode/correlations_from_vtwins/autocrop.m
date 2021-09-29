function [im2,bb] = autocrop(im,bg)
%autocrops zeros from image

if nargin < 2
  bg = 0;
end
try
  imbin = (im>bg);
  box = regionprops(imbin,'BoundingBox');
  bb = round(box(1).BoundingBox);
  bb(3) = bb(1)+bb(3);
  bb(4) = bb(2)+bb(4);
  if length(box)>=2
    for i = 2:length(box)
      bounds = round(box(i).BoundingBox);
      bb(1) = min(bb(1),bounds(1));
      bb(2) = min(bb(2),bounds(2));
      bb(3) = max(bb(3),bounds(1)+bounds(3));
      bb(4) = max(bb(4),bounds(2)+bounds(4));
    end
  end
  bb(1) = max(1,bb(1));
  bb(2) = max(1,bb(2));
  bb(3) = min(bb(3),size(im,2));
  bb(4) = min(bb(4),size(im,1));
  im2 = im(bb(2):bb(4),bb(1):bb(3));
catch
  im2 = im;
end