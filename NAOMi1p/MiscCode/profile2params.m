function [A_params,paramsM] = profile2params(A,z_tilt)

% profile2params(A)
%
% Function to extract the position parameters from spatial profiles. The
% inputs to this function are
%   A:        The 3D matrix of spatial profiles
%   z-tilt:   An optional input indicating the conversion ratio of px/um
%             for the vTwINS PSF, calculated as 2*tand(theta)*ss, where theta
%             is the angle of one half of the PSF to normal, and ss is the
%             scaling factor from px to um in the image
% The outputs of this function are
%   A_params: A 4 column vector corresponding to X-position (px),
%             Y-position (px), pair separation (px), and depth (um),
%             estimated using the centroid of each half of the spatial
%             profile.
%   paramsM:  Same as above except the center pixel position of each half
%             of the spatial profile is used instead.
% 2016 - Alexander Song

if nargin < 2
  z_tilt = 1;
end


A_params = nan(size(A,3),4);
paramsM = nan(size(A,3),4);
for i = 1:size(A,3)
  try
    im = A(:,:,i);
    imbin = (im>0);
    box = regionprops(imbin,'BoundingBox');
    bb = round(box(1).BoundingBox);
    bb(3) = bb(1)+bb(3);
    bb(4) = bb(2)+bb(4);
    if length(box)>=2
      for j = 2:length(box)
        bounds = round(box(j).BoundingBox);
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
    im = im(bb(2):bb(4),bb(1):bb(3));
    
    imL = im(:,1:ceil(size(im,2)/2));
    imR = im(:,floor((size(im,2)+1)/2):end);
    
    propsL = regionprops(imL>0);
    propsR = regionprops(imR>0);
    
    propsL = propsL(cat(1,propsL.Area)==max(cat(1,propsL.Area)));
    propsR = propsR(cat(1,propsR.Area)==max(cat(1,propsR.Area)));
    
    propsL = propsL(1);
    propsR = propsR(1);
    
    
    cenL = [bb(1) bb(2)]+propsL.Centroid;
    cenR = [bb(1) bb(2)]+propsR.Centroid+[size(im,2)-size(imR,2) 0];
    
    pLbb = propsL.BoundingBox;
    pRbb = propsR.BoundingBox;
    
    midL = [bb(1) bb(2)]+[pLbb(1)+pLbb(3)/2 pLbb(2)+pLbb(4)/2];
    midR = [bb(1) bb(2)]+[pRbb(1)+pRbb(3)/2 pRbb(2)+pRbb(4)/2]+[size(im,2)-size(imR,2) 0];
    
    A_params(i,:) = [(cenL(1)+cenR(1))/2 (cenL(2)+cenR(2))/2 cenR(1)-cenL(1) (cenR(1)-cenL(1))/z_tilt];
    paramsM(i,:) = [(midL(1)+midR(1))/2 (midL(2)+midR(2))/2 midR(1)-midL(1) (midR(1)-midL(1))/z_tilt];
  end
end