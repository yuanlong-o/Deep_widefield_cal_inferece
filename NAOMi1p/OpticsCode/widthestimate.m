function [width,pl,pr] = widthestimate(vector,fraction)
%calculates fwhm if fraction is 0.5, fw 1/e for 1/e, etc., uses linear
%interpolation

if nargin<2
  fraction = 0.5;
end

vector = vector/max(vector(:));
greater = vector>fraction;
flip = find(diff(greater));
s1 = (vector(flip(1)+1)-fraction)/(vector(flip(1)+1)-vector(flip(1)));
s2 = (vector(flip(2))-fraction)/(vector(flip(2))-vector(flip(2)+1));
width = s1+s2+sum(greater);
if(nargout>1)
  center = find(vector==1);
  center = center(1);
  pl = flip(1)-s1;
  pr = flip(2)+s2;
end