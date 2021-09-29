function width = widthestimate3D(matrix,fraction)
%calculates fwhm if fraction is 0.5, fw 1/e for 1/e, etc., uses linear
%interpolation

if nargin<2
  fraction = 0.5;
end

vec1 = squeeze(sum(sum(matrix,3),2));
vec2 = squeeze(sum(sum(matrix,3),1));
vec3 = squeeze(sum(sum(matrix,2),1));

width1 = widthestimate(vec1,fraction);
width2 = widthestimate(vec2,fraction);
width3 = widthestimate(vec3,fraction);
width = [width1 width2 width3];