function img_out = imgSubRowShift(img_in, buf_sz, x_off, y_off)

% img_out = imgSubRowShift(img_in, buf_sz, x_off, y_off)
%
% Subselect part of an image with per-row shifts. This sub-function was
% cluttering up the min function and was stand-alone in it's functionality.
%
% Inputs
%   - img_in      Input image to be sampled from
%   - buf_sz      Buffer on the edges of the image to avoid
%   - x_off       x-dimension offsets (per row)
%   - y_off       y-dimension offsets (per row)
%
% Outputs
%   - img_out     Output image sampled from img_in that consists of the 
%                 off-set rows starting from x_off down, with offests 
%                 give by y_off. 
%
% 2016 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate offsets

if(length(x_off)==1)
  x_off = x_off*ones(size(y_off));
end

img_tmp = zeros(size(img_in),'single');
img_out = zeros(size(img_in),'single');
imsz = size(img_in);

x_off = x_off-buf_sz;
y_off = y_off-buf_sz;
x_off2 = x_off+(0:length(x_off)-1)';

for kk = 1:size(x_off)
  if(mod(x_off2(kk),1)>0)
    if(floor(x_off2(kk))<=imsz(1) && floor(x_off2(kk))>=1)
      TMP1 = img_in(floor(x_off2(kk)),:);
    else
      TMP1 = nan(1,imsz(2),'single');
    end
    if(ceil(x_off2(kk))<=imsz(1) && ceil(x_off2(kk))>=1)
      TMP2 = img_in(ceil(x_off2(kk)),:);
    else
      TMP2 = nan(1,imsz(2),'single');
    end
    img_tmp(kk,:) = (TMP1*(1-mod(x_off2(kk),1)))+(TMP2*mod(x_off2(kk),1));
  else
    if(x_off2(kk)<=imsz(1) && x_off2(kk)>=1)
      img_tmp(kk,:) = img_in(x_off2(kk),:);
    else
      img_tmp(kk,:) = nan(1,imsz(2),'single');
    end
  end
end
offset = ceil(max(abs(y_off)));
img_tmp = padarray(img_tmp,[0 offset],nan,'both');

for kk = 1:size(img_out,1)  
  TMP1 = img_tmp(kk,(floor(y_off(kk))+offset)+(1:imsz(2)));
  TMP2 = img_tmp(kk,(ceil(y_off(kk))+offset)+(1:imsz(2)));
  img_out(kk,:) = (TMP1*(1-mod(y_off(kk),1)))+(TMP2*mod(y_off(kk),1));
end

img_out = img_out(buf_sz+1:end-buf_sz,buf_sz+1:end-buf_sz);

% img_out = zeros(size(img_in,1)-2*buf_sz,size(img_in,2)-2*buf_sz,'single'); % Initialize output image
%
% for kk = 1:size(img_out,1)
%   if(x_loc+kk-1>=1 && x_loc+kk-1<=size(img_out,1))
%     if(y_off(kk)<1)
%       img_out(kk, 2-y_off(kk):end) = img_in(x_loc+kk-1,1:(y_off(kk)+size(img_out,2)-1));        % Extract k^th row      
%     elseif((y_off(kk)+size(img_out,2)-1)>size(img_in,2))
%       img_out(kk, 1:end-y_off(kk)+1) = img_in(x_loc+kk-1,y_off(kk):end);        % Extract k^th row
%     else
%       img_out(kk, :) = img_in(x_loc+kk-1,y_off(kk):(y_off(kk)+size(img_out,2)-1));        % Extract k^th row      
%     end
%   end
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
