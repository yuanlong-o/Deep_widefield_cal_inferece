function frame_out = pixel_bleed(frame, p, b_max)

% function frame = pixel_bleed(frame, p, b_max)
%
% This function calculates the pixel-pixel bleed-through from the
% electronics setup. It samples the amount of the pulses from one
% integration/sampling period that leak into the next integration/sampling
% period. The inputs are:
%   - frame - Input video frame (fluorescence returns across the
%             field-of-view)
%   - p     - Probability of bleed-through occuring
%   - b_max - Maximum bleed-through 
%
% The output is
%   - frame_out - Output frame with pixel bleed-through
%
% 2017 - Adam Charles and Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(p>0)
x_bleed = b_max*max(rand(size(frame))-(1-p),0)/p;                          % Bleedthrough hapens with probability p and a max bleed-through of b_max
frame_out = frame - x_bleed.*frame + ...
                    [[0;x_bleed(1:(end-1),end)],x_bleed(:,1:(end-1))].*...
                            [[0;frame(1:(end-1),end)],frame(:,1:(end-1))]; % Actual amount of bleedthrough is a subtracted, shifted and scaled copy of the values
else
  frame_out = frame;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%