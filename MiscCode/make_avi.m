function make_avi(Fmov, video_name, varargin)

% make_avi(video_name, Fmov, contrast_val)
%
% Function to write the movie Fmov to an avi movie file. The inputs to this
% function are:
%   - Fmov         - 3D array of data to save (mov(:,:,kk) is the kk^th
%                    frame of the movie)
%   - video_name   - String containing the file-name to save the data as
%   - contrast_val - OPTIONAL input that indicates the contrast with which
%                    to plot each frame
%
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if ~ischar(video_name)
    write_to_file = false;                                                 % Check that the file-name is a string - if not do not write to a file
else
    write_to_file = true;
    if ~strcmp(video_name(end-3:end),'.avi')
        error('File name should have a .avi extension!')                   % Make sure to write to a .avi file
    end
end

if nargin > 2
    contrast_val = varargin{1};                                            % Check if the contrast value is provided
else
    contrast_val = 0.9;                                                    % Default contrast value is 0.9
end

if nargin > 3
    avg_val = varargin{2};                                                 % Check if the averaging frane size is provided
else
    avg_val = 5;                                                           % Default 5-frame average
end

if nargin > 4
    dff = varargin{3};                                                     % Check if the video is f or delta f/f 
else
    dff = false;                                                           % Default to video of f
end

if isempty(avg_val)
   avg_val = 5;
end

if isempty(dff)
   dff = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write video

if write_to_file
    writerObj = VideoWriter(video_name);                                   % Set up the video object to work with
    open(writerObj)                                                        % Open movie object 
end
if dff
    Fmed = median(Fmov,3);
    maxF = max(vec(bsxfun(@times,bsxfun(@plus,Fmov,-Fmed),1./Fmed)));      % Get the maximum fluorescence value
else
    maxF = max(Fmov(:));                                                   % Get the maximum fluorescence value
end

clims = [0*min(Fmov(:)), contrast_val*maxF];                               % Set constant color limits for the video frames

for kk = 1:(size(Fmov,3)-avg_val)
    if dff
        imagesc((mean(Fmov(:,:,kk:(kk+avg_val-1)),3) - Fmed)./Fmed,clims)  % Make the image of the kk^th frame
    else
        imagesc(mean(Fmov(:,:,kk:(kk+avg_val-1)),3),clims)                 % Make the image of the kk^th frame
    end
    colormap gray                                                          % Colorscale should be gray
    axis image                                                             % Make sure the axis sizes are reasonable  
    axis off                                                               % Remove axis numbering
    title(sprintf('Time: %3.3f', (kk-1)/30),'FontSize',20)                 % Set the title
    set(gcf,'color',[1,1,1])
    drawnow
    if write_to_file
        writeVideo(writerObj,getframe(gcf));                               % Write movie frame
    end
    pause(0.01)
end
if write_to_file
    close(writerObj);                                                      % Close movie object
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
