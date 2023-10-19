function  avi_to_tiff(input_path, output_path, ROI)
%AVI_TO_TIFF 此处显示有关此函数的摘要
%   此处显示详细说明
    vidObj = VideoReader(input_path);
    ROI = [350 450 1900 2050];
    width = ROI(4) - ROI(2) + 1;
    height = ROI(3) - ROI(1) + 1;
    frameNumber = vidObj.NumFrames;
    stack = zeros(height, width, frameNumber, 'uint8');
    i = 1;
       
    fprintf('loading images... Frame count: %d \n', frameNumber);
    while hasFrame(vidObj)
        vidFrame = readFrame(vidObj);
        frame = vidFrame(ROI(1):ROI(3), ROI(2):ROI(4));
        %frame = vidFrame;
        stack(:,:,i) = frame;
        i = i + 1;
        if mod(i,100) == 0
            disp(i);
        end
    
    end
   
    disp('saving images...');
    clear options;
    options.big = true;
    saveastiff(stack, output_path,options);
    fprintf('done. image saved to %s \n', output_path);

end

