function write_TPM_movie(mov, save_name, varargin)
    
% function write_TPM_movie(mov,save_name,{dataType,dataHeader,blockSize})
%
% Write the 3-D array in mov to a file. Can writ to either a tiff stack, a
% fits file, or a matlab '.mat' file. The mode is automatically chosen
% based on the extension of the file-name passed to the function
% (save_name). The inputs to this function are
%     - mov       - 3D array of data to save (mov(:,:,kk) is the kk^th
%                   frame of the movie)
%     - save_name - String containing the file-name to save the data as
%
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if ~ischar(save_name)
    error('Need to provide a string as the file-name (second input)!')     % Check that the file-name is a string
end

if nargin > 2
    dataType = varargin{1};
else
    dataType = 'single';
end    

if nargin > 3
    dataHeader = varargin{2};
else
    dataHeader = 'empty';
end    

if nargin > 4
    blockSize = varargin{3};
else
    blockSize = 2500;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Data

switch save_name(end-3:end)                                                % Check the file type that the save-name string implies
    case 'fits'
        fprintf('Writing movie to fits...')
        fitswrite(mov, save_name);                                         % Write as a fits file OR
        fprintf('done.\n')
    case '.tif'
        fprintf('Writing movie to tif...\n')
%         tiff_writer(save_name,mov);                                        % Write as a tif file
        tifwriteblock(mov,save_name,dataHeader,dataType,blockSize);        % Write as a tif file
        fprintf('done.\n')
    case '.mat'
        fprintf('Writing movie to a mat file...\n')
        save(save_name,'mov','-v7.3')                                      % Write as a mat file
        fprintf('done.\n')
    otherwise
        error('Unknown file type. Can only save to .fits, .mat, or .tif file-types!')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
