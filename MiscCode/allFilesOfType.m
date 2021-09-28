function filename = allFilesOfType(dataPath, varargin)

% Function to extract all .mat filenames in the fiven folder
%
% 2020 Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing
if nargin > 1;     fileExt = varargin{1};
else         ;     fileExt = '.mat';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check all the files in the directory

listing   = dir(dataPath);                                                 % Get list of files in the directory
filename  = cell(numel(listing),1);                                        % Initialize the cell array of files to segment
num_files = 0;                                                             % Initialize a counter for the number of files to segment
for kk = 1:numel(filename)
    TMP_name = listing(kk).name;                                           % Filename to be processed
    if (~listing(kk).isdir)&&(strcmp(TMP_name(end-3:end),fileExt))          % If not a directory and ends in .mat, add it to the list
        num_files           = num_files + 1;                               % Increase the number of found mat files
        filename{num_files} = TMP_name;                                    % Save the filename
    end
end
filename = filename(1:num_files);                                          % Store only the relevant file-names
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%