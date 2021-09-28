function reSampleMovie(fullPath, varargin)


if nargin > 1;    boostLvl = varargin{1};
else;             boostLvl = [];
end

fprintf('Loading movie parts...')
[save_path,save_filename,~] = fileparts(fullPath);                         % Get parts of filename
m = matfile(fullPath,'Writable',true);                                     % Access the mat file with writable permission

if isempty(boostLvl)
    fprintf('No boost level provided. Auto-normalizing to ~1 photon/pixel.')
    vol_params = m.vol_params;
    boostLvl   = prod(vol_params.vol_sz(1:2))/(250)^2;  
    fprintf('done.\n')
end

if boostLvl ~= 1
    noise_params = m.noise_params;                                             % Load the parameter struct that needs to be adjusted
    Fsim_clean   = m.Fsim_clean;                                               % Load the clean video
end

fprintf('done.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Actually adjust the power

if boostLvl ~= 1
    fprintf('Boosting clean movie and updating signal scale parameter...')
    noise_params.sigscale = boostLvl*noise_params.sigscale;                    % Boost the signal scale parameter
    Fsim_clean            = boostLvl*Fsim_clean;                               % Boost the actual signal
    fprintf('done.\n')
    fprintf('Resampling from the noise model with new signal scale...')
    Fsim                  = applyNoiseModel(Fsim_clean, noise_params);         % Rerun the noise model with the boosted singal
    fprintf('done.\n')

    fprintf('Saving new parameters and boosted movies...')
    m.Fsim_clean          = Fsim_clean;                                        % Write back to the mat file the adjusted clean video 
    m.Fsim                = Fsim;                                              % Write back to the mat file the adjusted sampled video 
    m.noise_params        = noise_params;                                      % Write back the parameter struct that needs to be adjusted
    fprintf('done.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clean the directory of older tif and avi files

tifFiles = strcat(save_path,'/',allFilesOfType(save_path, '.tif'));        % Find all the TIF files to delete for replacement
aviFiles = strcat(save_path,'/',allFilesOfType(save_path, '.avi'));        % Find all the AVI files to delete for replacement
if ~isempty(tifFiles); delete(tifFiles{:}); end                            % Delete all the TIF files
if ~isempty(aviFiles); delete(aviFiles{:}); end                            % Delete all the AVI files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make a movie of the results
% The output can also be saved as an avi movie for display purposes. As an
% example, the following line of code saves the video to such a file.

fprintf('Creating movies...')
if boostLvl ~= 1
    make_avi(Fsim,       [save_path,'/',save_filename, '.avi'],      0.2); % Make an avi of the noisy video
    make_avi(Fsim_clean, [save_path,'/',save_filename, 'clean.avi'], 0.2); % Make an avi of the clean video
else
    make_avi(m.Fsim,       [save_path,'/',save_filename, '.avi'],     0.2);% Make an avi of the noisy video
    make_avi(m.Fsim_clean, [save_path,'/',save_filename, 'clean.avi'],0.2);% Make an avi of the clean video
end
fprintf('done.\n')

fprintf('Saving as TIF stacks...')
if boostLvl ~= 1
    write_TPM_movie(Fsim,       [save_path,'/',save_filename, '_mov.tif'])
    write_TPM_movie(Fsim_clean, [save_path,'/',save_filename, '_movClean.tif'])
    clear Fsim Fsim_clean
else
    write_TPM_movie(m.Fsim,       [save_path,'/',save_filename, '_mov.tif'])
    write_TPM_movie(m.Fsim_clean, [save_path,'/',save_filename, '_movClean.tif'])
end
fprintf('done.\n')


fprintf('Creating separate model parts...')
saveSimulationParts([save_path,'/',save_filename,'.mat'])
fprintf('done.\n')
fprintf('Finished!\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end