function saveSimulationParts(fullFileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get basic data and load the .mat file

fprintf('Loading mat file...')
fileBase = fullFileName(1:end-4);
m        = matfile(fullFileName);
fprintf('done.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save parameters

fprintf('Saving parameter set...')
axon_params  = m.axon_params;
bg_params    = m.bg_params;
dend_params  = m.dend_params;
neur_params  = m.neur_params;
noise_params = m.noise_params;
psf_params   = m.psf_params;
scan_params  = m.scan_params;
spike_opts   = m.spike_opts;
tpm_params   = m.tpm_params;
vasc_params  = m.vasc_params;
vol_params   = m.vol_params;
save([fileBase,'_allParameters.mat'], 'axon_params', 'bg_params', ...
          'dend_params', 'neur_params', 'noise_params', 'psf_params',...
          'scan_params', 'spike_opts', 'tpm_params', 'vasc_params', ...
                                                    'vol_params','-v7.3'); % Save all the parameters
fprintf('done.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save volume

fprintf('Saving volume...')
vol_out = m.vol_out;                                                       % Load the volume
save([fileBase,'_volume.mat'], 'vol_out', 'axon_params', 'bg_params', ...
       'dend_params', 'neur_params', 'vasc_params', 'vol_params','-v7.3'); % Save the volume and associated parameters

clear vol_out                                                              % Keep the workplace clean
fprintf('done.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save PSF

fprintf('Saving psf...')
PSF_struct = m.PSF_struct;
save([fileBase,'_psf.mat'], 'PSF_struct', 'noise_params', ...
         'psf_params', 'scan_params', 'tpm_params', 'vol_params','-v7.3'); % Save the PSF and related parameters

clear PSF_struct                                                           % Keep the workplace clean
fprintf('done.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Activity

fprintf('Saving activity...')
neur_act = m.neur_act;                                                     % Load the neural activity
spikes   = m.spikes;                                                       % Load the underlying spiking data
save([fileBase,'_activity.mat'], 'spikes', 'neur_act', 'spike_opts', ...
                                                                 '-v7.3'); % Save all the parameters

clear spikes neur_act                                                      % Keep the workplace clean
fprintf('done.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save ideal profiles 

fprintf('Saving ideal components...')
ideal       = m.ideal;
comps       = m.comps;
idealTraces = m.idealTraces;
save([fileBase,'_idealGT.mat'], 'comps', 'ideal', 'idealTraces', '-v7.3'); % Save all the parameters

clear spikes neur_act                                                      % Keep the workplace clean
fprintf('done.\n')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
