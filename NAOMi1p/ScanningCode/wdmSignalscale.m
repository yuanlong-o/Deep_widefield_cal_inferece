function Ftavg = wdmSignalscale(wdm_params,psf_params, vol_params)

%  modified for 1p version.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nidx   = double(wdm_params.nidx);                                        % refractive index []
nac    = double(wdm_params.nac);                                         % collection NA []
lambda = double(wdm_params.lambda);    
phi   = double(wdm_params.phi);                                            % collection efficiency [], calcualted by objective model
eta   = double(wdm_params.eta);                                            % fluorophore QE []
conc  = double(wdm_params.conc);                                           % fluorophore concentration [uM]

epsilon = double(wdm_params.epsilon);                                          % one-photon exinction coefficient [M-1cm-1]

pavg  = double(wdm_params.pavg);                                           % average power [mW]

vol = 1e3 * 1e3; %% 1mm2, in micron^2
%% change unit
conc   = conc*1e-6; % molar concentration to number of molecules

lambda = lambda*1e-6;
pavg   = 1e-3*pavg/(6.626e-34 * 3e8/lambda); % times hc/lambda

Ftavg  = phi*eta*conc*epsilon* nidx* pavg * lambda / vol * vol_params.vres^2; % photons per micron^2
% add 2 multiplier to match the real situations
% about 1 photon per pulse (close to 80M photons/s) - or average 10 photons/pixel
