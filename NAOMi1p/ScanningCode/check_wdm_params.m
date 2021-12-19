function wdm_params = check_wdm_params(wdm_params)

% widefield imaging parameters
% last update: 9/18/2021. YZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(wdm_params)                                                     % Make sure that tpm_params is a struct
    clear wdm_params
    wdm_params = struct;
end

dParams.nidx   = 1;                                                        % 1
dParams.nac    = 0.3;                                                      % Default objective NA
dParams.phi    = [];                                                       % if this is empty, it will default to 80 percent transmission and 0.8NA light all on detector (no strong scattering) with average 40 perfect PMT QE
dParams.eta    = 0.6;                                                      % eGFP quantum yield

dParams.conc   = 10;                                                       % extinction value
dParams.epsilon  = 5.6e5;                                                  
dParams.qe = 0.7; % camera QE
dParams.lambda = 0.48;                                                     % 0.48um is the default for excitation wavelength

wdm_params = setParams(dParams, wdm_params);

% calculate collection efficiency of the system
if isempty(wdm_params.phi)
    wdm_params.phi = 0.8*((1-sqrt(1-(wdm_params.nac/wdm_params.nidx)^2))/2)*0.4;   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
