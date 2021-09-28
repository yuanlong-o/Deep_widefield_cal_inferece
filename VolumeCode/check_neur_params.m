function neur_params = check_neur_params(neur_params)

% neur_params = check_neur_params(neur_params)
%  
% This function checks the elements of the struct neur_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - neur_params - Struct containing parameters for neuron generation
%       .n_samps     - Number of sphere samples to use in creating the mesh
%                      for generating soma and nucleus shapes (default =
%                      200) 
%       .l_scale     - length-scale for the isotropic GP of the soma
%                      shapes. This controls the shape `bumpiness' (default
%                      = 90) 
%       .p_scale     - Overall variance of the isotropic GP of the soma 
%                      shape. (default = 90) 
%       .avg_rad     - Average radius of each neuron in um (default =
%                      6.6 um) 
%       .nuc_fluorsc - Potential fluorescence in the nucleus (default =
%                      0)
%       .min_thic    - Minimum cytoplasmic thickness (default = 0.8)
%       .eccen       - Maximum eccentricity of neuron (default = 0.35)
%       .exts        - Parameters dictating the max/min of the soma radii
%                      (Default = [0.75,1.7])
%       .nexts       - Parameters dictating the extent to shrink and smooth
%                      the nucleus (Default = [0.5, 1])
%       .neur_type   - Option for neuron type (Default 'pyr')
%
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(neur_params)                                                    % Make sure that neur_params is a struct
    clear neur_params
    neur_params = struct;
end

dParams.n_samps      = 200;                                                % Default numbers of samples on the sphere is 200             
dParams.l_scale      = 90;                                                 % Default length-scale is 90     
dParams.p_scale      = 1000;                                               % Default correlation scaling is 1000
dParams.avg_rad      = 5.9;                                                % Default average radius is 5.9     
dParams.nuc_rad      = [5.65 2.5];                                         % Default nuclear radius    
dParams.max_ang      = 20;                                                 % Default maximum angle tilt    
dParams.plot_opt     = false;                                              % Default plotting setting is to plot (Default FALSE)
dParams.dendrite_tau = 50;                                                 % Default
dParams.nuc_fluorsc  = 0;                                                  % Nuclear fluoresence level
dParams.min_thic     = [0.4 0.4];                                          % Minimum cytoplasmic thickness 
dParams.eccen        = [0.35 0.35 0.5];                                    % Default maximum eccentricity of neuron is [0.35 0.35 0.5]
dParams.exts         = [0.75,1.7];                                         % Parameters dictating the max/min of the soma radii
dParams.nexts        = [0.5,1];                                            % Parameters dictating the extent to shrink and smooth a nucleus
dParams.neur_type    = 'pyr';                                              % Option for neuron type (Default 'pyr')
dParams.fluor_dist   = [1 0.2];                                            % Somatic neural fluoresence distribution (mean length [um], mean weight [0-1])


neur_params = setParams(dParams,neur_params);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
