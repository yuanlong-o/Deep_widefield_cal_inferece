function [comps,baseim,ideal] = calculateIdealComps(vol_out,PSF_struct,neur_act,scan_params,noise_params, spike_opts, tpm_params, num_comps)

% function [comps,baseim,ideal] = calculateIdealComps(vol_out,PSF_struct,neur_act,scan_params,noise_params, spike_opts, tpm_params, num_comps)
%
% Function to calculate the projection of each cell 
%
% 2018 - Alex Song & Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scan for components


conds = fieldnames(neur_act);                                              % 
if(isstruct(neur_act.(conds{1})))
    neur_base.soma = ones(size(neur_act.(conds{1}).soma,1));               % 
    neur_base.dend = ones(size(neur_act.(conds{1}).dend,1));               % 
    neur_base.bg = ones(size(neur_act.(conds{1}).bg,1));                   % 
    scan_params.motion = 0;                                                % No motion in getting the scanning fpr the componentes 

    if nargin>5
        compIdxs       = 1:num_comps;                                      % Extract component indecies
        neur_base.soma = neur_base.soma(:,compIdxs);                       % Pull out the corresponding soma time-traces...
        neur_base.dend = neur_base.dend(:,compIdxs);                       %   ... and dendrite time-traces...
        neur_base.bg   = neur_base.bg(:,compIdxs);                         %   ... and background time-traces
    end
  
    [~,TMP] = scan_volume(vol_out, PSF_struct, neur_base, scan_params, ...
                                    noise_params, spike_opts, tpm_params); % Perform the scanning simulation
  
    for i = 1:length(conds)
        weights          = min(neur_act.(conds{i}).soma,[],2);             % -----
        comps.(conds{i}) = bsxfun(@times,TMP,reshape(weights,1,1,[]));     % 
    end
    clear TMP                                                              % Clear temporary variables to reduce RAM
  
    % Scan for base image
    if nargout > 1
        for i = 1:length(conds)
            f0.soma       = min(neur_act.(conds{i}).soma,[],2);            % -----
            f0.dend       = min(neur_act.(conds{i}).dend,[],2);            %   | - Create temporary "activity" with one time point - the minimum values
            f0.bg         = min(neur_act.(conds{i}).bg,[],2);              % -----
            [~,TMPbaseim] = scan_volume(vol_out, PSF_struct, f0, ...
                       scan_params, noise_params, spike_opts, tpm_params); % Perform the scanning simulation
            baseim.(conds{i}) = TMPbaseim;
        end
    end
  
    % Calculate SNR-adjusted ideal components
    if nargout > 2
        for i = 1:length(conds)
            ideal.(conds{i}) = comps2ideals(comps.(conds{i}), ...
                                                       baseim.(conds{i})); % Calculate the SNR-adjusted ideal components from the full components
        end
    end
  
else
    % Initializations
    n0 = min(neur_act.soma,[],2);                                          % Get minimum value for each neuron's soma activity
    n1 = min(neur_act.dend,[],2);                                          % Get minimum value for each neuron's dendrite activity
    n2 = min(neur_act.bg,[],2);                                            % Get minimum value for each background component's activity
  
    gIdxs          = (n0+n1+n2)>0;
    neur_base.soma = zeros(size(neur_act.soma,1),sum(gIdxs),'single');
    neur_base.dend = zeros(size(neur_act.dend,1),sum(gIdxs),'single');
    neur_base.bg   = zeros(size(neur_act.bg,1),sum(gIdxs),'single');

    neur_base.soma(sub2ind(size(neur_base.soma),find(gIdxs)',1:sum(gIdxs))) = n0(gIdxs);
    neur_base.dend(sub2ind(size(neur_base.dend),find(gIdxs)',1:sum(gIdxs))) = n1(gIdxs);
    neur_base.bg(sub2ind(size(neur_base.bg),find(gIdxs)',1:sum(gIdxs))) = n2(gIdxs);

    scan_params.motion = 0;                                              % No motion in getting the scanning fpr the componentes 

    if nargin>5
        compIdxs = 1:sum(gIdxs(1:num_comps));                            % Extract component indecies
    else
        compIdxs = 1:sum(gIdxs);
    end
    neur_base.soma = neur_base.soma(:,compIdxs);                         % Pull out the corresponding soma time-traces...
    neur_base.dend = neur_base.dend(:,compIdxs);                         %   ... and dendrite time-traces...
    neur_base.bg   = neur_base.bg(:,compIdxs);                           %   ... and background time-traces
  
    [~,TMP] = scan_volume(vol_out, PSF_struct, neur_base, scan_params,...
                                  noise_params, spike_opts, tpm_params); % Perform the scanning simulation
  
    if nargin>5
        comps = zeros(size(TMP,1),size(TMP,2),num_comps);
        comps(:,:,gIdxs(1:num_comps)) = TMP;
    else
        comps = zeros(size(TMP,1),size(TMP,2),length(gIdxs));
        comps(:,:,gIdxs) = TMP;
    end
    clear TMP;
  
    % Scan for base image
    if nargout > 1
        f0.soma = n0;                                                    % -----
        f0.dend = n1;                                                    %   | - Create temporary "activity" with one time point - the minimum values
        f0.bg   = n2;                                                    % -----
        [~,baseim] = scan_volume(vol_out, PSF_struct, f0, scan_params,...
                                  noise_params, spike_opts, tpm_params); % Perform the scanning simulation
    end
  
    % Calculate SNR-adjusted ideal components
    if nargout > 2
        ideal = comps2ideals(comps, baseim);                             % Calculate the SNR-adjusted ideal components from the full components
    end
  
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
