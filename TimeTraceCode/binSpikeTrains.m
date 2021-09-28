function S = binSpikeTrains(evt, evm, N_node, dt, T)

% function S = binSpikeTrains(evt, evm, N_node)
% 
% Function that takes a list of marked events (times in evt and marks that
% are integers in [1, N_node] in evm) and makes a spike count matrix S
% where each element is the number of spikes for each neuron (the rows) and
% a time-bin of size dt for each of the T times (the columns).
% Inputs:
%  - evt    - List of event times
%  - evm    - List of event marks (neuron numbers)
%  - N_node - Total number of possible marks
%  - dt     - Bin size (time-span) of each time-pooint
%  - T      - Total number of time-points
%
% Outputs:
%  - S      - N_node-by-T matrix of binned spike counts
%
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if numel(evt)~=numel(evm)
    error('Number of events must be the same as the number of marks!')
end

if max(evm) > N_node
    error('Total number of marks must be greater than the highest seen mark!')
end

if ceil(max(evt)/dt) > T
    error('Total time-span must be longer than the latest seen event!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get bin counts

S = zeros(N_node, T);                                                      % Initialize count matrix
for kk = 1:numel(evt)                                                      % Iterate over events
    S(evm(kk), ceil(evt(kk)/dt)) = S(evm(kk), ceil(evt(kk)/dt)) + 1;       % Add a spike to the location/mark of the event
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%