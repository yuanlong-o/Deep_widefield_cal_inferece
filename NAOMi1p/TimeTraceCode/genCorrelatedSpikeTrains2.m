function [S, varargout] = genCorrelatedSpikeTrains2(spike_opts, varargin)

% function S = genCorrelatedSpikeTrains(spike_opts) 
% 
% Function to generage a set of correlated spike trains by using a Hawkes
% process, or a discrete approximation of a Hawkes process. 
% 
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

spike_opts = check_spike_opts(spike_opts);                                 % Check the spiking option parameter struct
tmax       = spike_opts.dt*spike_opts.nt;                                  % Get the max time for the spiking simulation
N_node     = spike_opts.K;                                                 % Get the number of neurons
N_bg       = spike_opts.N_bg;                                              % Get the number of background processes

if nargin > 1
    batch_sz = varargin{1};
    if isempty(batch_sz)
        batch_sz = N_node;
    end
else
    batch_sz = N_node;                                                     % Default batch size is N_node (all nodes simulated together)
end

if nargin > 2
    n_locs = varargin{2};
else
    n_locs = [];                                                           % Default to no text spatial information for the traces
end

if nargin > 3
    discrete_flag = varargin{3};
else
    discrete_flag = true;                                                  % Default to running the faster, less exact, discrete approximation of teh hawkes process
end

if nargin > 4
    verbose = varargin{4};
else
    verbose = false;                                                        % Default to no text status updates
end

if (~isfield(spike_opts,'selfact'))||isempty(spike_opts.selfact)   
    selfAct = 1.26;                                                        % A scaling activity
else
    selfAct = spike_opts.selfact;
end


batch_num = ceil(N_node/batch_sz);                                         % Calculate number of batches
batch_bg  = ceil(N_bg/batch_num);                                          % Calculate number of background components per batch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample connectivity for Hawkes process
% spike_opts.rate = 0.2;
ascale = 4;
bscale = 2;
if verbose
    fprintf('Generating small world connectivity...')
end
MU = cell(batch_num,1);                                                    % Initialize base rate cell array
A  = cell(batch_num,1);                                                    % Initialize connectivity cell array
B  = cell(batch_num,1);                                                    % Initialize self-inhib cell array
for kk = 1:batch_num
    N_now  = min(batch_sz, N_node-(kk-1)*batch_sz);                        % Number of neurons to include in this batch
    N_nowB = min(batch_bg, N_bg-(kk-1)*batch_bg);                          % Number of background componants to include in this batch
    if discrete_flag
        A{kk}  = sampSmallWorldMat([N_now,N_nowB], 10, 0.3, 0.9, ...
                                           spike_opts.burst_mean, n_locs); % Generate random small-world matrix of excitation for the Hawkes process
    elseif ~discrete_flag
        A{kk}  = sampSmallWorldMat([N_now,N_nowB], 10, 0.3, 0.9, ...
                                           spike_opts.burst_mean, n_locs); % Generate random small-world matrix of excitation for the Hawkes process
    end
%     A{kk}  = 0.98*A{kk}/max(abs(eig(A{kk})));                              % Normalize excitation matrix 
    A{kk}  = ascale*A{kk}/mean(sum(A{kk}));                                  % Normalize excitation matrix 
    MU{kk} = [gamrnd(1,spike_opts.rate,[N_now,1]); 
              gamrnd(1,spike_opts.rate,[N_nowB,1])];
    B{kk}  = [gamrnd(3,bscale,[N_now,1]); 
              gamrnd(3,bscale,[N_nowB,1])];
    if verbose
        fprintf('.')
    end
end

if verbose
    fprintf('done.\n')
end
% B{1} = 0*B{1}+bscale/2;
A{1}(logical(eye(length(A{1})))) = selfAct*B{kk};
% A{1}(logical(eye(length(A{1})))) = 1.26*B{kk};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up and run marked point-process

evt = cell(batch_num,1);                                                   % Initialize event vector
evm = cell(batch_num,1);                                                   % Initialize mark vector
% extSc = 3;
% inbSc = 1;
% extSc = max(1,3+randn(N_now,1));
extSc = max(0.3,1+0.3*randn(N_now,1));
% extSc = 4.5;
% inbSc = 1+0.1*randn(N_now,1);
inbSc = extSc/2;
% inbSc = extSc./(2+0.2*randn(N_now,1));
gamma    = @(t) exp(-t);                                                   % Set up an anonymous function for the exponential function
gammalen = 10;                                                             % Allow gamma(t)=0 if t>gammalen to save time (truncates history)
wrand    = @(w) sum([0;cumsum(w(:))]<rand(1)*sum(w(:)));                   % Weighted random integer generation
genMKF   = @(m,a) @(t,ht,hm) wrand(m + a(:,hm)*gamma(t-ht));
genCIF   = @(summ,suma) @(tcurr,ht,hm) summ+suma(hm(:,1))*gamma(tcurr-ht);

if verbose
    fprintf('Simulating Hawkes process...')
end

alpha = 3;
rectfun = @(z) log(1 + exp(alpha*z));

for kk = 1:batch_num
    if discrete_flag
        yt      = zeros(size(MU{kk}));                                     % Initialize the vector of rate deviations
        zt      = zeros(size(MU{kk}));                                     % Initialize the vector of rate deviations
        evt{kk} = [];                                                      % Initialize the event times
        evm{kk} = [];                                                      % Initialize the event marks
        yt      = yt+5;
        for tt = 1:spike_opts.nt
%             xt = rand(size(MU{kk}))<1-exp(-rectfun(((zt-yt)*mean(MU{kk})+1).*MU{kk}).*spike_opts.dt);    % See which neurons fired currently
            xt = rand(size(MU{kk}))<1-...
                          exp(-rectfun((zt-yt+1)).*MU{kk}.*spike_opts.dt); % See which neurons fired currently
            zt = gamma(extSc*spike_opts.dt).*zt + A{kk}*(xt);              % Update rate variables            
            yt = gamma(inbSc*spike_opts.dt).*yt + B{kk}.*xt;
            evt{kk} = cat(1,evt{kk},(tt-1+rand(sum(xt),1))*spike_opts.dt); % Spike times uniformly spread out in the time bin
            evm{kk} = cat(1,evm{kk},vec(find(xt)));                        % Store which neurons fired
        end
    elseif ~discrete_flag
        MUs = sum(MU{kk});                                                 % Total base rate
        As  = sum(A{kk},1);                                                % Sum of potential excitation along the rows
        CIF = genCIF(MUs,As);
        MKF = genMKF(MU{kk},A{kk});                                        
        [evt{kk},evm{kk}] = markpointproc(CIF,[],MKF,tmax,inf,1,gammalen); % Run marked point process
    end
    if verbose
        fprintf('.')
    end
end

if verbose
    fprintf('done.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Concatenate all events

evt_cell = [];
evm_cell = [];
evt_bg   = [];
evm_bg   = [];

if verbose
    fprintf('Concatenating events...')
end

for kk = 1:batch_num
    N_now  = min(batch_sz, N_node-(kk-1)*batch_sz);                        % Number of neurons to include in this batch
    
    evt_cell = cat(1, evt_cell, evt{kk}(evm{kk}<=N_now));                  % Concatenate in soma events
    evm_cell = cat(1, evm_cell, (kk-1)*batch_num+evm{kk}(evm{kk}<=N_now)); % Concatenate in soma marks
    evt_bg   = cat(1, evt_bg, evt{kk}(evm{kk}>N_now));                     % Concatenate in background events
    evm_bg   = cat(1, evm_bg,(kk-1)*batch_bg+evm{kk}(evm{kk}>N_now)-N_now);% Concatenate in background marks
        
    if verbose
        fprintf('.')
    end
end

if verbose
    fprintf('done\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Partition events into bins

if verbose
    fprintf('Turning events to bin counts...')
end

S.soma = binSpikeTrains(evt_cell, evm_cell, N_node, spike_opts.dt, ...
                                                            spike_opts.nt);% Bin the continuous-time marked spike events into a spike-count matrix
S.bg   = binSpikeTrains(evt_bg, evm_bg, N_bg, spike_opts.dt,spike_opts.nt);% Bin the continuous-time marked spike events into a spike-count matrix

if verbose
    fprintf('done.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

if nargout > 1
    net_params.A  = A;
    net_params.mu = MU;
    varargout{1}  = net_params;
end

if nargout > 2
    ev.evt       = evt;
    ev.evm       = evm;
    varargout{2} = ev;
end

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
