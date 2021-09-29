function adj_mat = sampSmallWorldMat(N_node, K_conn, beta, varargin)

% adj_mat = sampSmallWorldMat(N_nodes, K_conn, beta {rand_opt, self_ex,
%                                                                  n_locs}) 
% 
% Create a random small-world graph using the Watts-Strogatz model for
% small world graphs. There is the extra options to make the connections
% scaled random numbers, to emphasise self-excitation (creating bursting
% behavior) and orienting the connectivity using the neural locations
% (connecting spatially neighboring neurons). 
% Inputs:
%  - N_node   - Number of nodes to simulate. If a scalar, only simulates
%               neural nodes. If a size-2 array, the second element
%               indicates the number of 'background processes' that are
%               correlated to all neurons. 
%  - K_conn   - Number of connections for each neuron
%  - beta     - Probability that a connection is "long" rather than local
%  - rand_opt - OPTIONAL: If a scalar this adds a random component to the
%               connection weights
%  - self_ex  - OPTIONAL: If a scalar, this value is added to all
%               self-excitation connections, increasing the bursting
%               behavior
%  - n_locs   - OPTIONAL: An array of locations for the N_node neurons. If
%               provided, the K_conn initial 'local' connections are chosen
%               based on the nearest spatial neighbors.
%
% Output:
%  - adj_mat  - An (N_node)x(N_node) (or sum(N_node)xsum(N_node)) array
%               connectivity matrix for a 'small-world' type matrix, as per
%               the Watts-Strogatz model.
% 
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if isempty(N_node)
    N_node = 100;                                                          % For testing, an empty node input makes a 100-node network w/ no background
end

if numel(N_node) > 1
    N_bg   = N_node(2);                                                    % Get number of background nodes
    N_node = N_node(1);                                                    % Get number of soma nodes
else
    N_bg   = 0;                                                            % If only one node input, then no background nodes
    N_node = N_node(1);                                                    % Get number of soma nodes
end

if isempty(K_conn)
    K_conn = ceil(2*log(N_node));                                          % If not provided, the sparseness of connectivity is dictated by the log of the number of nodes
end

if isempty(beta)
    beta = 0.2;                                                            % If beta is not provided then the default connection transfer probability is 0.2
end

if nargin > 3
    rand_opt = varargin{1};
else
    rand_opt = 0;                                                          % If no random variance provided, then the connections all have unit weight
end

if nargin > 4
    self_ex = varargin{2};
else
    self_ex = 4;                                                           % If no self-excitation value provided, then set to 4
end

if nargin > 5
    n_locs = varargin{3};
else
    n_locs = [];                                                           % If no spatial information, don't use it?
end

if isempty(n_locs)||(size(n_locs,1) ~= N_node(1))
    use_locs = false;
elseif size(n_locs,1) == N_node(1)
    use_locs = true;
else
    error('Input an incompatible location matrix')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run algorithm

if use_locs
    adj_mat = 0;
    for ll = 1:size(n_locs,2)
        adj_mat = adj_mat +bsxfun(@plus, n_locs(:,ll), -n_locs(:,ll).').^2;% Seed adjacency matrix by claculating distances between neurons one dimension at a time
    end
    adj_mat = sqrt(adj_mat);                                               % Doesn't relly affect much, but make sure it's the distance is the l_2 distance
    for kk = 1:size(adj_mat,1)
        [~,IX]        = sort(adj_mat(kk,:),'ascend');                      % Sort the distances to the kk^th neueon
        adj_mat(kk,:) = 0;                                                 % Reset the connections
        adj_mat(kk,IX(1:K_conn)) = 1;                                      % Connect the closest K_conn neurons to that neuron
    end
else
    adj_mat = toeplitz([ones(1,K_conn/2), zeros(1,N_node-K_conn/2)]);          % Initialize the adjacency matrix for a lattice graph
end

for kk = 1:N_node                                                          % Interate through nodes
    switch_flag = rand(K_conn,1)<beta;                                     % Sample K random numbers to test which connections change
    n_switch    = sum(switch_flag);                                        % Number of connections to change 
    new_cons    = randsample(sum(adj_mat(kk,:) == 0), n_switch);           % Randomly select where the new connections should go
    nc_1locs    = find(adj_mat(kk,:) == 1);                                % Get location of connections from the current node
    nc_0locs    = find(adj_mat(kk,:) == 0);                                % Get location of NO connections from the current node
    adj_mat(kk, nc_0locs(new_cons))    = 1;                                % Turn on new connections
  try
    adj_mat(kk, nc_1locs(switch_flag)) = 0;                                % Turn off where those connections used to be
  end
end

if N_bg > 0
    adj_mat = [[adj_mat, zeros(N_node,N_bg) ];
               [ones(N_bg,N_node), eye(N_bg)]];                            % Make the background nodes more generally influenced by the other neurons to induce more correlations
end

adj_mat = adj_mat + self_ex*eye(size(adj_mat));                            % Add a diagonal term that encourages bursts
% adj_mat = max(adj_mat + self_ex*eye(size(adj_mat)),0);                     % Add a diagonal term that encourages bursts
adj_mat = adj_mat.*(1-rand_opt + rand_opt.*(0.1+0.9*rand(size(adj_mat)))); % Make the actual weights random between 0.1 and 1

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%