function [sp,evm] = genNextSpikeTimepoint(spike_opts, n_locs, sp)

% function [sp,evm] = genNextSpikeTimepoint(spike_opts, n_locs, sp)
%
% Function to generate a set of correlated spike trains by using a Hawkes
% process, or a discrete approximation of a Hawkes process.
%
% 2019 - Alex Song and Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample connectivity for Hawkes process

if (~isfield(sp,'A'))
  ascale = 4;
  bscale = 2;
  if (~isfield(spike_opts,'selfact'))||isempty(spike_opts.selfact)
    selfact = 1.26;
  else
    selfact = spike_opts.selfact;                                                        % A scaling activity
  end
  sp.A  = sampSmallWorldMat(spike_opts.K, 10, 0.3, 0.9, spike_opts.burst_mean, n_locs); % Generate random small-world matrix of excitation for the Hawkes process
  sp.A  = ascale*sp.A/mean(sum(sp.A));                                  % Normalize excitation matrix
  sp.MU = single(gamrnd(1,spike_opts.rate,[spike_opts.K,1]));
  sp.B  = single(gamrnd(3,bscale,[spike_opts.K,1]));
  sp.A(logical(eye(length(sp.A)))) = selfact*sp.B;
  sp.A = sparse(sp.A);
end

% Set up and run marked point-process
if (~isfield(sp,'extSc'))
  sp.extSc = max(0.3,1+0.3*randn(spike_opts.K,1));
  sp.inbSc = sp.extSc/2;
  sp.gamma    = @(t) exp(-t);                                              % Set up an anonymous function for the exponential function
  alpha = 3;
  sp.rectfun = @(z) log(1 + exp(alpha*z));
  sp.yt      = zeros(size(sp.MU))+5;                                       % Initialize the vector of rate deviations
  sp.zt      = zeros(size(sp.MU));                                         % Initialize the vector of rate deviations
end

xt = rand(size(sp.MU))<1-exp(-sp.rectfun((sp.zt-sp.yt+1)).*sp.MU.*spike_opts.dt);    % See which neurons fired currently
sp.zt = sp.gamma(sp.extSc*spike_opts.dt).*sp.zt + sp.A*(xt);               % Update rate variables
sp.yt = sp.gamma(sp.inbSc*spike_opts.dt).*sp.yt + sp.B.*xt;
evm = xt;                                               % Store which neurons fired
