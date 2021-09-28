function trace_mat = expression_variation(trace_mat, p_off, min_mod)

% function trace_mat = expression_variation(trace_mat, p_off, min_red)
% 
% Function to take a matrix of traces (N neurons by T time-points) and to
% modulate each neuron's time trace by a value that is 0 with probability
% p_off (that cell has no expression and is invisible to the optics) or by
% a value that is randomly chosen at uniform between min_mod and 1. 
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if (min_mod(1) < 0)||(min_mod(1) > 1)
    error('minimum modulation level must be between zero and one.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modulate the per-trace levels

if numel(trace_mat) > 1
    N = size(trace_mat,1);                                                 % Extract number of cells
elseif numel(trace_mat) == 1
    N = trace_mat;                                                         % Extract number of cells
end
if(length(min_mod)==1)
  x = min_mod + (1-min_mod)*rand(N,1);                                     % Initialize modulation variables as uniform random between min_red and 1
else
  x = gamrnd(min_mod(2),min_mod(1),N,1);
%   x = min_mod(1) + abs(min_mod(2)*randn(N,1)+1-min_mod(1));                % Initialize modulation variables as uniform random between min_red and 1  
end

x = x.*(p_off < rand(N,1));                                                % Each cells to zero fluorescence with probability p_off
if numel(trace_mat) > 1
    trace_mat = bsxfun(@times, trace_mat, x);                              % Modulate the activity of each cell (the rows of trace_mat)
elseif numel(trace_mat) == 1
    trace_mat = x;                                                         % Otherwise just output the multiplicative values
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%