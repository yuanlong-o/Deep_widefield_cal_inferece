function out = veclinspace(v1,v2,n)

% out = veclinspace(v1,v2,n)
% 
% A vectorized version of linspace that creates a sequence of n vectors
% between v1 and v2. 
% 
% 2017 - Adam Charles and Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = zeros(length(v1),n);                                                 % Initialize the output at all zeros
for i = 1:length(v1)
    out(i,:) = linspace(v1(i),v2(i),n);                                    % Perform linspace on the i^th component of the vectors
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%