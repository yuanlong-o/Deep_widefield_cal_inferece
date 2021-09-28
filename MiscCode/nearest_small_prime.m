function N = nearest_small_prime(N,P)

% N = nearest_small_prime(N,P)
% 
% Find the closest, larger number to N with no factors greater than P.
% 
% 2018 - Adam Charles (from code by Chris Turnes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check inputs

if isempty(P)
    P = 7;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find next largest value with a good factorization

if numel(N) > 1
    for kk = 1:numel(N)
        N(kk) = nearest_small_prime(N(kk),P);                              % If a vector of values is given, find the 
    end
else
    if abs(N-round(N))>1e-3                                                % Test that N is an integer (or close to it)
        error('N is not even close to an integer!')                        % If not, throw an error
    else                   
        N = round(N);                                                      % Otherwise make sure that N is technically an integer
    end
    if N > 0
        while max(factor(N)) > P
            N = N + 1;                                                     % Iteratively increase N until the maximum factor is P
        end
    else
        warning('N is not a positive number')
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
