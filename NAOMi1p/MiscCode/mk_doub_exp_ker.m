function h = mk_doub_exp_ker(t_on, t_off, A, dt, varargin)

% h = mk_doub_exp_ker(t_on, t_off, A, dt)
% 
% Function to create a double-exponential kernel. The kernel is defined as
%
% h(t) = A*(1-exp(-t/t_on))*exp(-t/t_off)
%
% The function is evaluated at increments of dt from t=0 to t=t_max defined
% as the time where the kernel falls below 1e-3 of the max value of the
% kernel.
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs

if nargin > 4
    dext_type = varargin{1};
else
    dext_type = 'mult';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make kernel
% 
% t_on    = 1./t_on;                                                         % Invert the first constant t_on to make the rest of the code simpler
% t_off   = 1./t_off;                                                        % Invert the second constant t_off to make the rest of the code simpler

switch dext_type
    case 'times'
        loc_max = log((t_off+t_on)./t_off)./t_on;                          % Find the location of the peak of the kernel
        dexp    = @(z) A*(1 - exp(-t_on*z)).*exp(-t_off*z);                % Create an anonymous function that evaluates the kernel at a set of points z
    case 'plus'
        if numel(A) == 1
            B = A;
        else
            B = A(2);
            A = A(1);
        end
        t1      = min(t_off,t_on);                                         % Make sure that the first time constant is faster
        t2      = max(t_off,t_on);                                         % Make sure that the second time constant is slower
        loc_max = log(A*t1/(B*t2))/(t1 - t2);                              % Find the location of the peak of the kernel
        dexp    = @(z) A*exp(-t1*z)-B*exp(-t2*z);                          % Create an anonymous function that evaluates the kernel at a set of points z
    case 'min'
        if numel(A) == 1
            B = A;
        else
            B = A(2);
            A = A(1);
        end
        dexp    = @(z) min(A*(1 - exp(-t_on*z)),B*exp(-t_off*z));          % Create an anonymous function that evaluates the kernel at a set of points z
        loc_max = fmincon(@(z)-dexp(z),0.5,[],[],[],[],0,[]);              % Find the location of the peak of the kernel
        max_val = dexp(loc_max);  
        t_max = -log(max_val*1e-3/A)./t_off;  
        tz = vec(0:dt:(t_max(1)+dt));
%         figure(2); plot(tz, [A*(1 - exp(-t_on*tz)),B*exp(-t_off*tz),dexp(tz)])
    otherwise
        loc_max = log((t_off+t_on)./t_off)./t_on;                          % Find the location of the peak of the kernel
        dexp    = @(z) A*(1 - exp(-t_on*z)).*exp(-t_off*z);                % Create an anonymous function that evaluates the kernel at a set of points z
end

max_val = dexp(loc_max);                                                   % Find the value at the maximum of the kernel


t_max = -log(max_val*1e-3/A)./t_off;                                       % Find approximately where the kernel dips below 1e-3 of its maximum to determine the extent of the kernel
h     = dexp(0:dt:(t_max(1)+dt));                                          % Evaluate the kernel at each dt up to the maximum time where the kernel is above 1e-3 of its max

% figure(3); plot(0:dt:(t_max(1)+dt),h)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%