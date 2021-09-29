function [evt,evm] = markpointproc(CIF,CIFMAX,MKF,TIMEMAX,NUMMAX,MARKDIM,HISTLEN)
% MARKPOINTPROC
%  [T,M] = MARKPOINTPROC(CIF,CIFMAX,MKF,TIMEMAX,NUMEVMAX,MARKDIM,HISTLEN)
%  generates a marked point process using the conditional intensity
%  function CIF and mark generating function MKF.  Events are generated
%  until either the time TIMEMAX is reached or NUMMAX events are generated.
%  It returns times T and marks M of the generated events, which
%  respectively have dimensions NUMEVTSx1 and NUMEVTSxMARKDIM.
%  
%  The functions CIF, CIFMAX, and MKF should have the following formats:
% 
%  R = CIF(T,PASTT,PASTM) returns the CIF evaluated at time T, given past
%  event times PASTT and past marks PASTM.
% 
%  [R,S] = CIFMAX(T,PASTT,PASTM) returns R, an upper bound on the CIF over
%  the interval (T,T+S), given past event times PASTT and past marks PASTM.
%  If CIF is monotonically nonincreasing between event times then set
%  CIFMAX = [], which is equivalent to R = CIF(T,PASTT,PASTM) and S =
%  HISTLEN.
% 
%  M = MKF(T,PASTT,PASTM) returns a mark for an event at time T, given past
%  event times PASTT and past marks PASTM.  M should be a 1xMARKDIM vector.
%  MKF can be omitted when MARKDIM==0.
%  
%  The optional parameter HISTLEN can be used to limit the length (in time)
%  of history to use for CIF and MKF.  HISTLEN should be sufficiently large
%  to ensure that the effects of an event that far in the past are
%  negligible.  Default is HISTLEN=inf.
% 
% Originally by Michael Moore, modified by Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin < 7 || isempty(HISTLEN)
	HISTLEN = inf;
end

if isempty(CIFMAX)
	CIFMAX = @(t,ht,hm) deal(CIF(t+eps(t),ht,hm),HISTLEN);
end

if isinf(TIMEMAX) && isinf(NUMMAX)
	error('No stopping criteria - both TIMEMAX and NUMMAX are infinite');
end

if isfinite(NUMMAX)
	evt = zeros(NUMMAX,1);                                                 % can pre-allocate
	evm = zeros(NUMMAX,MARKDIM);
else
	evt = zeros(100,1);                                                    % may need more, but start with something
	evm = zeros(100,MARKDIM);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run marked point process

histstart = 1;                                                             % Initialize the history vector
t         = 0;                                                             % Initialize time progress
i         = 0;                                                             % Initialize number of events

evtp = evt(histstart:i);                                                   % 
evmp = evm(histstart:i,:);                                                 % 

while i < NUMMAX
	[cifmaxt,cifmaxintvl] = CIFMAX(t,evtp,evmp);
	tstep = -log(rand)/cifmaxt;                                            % Sample next candidate time-difference
	t1    = t+tstep;                                                       % Calculate new time
	if t1 > TIMEMAX
		break;                                                             % If over the max time, the sampling is over
	elseif tstep >= cifmaxintvl
		t = t+cifmaxintvl;                                                 % 
	else
		t = t1;
		cift = CIF(t,evtp,evmp);
		ratefrac = cift/cifmaxt;

		if ratefrac > 1 % cheap enough to check
			error('CIF evaluated to a value greater than CIFMAX')
		end

		if rand < ratefrac % do we keep it?
			i = i+1;
			evt(i) = t;
			if MARKDIM > 0
				evm(i,:) = MKF(t,evtp,evmp);
			else
				evm(i,:) = zeros(1,0);
			end
			% update history
			while t-HISTLEN > evt(histstart)
				histstart = histstart+1;
			end
			evtp = evt(histstart:i);
			evmp = evm(histstart:i,:);
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

evt = evt(1:i);                                                            % trim excess entries, if necessary (for spike times)
evm = evm(1:i,:);                                                          % trim excess entries, if necessary (for marks)

if nargout <= 1
	evt = [evt,evm];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
