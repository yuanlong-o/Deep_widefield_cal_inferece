function spikes = spikesCellToVec(s,dt,numt)


if(nargin<2)
  dt = 0.01;
end

if(nargin<3)
  numt = ceil(max([s{:}])/dt);
end

spikes = false(length(s),numt);
for i = 1:length(s)
  spikes(i,ceil(s{i}/dt)) = 1;
end