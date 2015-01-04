function spikeTimes = get_spike_times(data,win,varargin)
% function spikeTimes = get_spike_times(pData,win,[repNums])
%
% This function will output a cell structure containing the spike times
% within the designated time window. These data are not sorted (e.g. by
% frequency or attenuation).
%
% Input
% -----
% data:		Input data structture. See getTonePsthClassMetrics() for a
%				description.
% win:		two element vector defining the window (in seconds)
% repNums:	array of repetition numbers for which spike times are extracted (optional)
%
% Output
% ------
% spikeTimes:	cell with spike times for each line in the picture data

if ~isempty(varargin)
	repNums = varargin{1};
	if ~isempty(repNums)
		spikeTimes = cell(1,length(repNums));
	else
		spikeTimes = cell(1,max(data.spikes(:,1)));
		repNums = 1:length(spikeTimes);
	end
else
	spikeTimes = cell(1,max(data.spikes(:,1)));
	repNums = 1:length(spikeTimes);
end

for iRep = 1:length(repNums)
	ixLine = find(data.spikes(:,1) == repNums(iRep));
	ixTime = find(data.spikes(ixLine,2) >= win(1) & data.spikes(ixLine,2) < win(2));
	spikeTimes{iRep} = data.spikes(ixLine(ixTime),2);
end
