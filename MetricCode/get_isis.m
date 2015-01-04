function [isi,times,singletons] = get_isis(spiketimes,fMeanTime)
% function [isi,times,singletons] = get_isis(spiketimes,[fMeanTime])
%
% This function will take a cell array of spike trains and calculate an ISI
% distribution for all arrays combined.  spiketimes is in the form
% output by get_spike_times().
%
% input
% ------
% spiketimes:	The spiketimes, as formatted by get_spike_times()
% fMeanTime:	If false, this function returns the times of the first spike
%					of the pair. If true it returns the mean times.
%
% output
% ------
% isi:	     The ISI durations
% times:      The time of the initial spike for each ISI value, or the mean
%             time of the pair if fMeanTime is true.
% singletons: The number of picture lines (presentations) in which one and
%             only one spike occured.

if nargin < 2
	fMeanTime = false;
end

isi = [];
times = [];
singletons = 0;

for i = 1 : length(spiketimes)    % find the isi's of the spikes gathered from above
	if length(spiketimes{i}) > 1
		int = diff(spiketimes{i})';
		isi = [isi int];
		if ~fMeanTime
			times = [times spiketimes{i}(1:end-1)'];    % the time of the initial spike for each interval
		else
			times = [times (spiketimes{i}(1:end-1)+spiketimes{i}(2:end))'/2];    % the average time of the spikes for each interval
		end
	elseif length(spiketimes{i}) == 1
		singletons = singletons + 1;
	end
end
