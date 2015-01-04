function latency = latencyChase2007(data,varargin)
% function latency = latencyChase2007(data,[threshold],[fEachLine])
%
% Estimates the onset latency of a neuron to a particular stimulus using
% the probability techniques of Steve Chase, 2007 (also described in his
% 2006 Hopkins PhD thesis). If no onset is found within 30 ms, the
% function returns NaN. If fEachLine is true then the latency is estimated
% for each repetition of the stimulus.
%
%
% Input
% -----
% data:			Input data structture. See getTonePsthClassMetrics() for a
%					description.
% threshold:	The decision probability threshold.  When probability falls
%					below this value, onset is declared.  Default is 10^-6.
% fEachLine:	If false (the default) all repetitions of the stimuli are
%					assumed to have been generated by identical stimuli and a
%					single latency value is returned (as in Chase 2007). If
%					true, latency estimates are made for each repetitions using
%					the spontaneous rate from all lines.
%
% Output
% ------
% latency:		The estimated onset latency (latencies) to the stimuli.
%
%
% Author: Matt Roos, 12/10/09
%
% Reference:
% Chase and Young, First-spike latency information in single neurons
% increases when referenced to population onset. PNAS, 2007


if nargin > 3
	error('Too many input arguments.');
end
if nargin < 2
	threshold = 10^(-6);
end
if nargin >= 2
	if ~isempty(varargin{1})
		threshold = varargin{1};
	else
		threshold = 10^(-6);
	end
end
if nargin > 2
	fEachLine = varargin{2};
else
	fEachLine = false;
end


nReps = max(data.spikes(:,1)) - min(data.spikes(:,1)) + 1;
spikeTimes = sort(data.spikes(:,2));

	% compute intertrial interval in seconds
	ITI = (data.stimulusOnDur + data.stimulusOffDur);

	if (ITI < data.stimulusOnDur + 0.110)
		delta = ITI - data.stimulusOnDur;
		warning(['Intertrial interval is only ' num2str(delta)  'ms longer than stimulus duration']);
	end

	% get spont rate estimate from last 100 ms of each line
	ixSpontSpikes = find(data.spikes(:,2) >= ITI-0.1);
	spontRate = length(ixSpontSpikes)/nReps/0.1;

	if ~fEachLine
		
		% Compute single latency from all stimulus repetitions...
		
		fDone = false;
		minSpikes = 5;		% number of spikes in shortest analysis window. Steve Chase set it to 5.
		spikeNum = minSpikes+1;
		totalNumSpikes = length(spikeTimes);
		spikeTimes = [0; spikeTimes];	% helps with looping below

		while (~fDone && spikeTimes(spikeNum) < 0.30)
			probabilities = ones(spikeNum-1,1);
			for spikesInWin = minSpikes:(spikeNum-1)
				timeDur = spikeTimes(spikeNum) - spikeTimes(spikeNum-spikesInWin);
				term1 = (nReps*spontRate*timeDur).^(0:(spikesInWin-1));
				probabilities(spikesInWin) = 1 - exp(-nReps*spontRate*timeDur)*sum(term1./factorial(0:(spikesInWin-1)));
			end

			if min(probabilities) <= threshold;
				fDone = true;
				latency = spikeTimes(spikeNum);
			end

			spikeNum = spikeNum + 1;
			if (spikeNum > totalNumSpikes)
				break;
			end
		end

		if ~fDone
			% never found onset
			latency = NaN;
		end

	else
		
		% Compute latencies from each repetition...

		latency = nan(nReps,1);
		for iRep = 1:nReps
			ixSpikes = find(data.spikes(:,1) == iRep);
			spikeTimes = sort(data.spikes(ixSpikes,2));
			
			if isempty(spikeTimes)
				latency(iRep) = NaN;
				continue;
			end
			
			fDone = false;
			minSpikes = 1;		% number of spikes in shortest analysis window. Steve Chase set it to 5.
			spikeNum = minSpikes+1;
			totalNumSpikes = length(spikeTimes);
			spikeTimes = [0; spikeTimes];	% helps with looping below

			while (~fDone && spikeTimes(spikeNum) < 0.30)
				probabilities = ones(spikeNum-1,1);
				for spikesInWin = minSpikes:(spikeNum-1)
					timeDur = spikeTimes(spikeNum) - spikeTimes(spikeNum-spikesInWin);
					term1 = (spontRate*timeDur).^(0:(spikesInWin-1));
					probabilities(spikesInWin) = 1 - exp(-spontRate*timeDur)*sum(term1./factorial(0:(spikesInWin-1)));
				end

				if min(probabilities) <= threshold;
					fDone = true;
					latency(iRep) = spikeTimes(spikeNum);
				end

				spikeNum = spikeNum + 1;
				if (spikeNum > totalNumSpikes)
					break;
				end
			end

			if ~fDone
				% never found onset
				latency(iRep) = NaN;
			end
		end		
	end