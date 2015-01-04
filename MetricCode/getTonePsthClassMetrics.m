function metricsClassify = getTonePsthClassMetrics(data,nAdapReps)
% function metricsClassify = getTonePsthClassMetrics(data,[nAdapReps])
%
% Get tone PSTH classification metrics as described in Roos and May, 2012.
%
%
% Given the spike times in response to multiple repetitions of a tone
% stimulus this function computes metrics used to classify VCN units as
% primarylike or chopper, and their subtypes.
%
% INPUTS
% ------
% data.spikes	An Nx2 matrix defining the spike times of N spikes. The
%					first column contains the number of the stimulus repetition
%					(typically between 1 and several hundred) and the second
%					column contains the spike time in seconds.
% data.stimulusOnDur		The duration of the stimulus in seconds. Stimuli
%								are expected to be at the beginning of each
%								repetition/trial.
% data.stimulusOffDur	The duration silent period that follows the
%								stimulus in each repetition/trial, in seconds.
% nAdapReps		The first nAdapReps repetitions/trials are not used for
%					calculating metrics, to mitigate adaptation effects. The
%					default is 30 repetitions (10% of 300).
%
%
% OUTPUTS
% -------
% metricsCalssify		A structure containing metrics needed to navigate the
%							classifation tree in Roos and May, 2012.
%
%
% Note: ISI and CV values are computed using 27 spikes per "bin" (variable
% duration bins), which is 10% of 300-30=270 lines.  Should probably modify
% this to use 10% of however many repetitions there are (minus the
% adaptation repetitions).
%
% Metrics computed here match those described in the following paper:
%
% Roos MJ and May BJ, 2012. Classification of unit types in the
% anteroventral cochlear nucleus of laboratory mice. Hearing Research. 289,
% 13-26.
%
%
% Matt Roos, 12/24/2014


% Determine how many early stimulus presentations to ignore due to
% adaptation effects.
if nargin < 2
	nAdapReps = 30;
end

% Get first and second spikes for each stimulus repitition
latChase = latencyChase2007(data,1e-4);	% Chase default is 1e-6
range = [latChase-0.0005 latChase+0.01];	% time range over which to search for first spikes
spikes = data.spikes;	% for convenience
nTotalReps = max(spikes(:,1));
ixFirstSpikes = nan(nTotalReps,1);
ixSecondSpikes = nan(nTotalReps,1);
for lineCnt = nAdapReps+1:nTotalReps
	ixSpikes = find(spikes(:,1)==lineCnt);
	if isempty(ixSpikes)
		continue;
	end
	ixInWindow = find(spikes(ixSpikes,2) >= range(1) & ...
							spikes(ixSpikes,2) <= range(2));
	if ~isempty(ixInWindow)
		[~,ixClosestSpike] = min(abs(spikes(ixSpikes(ixInWindow),2)-latChase));
		ixFirstSpikes(lineCnt) = ixSpikes(ixInWindow(ixClosestSpike));
		ix = find(spikes(ixSpikes,2) > spikes(ixFirstSpikes(lineCnt),2));
		if ~isempty(ix)
			ixSecondSpikes(lineCnt) = ixSpikes(ix(1));
		end
	end
end

% Compute FSL metrics
ixFirstSpikes = ixFirstSpikes(~isnan(ixFirstSpikes));	% remove NaNs
ixSecondSpikes = ixSecondSpikes(~isnan(ixSecondSpikes));	% remove NaNs
fslMed = median(spikes(ixFirstSpikes,2));	% shorthand for below
fslSD = std(spikes(ixFirstSpikes,2));

% Jensen-Shannon (Kullback-Leibler) divergence test for PDF differences.
% First build the PDFs.
xmin = min([spikes(ixFirstSpikes,2); spikes(ixSecondSpikes,2)]);
xmax = max([spikes(ixFirstSpikes,2); spikes(ixSecondSpikes,2)]);
klStepSize = 0.0002;
x = (xmin:klStepSize:(xmax+klStepSize))';
pdf1 = (hist(spikes(ixFirstSpikes,2),x))';
pdf2 = (hist(spikes(ixSecondSpikes,2),x))';
try
	jsDivergence = kldiv(x,pdf1/sum(pdf1)+eps,pdf2/sum(pdf2)+eps,'js');
catch
	jsDivergence = NaN; % KL (JS) divergence
end

% Compute CV and ISI functions. Use fixed number of ISI samples rather
% than fixed bin width.
[cv,t,isi,~,~,~] = get_CV(data,[0.00 0.07],0.00015,0.00005,(nAdapReps+1:300),false,27);
% Remove NaNs
ix = find(~isnan(isi));
t = t(ix);
isi = isi(ix);
cv = cv(ix);

% Extract CV and ISI metrics and compute ISI slope/transience metrics.
if isempty(cv(~isnan(cv)))
	% get_CV() failed. this might be an onset unit
	cvEarly = NaN;
	cvLate = NaN;
	isiAvgSlope = NaN;
	isiMaxSlope = NaN;
	isiAdapRatio = NaN;
else
	% Get estimates of early CV and ISI
	ix = find(t>=fslMed-0.0001 & t<=fslMed+0.0001);	
	cvEarly =  nanmedian(cv(ix));			% CV at onset +/- 0.1 ms	
	isiEarly = nanmedian(isi(ix));	% ISI at onset +/- 0.1 ms
	
	% Get estimates of middle ISI and ISI adaptation ratio
	ix = find(t>=fslMed+0.0010 & t<=fslMed+0.020);
	isiMiddle = nanmedian(isi(ix));	% ISI at onset + 10-20 ms	
	isiAdapRatio = isiMiddle/isiEarly;	% ISI Adaptation Ratio
		
	% Get estimates of late CV and ISI, and ISI average slope
	ix = find(t>=fslMed+0.020 & t<=fslMed+0.030);
	cvLate = nanmedian(cv(ix));	% Delayed CV (onset + 20-30 ms)
	isiLate = nanmedian(isi(ix));			% ISI at onset + 20-30 ms	
	isiAvgSlope = (isiLate-isiEarly)/0.025;
	
	% Get estimate of max ISI slope (over 1 ms duration, from onset to 15 ms)
	try
		Ts = min(diff(t));	% sampling rate
		N = round(0.001/Ts);			% number of samples in 1 ms averaging filter
		hFilt1ms = ones(N,1)/N;
		isiFilt = filter(hFilt1ms,1,isi);
		[~,startSample] = min(abs(fslMed-t));
		[~,stopSample] = min(abs(0.015-t));
		startSample = startSample + (N-1);	% compensate for filter delay
		stopSample = stopSample + (N-1);	% compensate for filter delay
		hFiltDiff = [1; zeros(N-1,1); -1];	% difference filter with span of 1 ms
		isiDiff = filter(hFiltDiff,1,isiFilt(startSample:stopSample));
		isiMaxSlope = max(isiDiff(N+1:end))/(N*Ts);
		if isempty(isiMaxSlope)
			isiMaxSlope = NaN;	% unit latency was too high?
		end
	catch
		isiMaxSlope = NaN;
	end
end

% Put all metrics used for classification tree into one structure.
metricsClassify.CV_E = cvEarly;
metricsClassify.CV_L = cvLate;
metricsClassify.JSD = jsDivergence;
metricsClassify.FSL_SD = fslSD * 1000;	% convert from sec to ms
metricsClassify.ISI_AR = isiAdapRatio;
metricsClassify.ISI_MS = isiMaxSlope;
metricsClassify.ISI_AS = isiAvgSlope;

