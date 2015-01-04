function [CV,time,isiMean,isiStd,nIsi,vecStr] = get_CV(data,range,binDur,varargin)
% function [CV,time,isiMean,isiStd,nIsi,vecStr] = get_CV(data,range,binDur,[stepDur],[repNums],[bFixBinDur],[nIsiSamp])
%
% This function will calculate the coefficient of variation according to the
% algorithm described in Young et al 1988.  Young recommend a binDur of
% 0.0001 or 0.0002 s.  CV is computed from data within the specified range.
% (two-element vector of start and stop isiTimes).
%
% Input
% -----
% data:				Input data structture. See getTonePsthClassMetrics() for a
%						description.
% range:				Two-element vector defining time range over which CV is
%						computed.
% binDur:			Duration of bins (in seconds)
% stepDur:			Temporal step size between bin time. Default is
%						binDur.
% repNums:			Array of repetition numbers from which CV is computed.
% bFixBinDur:		If true, uses a fixed bin duration.
%						If false, uses a fixed number of ISI samples (nIsiSamp).
%						In this case the average of nIsiSamp consecutive
%						ISIs are associated with the average time of those ISIs.
%						The entire ISI function is then binned and interpolated
%						to achieve uniform sample spacing defined by stepDur.
%						Default value is true.
% nIsiSamp:			Number of ISI samples to use when bFixBinDur is false.
%						*This must be an ODD number.* Default is 27 (10% of 270,
%						which assumes 300 reps with the first 30 discarded).
%
% Output
% ------
% CV:			ISI coefficient of variation (=isiStd/isiMean)
% time:		Time at center of each bin
% isiMean:	ISI mean
% isiStd:	ISI standard deviation
% nIsi:		number of ISIs with first spike in the bin
% vecStr:	Raster vector strength
%
% ISI-depedent values are computed from intervals in which the first spike
% occurs within the bin, i.e., two or more spikes are needed to compute a
% non-zero ISI standard deviation and CV.

% TODO: make these input parameters?
bRemoveExtremeISIs = true;	% currently removes values > 3 SDs from the mean
bUseMedianISIs = false;		% use median ISI for CV and vecStr calculations (rather than mean ISI)


if (nargin < 3 || nargin > 7)
	error('Incorrect number of input argument.');
end

switch nargin
	case 3
		stepDur = binDur;
		repNums = [];
		bFixBinDur = true;
		nIsiSamp = NaN;
	case 4
		stepDur = varargin{1};
		repNums = [];
		bFixBinDur = true;
		nIsiSamp = NaN;
	case 5
		stepDur = varargin{1};
		repNums = varargin{2};
		bFixBinDur = true;
		nIsiSamp = NaN;
	case 6
		stepDur = varargin{1};
		repNums = varargin{2};
		bFixBinDur = varargin{3};
		nIsiSamp = 27;
	case 7
		stepDur = varargin{1};
		repNums = varargin{2};
		bFixBinDur = varargin{3};
		nIsiSamp = varargin{4};
end

% Note on get_isis() below:  The ISI time for each pair of spikes is the
% time of the first spike of the pair, if the last input variable is
% 'false'.  Otherwise the ISI time is the mean of the two spike times. It's
% generally preferable to use the first spike time as this help
% differentiate between PriN and ChT units when computing CV. (However, when
% using ISIs to compute raster vector strength as done below, it may be
% better to use the mean spike time?) --MJ Roos
[isi, isiTimes, singles] = get_isis(get_spike_times(data,range,repNums),false);
spikeTimes = data.spikes(:,2);

edges = (range(1):stepDur:range(2))';
time = edges(1:end-1) + stepDur/2;
isiMean = nan(size(time));
isiMedian = nan(size(time));
isiStd = nan(size(time));
CV = nan(size(time));
nIsi = nan(size(time));
vecStr = nan(size(time));

halfBinSize = binDur/2;

if (bFixBinDur)
	for iTime = 1:length(time)
		% find the indices of the intervals beginning within each bin regardless of presentation number
		ind_times = find( isiTimes >= (time(iTime)-halfBinSize) & isiTimes < (time(iTime)+halfBinSize) );
		nIsi(iTime) = length(ind_times);

		% pull out those indices that refer to negative intervals or those > 0.05 s long
		ind = find(isi(ind_times) > 0 & isi(ind_times) < 0.05);
		%ind = find(isi(ind_times) > 0 & isi(ind_times) < 0.02);  % a single long ISI can give very misleading CV value
		if length(ind) >= 3	% requiring 3 intervals per bin (Young et al 1988)
			% Compute mean, std, and CV for ISIs beginning in each bin
			
			if(bRemoveExtremeISIs)
				%---- Option 1: Remove the top 5% of ISIs since these may highly skew the CV estimate
% 				isi2 = sort(isi(ind_times(ind)));
% 				isiLength = length(isi2);
% 				isi2(round(0.95*isiLength):end) = [];
% 				isiMean(iTime) = mean(isi2);
% 				isiStd(iTime) = std(isi2);
				
				%----Option 2: Remove values greater than 3 SDs from the mean
				isi2 = isi(ind_times(ind));
				isiMean(iTime) = mean(isi2);
				isiStd(iTime) = std(isi2);
				bExtreme = abs((isi2-isiMean(iTime))/isiStd(iTime))>3;
				isiMean(iTime) = mean(isi2(~bExtreme));
				isiStd(iTime) = std(isi2(~bExtreme));
				
				%---- Option 3: Define allowable range based on interquartile range
				% Not yet implemented, if ever.
		
				if bUseMedianISIs
					isiMedian(iTime) = median(isi2);	% not excluding extremes in median calculation
				end
			else
				isiMean(iTime) = mean(isi(ind_times(ind)));
				isiStd(iTime) = std(isi(ind_times(ind)));
				if bUseMedianISIs
					isiMedian(iTime) = median(isi(ind_times(ind)));
				end
			end
			
			if bUseMedianISIs
				CV(iTime) = isiStd(iTime) / isiMedian(iTime);
			else
				CV(iTime) = isiStd(iTime) / isiMean(iTime);
			end				

			% compute variability of spike timing using vector strength
			if bUseMedianISIs
				ixSpikes = find( spikeTimes >= (time(iTime)) & spikeTimes < (time(iTime)+isiMedian(iTime)) );
				vectors = exp(1i*2*pi*spikeTimes(ixSpikes)/isiMedian(iTime));
			else
				ixSpikes = find( spikeTimes >= (time(iTime)) & spikeTimes < (time(iTime)+isiMean(iTime)) );
				vectors = exp(1i*2*pi*spikeTimes(ixSpikes)/isiMean(iTime));
			end
			vecStr(iTime) = abs(mean(vectors));
			%vecStr(iTime) = 2*length(vectors)*vecStr(iTime).^2;	% Rayleigh statistic
		end
	end

else
	% Compute metrics using variable bin width and fixed number of ISIs per bin
	ix = isiTimes >= range(1) & isiTimes <= range(2);
	isiTimes = isiTimes(ix);
	isi = isi(ix);
	[isiTimes,ix] = sort(isiTimes);
	isi = isi(ix);
			
	% Compute mean ISI and ISI times...
	
	if length(isiTimes) < nIsiSamp + 3
		warning('get_CV(): Not enough ISIs to compute mean ISI and CV values.');
		return;
	end
	
	% Put time and ISI values into Toeplitz matrices to compute mean and
	% stardard deviation.  This method allows extreme values to be easily
	% removed if desired.
	isiTimesToeplitz = flipud(toeplitz(isiTimes(1:nIsiSamp),isiTimes));
	isiTimesToeplitz(:,1:nIsiSamp-1) = [];
	isiToeplitz = flipud(toeplitz(isi(1:nIsiSamp),isi));
	isiToeplitz(:,1:nIsiSamp-1) = [];
	
	if bRemoveExtremeISIs
		%---- Option 1: Remove the top 5% of ISIs since these may highly skew the CV estimate
%		nRemove = round(nIsiSamp*0.05);
% 		[isiToeplitz,ix] = sort(isiToeplitz);
% 		isiToeplitz(end-nRemove+1:end,:) = [];
% 		% remove corresponding ISI times (although probably unnecessary)
% 		for k = 1:length(ix)
% 			isiTimesToeplitz(:,k) = isiTimesToeplitz(ix(:,k),k);
% 		end
% 		isiTimesToeplitz(end-nRemove+1:end,:) = [];
		
		%---- Option 2: Remove values greater than 3 SDs from the mean
		isiFiltMean = mean(isiToeplitz);
		isiFiltStd = std(isiToeplitz);
		for iBin = 1:size(isiToeplitz,2)
			ixExtreme = abs((isiToeplitz(:,iBin)-isiFiltMean(iBin))/isiFiltStd(iBin))>3;
			isiToeplitz(ixExtreme,iBin) = NaN;
			isiTimesToeplitz(ixExtreme,iBin) = NaN;
		end
		
		%---- Option 3: Define allowable range based on interquartile range
		% Not yet implemented, if ever.
	end
	
	isiTimes = nanmean(isiTimesToeplitz);
	isiFiltMean = nanmean(isiToeplitz);
	isiFiltStd = nanstd(isiToeplitz);
	
	% Perform median filtering of ISIs. Assumes odd filter length. Median
	% filtering inherently removes extremes so not doing that here.
	isiFiltMedian = medfilt1(isi,nIsiSamp);
	isiFiltMedian(1:(nIsiSamp-1)/2) = [];
	isiFiltMedian(end-(nIsiSamp-1)/2+1:end) = [];
	
	% bin to get uniform sampling rate
	for iTime = 1:length(time)
		% find the indices of the intervals beginning within each bin regardless of presentation number
		ix = find( isiTimes >= (time(iTime)-halfBinSize) & isiTimes < (time(iTime)+halfBinSize) );
		isiMean(iTime) = mean(isiFiltMean(ix));
		isiStd(iTime) = mean(isiFiltStd(ix));
		isiMedian(iTime) = mean(isiFiltMedian(ix));
	end
	
	% Interpolate to replace NaN values
	ixNaN = isnan(isiMean);
	isiMean(ixNaN) = interp1(time(~ixNaN),isiMean(~ixNaN),time(ixNaN),'linear');
	
	ixNaN = isnan(isiStd);
	isiStd(ixNaN) = interp1(time(~ixNaN),isiStd(~ixNaN),time(ixNaN),'linear');
		
	if bUseMedianISIs
		ixNaN = isnan(isiMedian);
		isiMedian(ixNaN) = interp1(time(~ixNaN),isiMedian(~ixNaN),time(ixNaN),'linear');
		ix =	~isnan(isiMedian.*isiStd);
		CV(ix) = isiStd(ix)./isiMedian(ix);
	else
		ix =	~isnan(isiMean.*isiStd);
		CV(ix) = isiStd(ix)./isiMean(ix);
	end
	
	% Compute vector strength
	for iTime = 1:length(time)
		% compute variability of spike timing using vector strength
		if bUseMedianISIs
			ixSpikes = find( spikeTimes >= (time(iTime)) & spikeTimes < (time(iTime)+isiMedian(iTime)) );
			vectors = exp(1i*2*pi*spikeTimes(ixSpikes)/isiMedian(iTime));
		else
			ixSpikes = find( spikeTimes >= (time(iTime)) & spikeTimes < (time(iTime)+isiMean(iTime)) );
			vectors = exp(1i*2*pi*spikeTimes(ixSpikes)/isiMean(iTime));
		end
		vecStr(iTime) = abs(mean(vectors));
		%vecStr(iTime) = 2*length(vectors)*vecStr(iTime).^2;	% Rayleigh statistic
	end
end
