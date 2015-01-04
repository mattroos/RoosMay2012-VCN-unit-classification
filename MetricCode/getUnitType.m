function unitType = getUnitType(metrics)
% function unitType = getUnitType(metrics)
%
% Based on unit PSTH metrics, classify unit type according to Roos and May
% 2012.


unitType.genericType = 'NotClassified';
unitType.subType = 'NotClassified';

if metrics.CV_E > 0.25
	unitType.genericType = 'Primarylike';
	if metrics.JSD > 0.85 && metrics.FSL_SD < 0.5
		unitType.subType = 'Pri-N';
	elseif metrics.JSD <= 0.85 && metrics.FSL_SD >= 0.5
		unitType.subType = 'Pri';
	end	
elseif metrics.CV_E < 0.2
	unitType.genericType = 'Chopper';
	if metrics.ISI_AR > 2.3 && metrics.ISI_MS > 0.9
		unitType.subType = 'Ch-T';
	elseif metrics.ISI_AR <= 2.3 && metrics.ISI_MS <= 0.9
		unitType.genericType = 'Ch-S';		
		if metrics.ISI_AS < 0.04 && metrics.CV_L < 0.3
			unitType.subType = 'Ch-S';
		elseif metrics.ISI_AS >= 0.04 && metrics.CV_L >= 0.3
			unitType.subType = 'Ch-SA';			
		end
	end
end

	
	