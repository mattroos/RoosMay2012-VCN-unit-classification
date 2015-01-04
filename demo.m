% demo.m
%
% Demo units are those of Figure 1 in Roos and May, 2012.

clear

fs = filesep;
dirWork = pwd;
dirCode = [dirWork fs 'MetricCode' fs];
dirData = [dirWork fs 'DemoData' fs];

restoredefaultpath
addpath(dirCode,'-begin');


figure(1); clf;
binSize = 0.1;	% histogram bin size in ms


%% Pri Unit
load([dirData 'Pri' fs 'p0056_u6_01_PST_tone']);

% Compute metrics of Roos and May 2012 and classify unit
metrics = getTonePsthClassMetrics(data,30);
unitType = getUnitType(metrics);

% Plot spike raster
subplot(2,4,1);
plot(data.spikes(:,2)*1000,data.spikes(:,1),'k.','MarkerSize',1);
axis([0 70 0 300]);
ylabel('Presentation number');
title(unitType.subType);

% Plot spike PSTH
subplot(2,4,5);
ix = find(data.spikes(:,1)>30);
[cnt,bins] = hist(data.spikes(ix,2)*1000,0:binSize:100);
bar(bins,cnt/binSize/270*1000,1,'k');
axis([0 70 0 2500]);
ylabel('Spikes/s');
xlabel('Time re tone onset (ms)');


%% PriN Unit
load([dirData 'PriN' fs 'p0027_u4_01_PST_tone']);

% Compute metrics of Roos and May 2012 and classify unit
metrics = getTonePsthClassMetrics(data,30);
unitType = getUnitType(metrics);

% Plot spike raster
subplot(2,4,2);
plot(data.spikes(:,2)*1000,data.spikes(:,1),'k.','MarkerSize',1);
axis([0 70 0 300]);
title(unitType.subType);


% Plot spike PSTH
subplot(2,4,6);
ix = find(data.spikes(:,1)>30);
[cnt,bins] = hist(data.spikes(ix,2)*1000,0:binSize:100);
bar(bins,cnt/binSize/270*1000,1,'k');
axis([0 70 0 2500]);
xlabel('Time re tone onset (ms)');


%% ChS Unit
load([dirData 'ChS' fs 'p0001_u4_01_PST_tone']);

% Compute metrics of Roos and May 2012 and classify unit
metrics = getTonePsthClassMetrics(data,30);
unitType = getUnitType(metrics);

% Plot spike raster
subplot(2,4,3);
plot(data.spikes(:,2)*1000,data.spikes(:,1),'k.','MarkerSize',1);
axis([0 70 0 300]);
title(unitType.subType);

% Plot spike PSTH
subplot(2,4,7);
ix = find(data.spikes(:,1)>30);
[cnt,bins] = hist(data.spikes(ix,2)*1000,0:binSize:100);
bar(bins,cnt/binSize/270*1000,1,'k');
axis([0 70 0 2500]);
xlabel('Time re tone onset (ms)');


%% ChT Unit
load([dirData 'ChT' fs 'p0027_u5_02_PST_tone']);

% Compute metrics of Roos and May 2012 and classify unit
metrics = getTonePsthClassMetrics(data,30);
unitType = getUnitType(metrics);

% Plot spike raster
subplot(2,4,4);
plot(data.spikes(:,2)*1000,data.spikes(:,1),'k.','MarkerSize',1);
axis([0 70 0 300]);
title(unitType.subType);


% Plot spike PSTH
subplot(2,4,8);
ix = find(data.spikes(:,1)>30);
[cnt,bins] = hist(data.spikes(ix,2)*1000,0:binSize:100);
bar(bins,cnt/binSize/270*1000,1,'k');
axis([0 70 0 2500]);
xlabel('Time re tone onset (ms)');

