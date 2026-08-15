function metrics = run_burst_analysis(dataFile, outputDirectory, config)
%RUN_BURST_ANALYSIS Detect spikes, bursts, and return-map statistics.
%   METRICS = RUN_BURST_ANALYSIS(DATAFILE, OUTPUTDIRECTORY, CONFIG)
%   accepts numeric tabular data with time in column 1 and a signal in the
%   column selected by CONFIG.SignalColumn.

if nargin < 3
    config = struct();
end

config = applyDefaults(config, struct( ...
    'SignalColumn', 2, ...
    'SignalScale', 1, ...
    'DetectValleys', true, ...
    'ThresholdFactor', 1.5, ...
    'EnableAdaptiveThreshold', false, ...
    'AdaptiveThresholdFactors', [2.0, 1.5, 1.2, 1.0, 0.8, 0.6], ...
    'MinimumAdaptiveSpikes', 50, ...
    'MaximumAdaptiveSpikes', 1000, ...
    'MaxISIWithinBurstSeconds', 0.04, ...
    'MinimumSpikesPerBurst', 3, ...
    'BurstingMeanSizeThreshold', 3.0, ...
    'DisplayWindowSeconds', 10, ...
    'DatasetLabel', 'Signal', ...
    'OutputStem', 'burst_analysis'));

if ~isfile(dataFile)
    error('BNCAnalysis:MissingInput', 'Input file not found: %s', dataFile);
end
if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end

fprintf('Loading %s data from %s\n', config.DatasetLabel, dataFile);
data = readmatrix(dataFile);
if size(data, 2) < config.SignalColumn
    error('BNCAnalysis:InvalidInput', ...
        'Signal column %d is not present in the input file.', ...
        config.SignalColumn);
end

time = data(:, 1);
rawSignal = data(:, config.SignalColumn) * config.SignalScale;
validRows = isfinite(time) & isfinite(rawSignal);
time = time(validRows);
rawSignal = rawSignal(validRows);

if numel(time) < 10
    error('BNCAnalysis:InsufficientData', ...
        'At least ten finite samples are required.');
end

[time, uniqueIndex] = unique(time, 'stable');
rawSignal = rawSignal(uniqueIndex);
samplingRate = 1 / median(diff(time));
duration = time(end) - time(1);

missingMask = isnan(rawSignal);
if any(missingMask)
    rawSignal = fillmissing(rawSignal, 'linear');
end

centeredSignal = rawSignal - median(rawSignal);
standardDeviation = std(centeredSignal);
robustDeviation = median(abs(centeredSignal - median(centeredSignal))) / 0.6745;
positiveNoiseEstimates = [standardDeviation, robustDeviation];
positiveNoiseEstimates = positiveNoiseEstimates( ...
    isfinite(positiveNoiseEstimates) & positiveNoiseEstimates > 0);
if isempty(positiveNoiseEstimates)
    error('BNCAnalysis:ConstantSignal', ...
        'The input signal has no measurable variation.');
end
noiseLevel = min(positiveNoiseEstimates);

if range(centeredSignal) > 10 * noiseLevel
    lowerPercentile = prctile(centeredSignal, 0.1);
    upperPercentile = prctile(centeredSignal, 99.9);
    percentileRange = upperPercentile - lowerPercentile;
    if percentileRange > 0
        normalizedSignal = (centeredSignal - lowerPercentile) / ...
            percentileRange - 0.5;
    else
        normalizedSignal = centeredSignal / (3 * noiseLevel);
    end
else
    normalizedSignal = centeredSignal / (3 * noiseLevel);
end

if config.DetectValleys
    detectionSignal = -normalizedSignal;
else
    detectionSignal = normalizedSignal;
end

absoluteDeviation = abs(detectionSignal - median(detectionSignal));
medianAbsoluteDeviation = median(absoluteDeviation);
if medianAbsoluteDeviation > 0
    normalizedNoise = medianAbsoluteDeviation / 0.6745;
else
    normalizedNoise = std(detectionSignal);
end
backgroundSignal = detectionSignal(abs(detectionSignal) < 2 * normalizedNoise);
if ~isempty(backgroundSignal) && std(backgroundSignal) > 0
    normalizedNoise = min(normalizedNoise, std(backgroundSignal));
end

thresholdFactorUsed = config.ThresholdFactor;
detectionThreshold = thresholdFactorUsed * normalizedNoise;
minimumPeakDistance = max(1, round(0.002 * samplingRate));

[detectedPeaks, peakLocations] = findpeaks(detectionSignal, ...
    'MinPeakHeight', detectionThreshold, ...
    'MinPeakDistance', minimumPeakDistance);

if config.EnableAdaptiveThreshold && ...
        numel(peakLocations) < config.MinimumAdaptiveSpikes
    candidateCounts = zeros(size(config.AdaptiveThresholdFactors));
    for candidateIndex = 1:numel(config.AdaptiveThresholdFactors)
        candidateThreshold = config.AdaptiveThresholdFactors(candidateIndex) * ...
            normalizedNoise;
        [~, candidateLocations] = findpeaks(detectionSignal, ...
            'MinPeakHeight', candidateThreshold, ...
            'MinPeakDistance', minimumPeakDistance);
        candidateCounts(candidateIndex) = numel(candidateLocations);
    end
    acceptable = find(candidateCounts >= config.MinimumAdaptiveSpikes & ...
        candidateCounts <= config.MaximumAdaptiveSpikes, 1, 'first');
    if ~isempty(acceptable)
        thresholdFactorUsed = config.AdaptiveThresholdFactors(acceptable);
        detectionThreshold = thresholdFactorUsed * normalizedNoise;
        [detectedPeaks, peakLocations] = findpeaks(detectionSignal, ...
            'MinPeakHeight', detectionThreshold, ...
            'MinPeakDistance', minimumPeakDistance);
    end
end
spikeTimes = time(peakLocations);
if config.DetectValleys
    spikeAmplitudes = -detectedPeaks;
else
    spikeAmplitudes = detectedPeaks;
end

interSpikeIntervals = diff(spikeTimes);
[spikesPerBurst, burstStartTimes] = detectBursts( ...
    spikeTimes, interSpikeIntervals, ...
    config.MaxISIWithinBurstSeconds, config.MinimumSpikesPerBurst);

if isempty(interSpikeIntervals)
    meanISI = NaN;
    cvISI = NaN;
else
    meanISI = mean(interSpikeIntervals);
    cvISI = std(interSpikeIntervals) / meanISI;
end

numberOfBursts = numel(spikesPerBurst);
if numberOfBursts > 0
    meanBurstSize = mean(spikesPerBurst);
    burstFrequency = numberOfBursts / duration;
else
    meanBurstSize = 0;
    burstFrequency = 0;
end

intraBurstIntervals = interSpikeIntervals( ...
    interSpikeIntervals <= config.MaxISIWithinBurstSeconds);
interBurstIntervals = interSpikeIntervals( ...
    interSpikeIntervals > config.MaxISIWithinBurstSeconds);
if isempty(intraBurstIntervals)
    intraBurstFrequency = 0;
else
    intraBurstFrequency = mean(1 ./ intraBurstIntervals);
end
if isempty(interBurstIntervals)
    interBurstInterval = 0;
else
    interBurstInterval = mean(interBurstIntervals);
end

if isempty(spikeAmplitudes)
    meanAmplitude = NaN;
    standardDeviationAmplitude = NaN;
    medianAmplitude = NaN;
else
    meanAmplitude = mean(spikeAmplitudes);
    standardDeviationAmplitude = std(spikeAmplitudes);
    medianAmplitude = median(spikeAmplitudes);
end

isBursting = meanBurstSize >= config.BurstingMeanSizeThreshold && ...
    isfinite(cvISI) && cvISI > 1;

metrics = struct();
metrics.datasetLabel = config.DatasetLabel;
metrics.dataFile = dataFile;
metrics.samplingRateHz = samplingRate;
metrics.durationSeconds = duration;
metrics.totalSpikes = numel(spikeTimes);
metrics.meanFiringRateHz = numel(spikeTimes) / duration;
metrics.meanISISeconds = meanISI;
metrics.cvISI = cvISI;
metrics.numberOfBursts = numberOfBursts;
metrics.meanBurstSize = meanBurstSize;
metrics.burstFrequencyHz = burstFrequency;
metrics.intraBurstFrequencyHz = intraBurstFrequency;
metrics.interBurstIntervalSeconds = interBurstInterval;
metrics.meanNormalizedAmplitude = meanAmplitude;
metrics.standardDeviationAmplitude = standardDeviationAmplitude;
metrics.detectionThreshold = detectionThreshold;
metrics.thresholdFactorUsed = thresholdFactorUsed;
metrics.isBursting = isBursting;

fprintf('Detected %d spikes and %d bursts.\n', ...
    metrics.totalSpikes, metrics.numberOfBursts);

analysisFigure = figure('Name', [config.DatasetLabel, ' Burst Analysis'], ...
    'Color', 'w', 'Position', [100, 100, 1400, 800]);

windowStart = time(1);
windowEnd = min(windowStart + config.DisplayWindowSeconds, time(end));
windowMask = time >= windowStart & time <= windowEnd;
windowTime = time(windowMask);
windowSignal = normalizedSignal(windowMask);
windowSpikeMask = spikeTimes >= windowStart & spikeTimes <= windowEnd;

subplot(2, 3, 1);
plot(windowTime, windowSignal, 'b', 'LineWidth', 0.8); hold on;
plot(spikeTimes(windowSpikeMask), spikeAmplitudes(windowSpikeMask), ...
    'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
if config.DetectValleys
    plottedThreshold = -detectionThreshold;
else
    plottedThreshold = detectionThreshold;
end
yline(plottedThreshold, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Normalized Amplitude'); grid on;
title(sprintf('A. Signal and Detected Spikes (%g s)', ...
    config.DisplayWindowSeconds));

subplot(2, 3, 2);
if numel(interSpikeIntervals) >= 2
    histogram(interSpikeIntervals, 50, 'Normalization', 'pdf', ...
        'FaceColor', [0.2, 0.4, 0.8], 'EdgeColor', 'none'); hold on;
    xline(config.MaxISIWithinBurstSeconds, 'r--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log'); grid on;
    title(sprintf('B. ISI Distribution (CV = %.2f)', cvISI));
    xlabel('ISI (s, log scale)'); ylabel('Probability Density');
else
    displayUnavailable('B. ISI Distribution');
end

subplot(2, 3, 3);
if numel(interSpikeIntervals) >= 3
    scatter(interSpikeIntervals(1:end-1), interSpikeIntervals(2:end), ...
        20, 'k', 'filled', 'MarkerFaceAlpha', 0.3); hold on;
    set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; axis square;
    xlabel('ISI_n (s)'); ylabel('ISI_{n+1} (s)');
    title('C. ISI Return Map');
else
    displayUnavailable('C. ISI Return Map');
end

subplot(2, 3, 4);
if ~isempty(spikesPerBurst)
    histogram(spikesPerBurst, 'BinMethod', 'integers', ...
        'FaceColor', [0.8, 0.2, 0.2], 'EdgeColor', 'black'); hold on;
    xline(meanBurstSize, 'b--', 'LineWidth', 2);
    xline(config.BurstingMeanSizeThreshold, 'r--', 'LineWidth', 2);
    xlabel('Spikes per Burst'); ylabel('Count'); grid on;
    title(sprintf('D. Burst Size Distribution (N = %d)', numberOfBursts));
else
    displayUnavailable('D. Burst Size Distribution');
end

subplot(2, 3, 5);
if ~isempty(spikeTimes)
    plot(spikeTimes, 1:numel(spikeTimes), 'k.', 'MarkerSize', 8);
    xlabel('Time (s)'); ylabel('Spike Index'); grid on;
    title('E. Spike-Time Sequence');
else
    displayUnavailable('E. Spike-Time Sequence');
end

subplot(2, 3, 6);
if ~isempty(spikeAmplitudes)
    histogram(spikeAmplitudes, 30, 'Normalization', 'pdf', ...
        'FaceColor', [0.5, 0.2, 0.8], 'EdgeColor', 'none'); hold on;
    xline(meanAmplitude, 'r--', 'LineWidth', 2);
    xline(medianAmplitude, 'g--', 'LineWidth', 2);
    xlabel('Normalized Amplitude'); ylabel('Probability Density'); grid on;
    title('F. Spike-Amplitude Distribution');
else
    displayUnavailable('F. Spike-Amplitude Distribution');
end

sgtitle(sprintf('%s | Spikes: %d | Bursts: %d | Mean burst size: %.2f', ...
    config.DatasetLabel, metrics.totalSpikes, metrics.numberOfBursts, ...
    metrics.meanBurstSize), 'FontWeight', 'bold');

writeAnalysisOutputs(outputDirectory, config.OutputStem, metrics, time, ...
    normalizedSignal, spikeTimes, spikeAmplitudes, interSpikeIntervals, ...
    spikesPerBurst, burstStartTimes, windowTime, windowSignal, ...
    windowSpikeMask);

saveas(analysisFigure, ...
    fullfile(outputDirectory, [config.OutputStem, '.png']));
savefig(analysisFigure, ...
    fullfile(outputDirectory, [config.OutputStem, '.fig']));
fprintf('Analysis complete. Results saved to %s\n', outputDirectory);
end


function [spikesPerBurst, burstStartTimes] = detectBursts( ...
    spikeTimes, interSpikeIntervals, maximumISI, minimumSpikes)
spikesPerBurst = [];
burstStartTimes = [];
if isempty(spikeTimes)
    return;
end

currentSpikeCount = 1;
currentStartTime = spikeTimes(1);
for index = 1:numel(interSpikeIntervals)
    if interSpikeIntervals(index) <= maximumISI
        currentSpikeCount = currentSpikeCount + 1;
    else
        if currentSpikeCount >= minimumSpikes
            spikesPerBurst(end + 1, 1) = currentSpikeCount; %#ok<AGROW>
            burstStartTimes(end + 1, 1) = currentStartTime; %#ok<AGROW>
        end
        currentSpikeCount = 1;
        currentStartTime = spikeTimes(index + 1);
    end
end
if currentSpikeCount >= minimumSpikes
    spikesPerBurst(end + 1, 1) = currentSpikeCount;
    burstStartTimes(end + 1, 1) = currentStartTime;
end
end


function writeAnalysisOutputs(outputDirectory, outputStem, metrics, time, ...
    normalizedSignal, spikeTimes, spikeAmplitudes, interSpikeIntervals, ...
    spikesPerBurst, burstStartTimes, windowTime, windowSignal, ...
    windowSpikeMask)
writetable(table(windowTime, windowSignal, ...
    'VariableNames', {'Time_s', 'Normalized_Signal'}), ...
    fullfile(outputDirectory, 'figure_a_signal.csv'));

if any(windowSpikeMask)
    writetable(table(spikeTimes(windowSpikeMask), ...
        spikeAmplitudes(windowSpikeMask), ...
        'VariableNames', {'SpikeTime_s', 'Amplitude_Normalized'}), ...
        fullfile(outputDirectory, 'figure_a_detected_spikes.csv'));
end
if ~isempty(interSpikeIntervals)
    writetable(table(interSpikeIntervals, 'VariableNames', {'ISI_s'}), ...
        fullfile(outputDirectory, 'figure_b_isi.csv'));
end
if numel(interSpikeIntervals) >= 3
    writetable(table(interSpikeIntervals(1:end-1), ...
        interSpikeIntervals(2:end), ...
        'VariableNames', {'ISI_n_s', 'ISI_n_plus_1_s'}), ...
        fullfile(outputDirectory, 'figure_c_return_map.csv'));
end
if ~isempty(spikesPerBurst)
    writetable(table((1:numel(spikesPerBurst))', spikesPerBurst, ...
        burstStartTimes, ...
        'VariableNames', {'Burst_ID', 'Spikes_in_Burst', 'Start_Time_s'}), ...
        fullfile(outputDirectory, 'figure_d_bursts.csv'));
end
if ~isempty(spikeTimes)
    writetable(table(spikeTimes, (1:numel(spikeTimes))', ...
        spikeAmplitudes, ...
        'VariableNames', {'SpikeTime_s', 'Spike_Index', ...
        'Amplitude_Normalized'}), ...
        fullfile(outputDirectory, 'detected_spikes.csv'));
end
writetable(table(time, normalizedSignal, ...
    'VariableNames', {'Time_s', 'Normalized_Signal'}), ...
    fullfile(outputDirectory, 'normalized_signal.csv'));

save(fullfile(outputDirectory, [outputStem, '_metrics.mat']), 'metrics');
reportPath = fullfile(outputDirectory, [outputStem, '_report.txt']);
fileIdentifier = fopen(reportPath, 'w');
if fileIdentifier < 0
    error('BNCAnalysis:OutputError', 'Cannot write report: %s', reportPath);
end
cleanup = onCleanup(@() fclose(fileIdentifier));
fprintf(fileIdentifier, 'Burst Analysis Report\n');
fprintf(fileIdentifier, '=====================\n');
fprintf(fileIdentifier, 'Dataset: %s\n', metrics.datasetLabel);
fprintf(fileIdentifier, 'Total spikes: %d\n', metrics.totalSpikes);
fprintf(fileIdentifier, 'Number of bursts: %d\n', metrics.numberOfBursts);
fprintf(fileIdentifier, 'Mean burst size: %.6g\n', metrics.meanBurstSize);
fprintf(fileIdentifier, 'CV of ISI: %.6g\n', metrics.cvISI);
fprintf(fileIdentifier, 'Mean firing rate (Hz): %.6g\n', ...
    metrics.meanFiringRateHz);
clear cleanup;
end


function displayUnavailable(titleText)
text(0.5, 0.5, 'Insufficient data', ...
    'HorizontalAlignment', 'center');
title(titleText); axis off;
end


function config = applyDefaults(config, defaults)
names = fieldnames(defaults);
for index = 1:numel(names)
    name = names{index};
    if ~isfield(config, name) || isempty(config.(name))
        config.(name) = defaults.(name);
    end
end
end
