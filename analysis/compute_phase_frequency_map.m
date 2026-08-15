function results = compute_phase_frequency_map(dataFile, outputDir, options)
%COMPUTE_PHASE_FREQUENCY_MAP Analyze phase locking across frequency.
%   RESULTS = COMPUTE_PHASE_FREQUENCY_MAP(DATAFILE, OUTPUTDIR, OPTIONS)
%   reads a three-column matrix [time, signal_1, signal_2], estimates
%   phase coherence and phase-binned cross-power, and saves figures and a
%   MAT file in OUTPUTDIR.

if nargin < 3
    options = struct();
end

options = applyDefaults(options, struct( ...
    'FrequencyRange', 2:2:200, ...
    'RelativeBandwidth', 0.05, ...
    'PhaseBins', 100, ...
    'LoopLabel', 'Paired signals', ...
    'ExpectedPhaseRadians', 0, ...
    'OutputStem', 'phase_frequency'));

if ~isfile(dataFile)
    error('BNCAnalysis:MissingInput', 'Input file not found: %s', dataFile);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

data = readmatrix(dataFile);
if size(data, 2) < 3
    error('BNCAnalysis:InvalidInput', ...
        'Expected at least three columns: time, signal_1, signal_2.');
end

data = data(:, 1:3);
data = data(all(isfinite(data), 2), :);
if size(data, 1) < 10
    error('BNCAnalysis:InsufficientData', ...
        'At least ten finite samples are required.');
end

time = data(:, 1);
signal1 = data(:, 2) - mean(data(:, 2));
signal2 = data(:, 3) - mean(data(:, 3));
samplingRate = 1 / median(diff(time));

frequencyRange = options.FrequencyRange(:)';
phaseEdges = linspace(-pi, pi, options.PhaseBins + 1);
phaseCenters = (phaseEdges(1:end-1) + phaseEdges(2:end)) / 2;
powerMap = zeros(numel(frequencyRange), options.PhaseBins);
phaseCoherence = nan(numel(frequencyRange), 1);

fprintf('Computing phase-frequency map for %s...\n', options.LoopLabel);

for index = 1:numel(frequencyRange)
    centerFrequency = frequencyRange(index);
    bandwidth = centerFrequency * options.RelativeBandwidth;
    lowFrequency = max(0.1, centerFrequency - bandwidth);
    highFrequency = min(centerFrequency + bandwidth, 0.95 * samplingRate / 2);

    if lowFrequency >= highFrequency
        continue;
    end

    filtered1 = bandpass(signal1, [lowFrequency, highFrequency], samplingRate);
    filtered2 = bandpass(signal2, [lowFrequency, highFrequency], samplingRate);
    analytic1 = hilbert(filtered1);
    analytic2 = hilbert(filtered2);

    amplitude1 = abs(analytic1);
    amplitude2 = abs(analytic2);
    phaseDifference = angle(exp(1i * (angle(analytic1) - angle(analytic2))));
    phaseCoherence(index) = abs(mean(exp(1i * phaseDifference)));
    crossPower = amplitude1 .* amplitude2;

    [~, ~, phaseIndex] = histcounts(phaseDifference, phaseEdges);
    valid = phaseIndex > 0 & isfinite(crossPower);
    if any(valid)
        meanPower = accumarray(phaseIndex(valid), crossPower(valid), ...
            [options.PhaseBins, 1], @mean, 0);
        powerMap(index, :) = meanPower;
    end
end

maximumPower = max(powerMap, [], 'all');
if maximumPower > 0
    normalizedPowerMap = powerMap / maximumPower;
else
    normalizedPowerMap = powerMap;
end

modulationStrength = max(normalizedPowerMap, [], 2) - ...
    min(normalizedPowerMap, [], 2);
modulationStrengthSmoothed = movmean(modulationStrength, 3);
phaseCoherenceSmoothed = movmean(phaseCoherence, 3, 'omitnan');

validFit = isfinite(phaseCoherenceSmoothed);
if nnz(validFit) >= 4
    polynomialCoefficients = polyfit( ...
        frequencyRange(validFit), phaseCoherenceSmoothed(validFit)', 3);
    phaseCoherenceFit = polyval(polynomialCoefficients, frequencyRange);
    residuals = phaseCoherenceSmoothed(validFit) - ...
        phaseCoherenceFit(validFit)';
    residualSumSquares = sum(residuals .^ 2);
    totalSumSquares = sum((phaseCoherenceSmoothed(validFit) - ...
        mean(phaseCoherenceSmoothed(validFit))) .^ 2);
    if totalSumSquares > 0
        fitR2 = 1 - residualSumSquares / totalSumSquares;
    else
        fitR2 = NaN;
    end
else
    polynomialCoefficients = nan(1, 4);
    phaseCoherenceFit = nan(size(frequencyRange));
    fitR2 = NaN;
end

[peakPower, peakFrequencyIndex] = max(max(powerMap, [], 2));
bestFrequency = frequencyRange(peakFrequencyIndex);
bestPhaseDistribution = powerMap(peakFrequencyIndex, :);
if peakPower > 0
    bestPhaseDistribution = bestPhaseDistribution / peakPower;
end

[maximumCoherence, maximumCoherenceIndex] = max(phaseCoherenceSmoothed);

mainFigure = figure('Name', 'Phase-Frequency Analysis', 'Color', 'w', ...
    'Position', [100, 100, 1400, 700]);

subplot(2, 3, [1, 2]);
plot(time, signal1, 'b', 'LineWidth', 0.5); hold on;
plot(time, signal2, 'r', 'LineWidth', 0.5);
title(sprintf('A. Raw Signals (%s)', options.LoopLabel));
xlabel('Time (s)'); ylabel('Amplitude'); axis tight; grid on;
legend('Signal 1', 'Signal 2', 'Location', 'best');

subplot(2, 3, 3);
imagesc(phaseCenters, frequencyRange, normalizedPowerMap);
set(gca, 'YDir', 'normal'); colormap('turbo'); colorbar;
xlabel('Phase Difference (rad)'); ylabel('Frequency (Hz)');
title('B. Mean Power Correlation Map');
xline(options.ExpectedPhaseRadians, 'w--', 'LineWidth', 1.5, ...
    'Label', 'Expected locking');

subplot(2, 3, 4);
closedPhase = [phaseCenters, phaseCenters(1)];
closedPower = [bestPhaseDistribution, bestPhaseDistribution(1)];
polarplot(closedPhase, closedPower, 'r-', 'LineWidth', 2);
rlim([0, 1.1]);
title(sprintf('C. Phase Tuning at %.1f Hz', bestFrequency));

subplot(2, 3, [5, 6]);
plot(frequencyRange, phaseCoherenceSmoothed, 'b-', 'LineWidth', 2, ...
    'DisplayName', 'Phase coherence'); hold on;
plot(frequencyRange, phaseCoherenceFit, 'r--', 'LineWidth', 2, ...
    'DisplayName', sprintf('Cubic fit (R^2 = %.3f)', fitR2));
plot(frequencyRange(maximumCoherenceIndex), maximumCoherence, 'ro', ...
    'MarkerSize', 9, 'MarkerFaceColor', 'r', ...
    'DisplayName', sprintf('Peak at %.1f Hz', ...
    frequencyRange(maximumCoherenceIndex)));
xlabel('Frequency (Hz)'); ylabel('Phase Coherence');
title('D. Frequency-Dependent Phase Coherence');
grid on; xlim([min(frequencyRange), max(frequencyRange)]);
legend('Location', 'best');

secondaryFigure = figure('Name', 'Phase Modulation Summary', 'Color', 'w', ...
    'Position', [100, 100, 900, 400]);
subplot(1, 2, 1);
plot(frequencyRange, modulationStrengthSmoothed, 'g-', 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Phase Modulation Strength'); grid on;
title('Phase Modulation Strength');
subplot(1, 2, 2);
scatter(phaseCoherenceSmoothed, modulationStrengthSmoothed, 30, ...
    frequencyRange, 'filled');
xlabel('Phase Coherence'); ylabel('Phase Modulation Strength');
title('Coherence versus Modulation'); colorbar; grid on;

saveas(mainFigure, fullfile(outputDir, [options.OutputStem, '_map.png']));
savefig(mainFigure, fullfile(outputDir, [options.OutputStem, '_map.fig']));
saveas(secondaryFigure, ...
    fullfile(outputDir, [options.OutputStem, '_summary.png']));

results = struct();
results.dataFile = dataFile;
results.samplingRateHz = samplingRate;
results.frequencyRangeHz = frequencyRange;
results.phaseCentersRadians = phaseCenters;
results.powerMap = powerMap;
results.normalizedPowerMap = normalizedPowerMap;
results.phaseCoherence = phaseCoherence;
results.phaseCoherenceSmoothed = phaseCoherenceSmoothed;
results.phaseCoherenceFit = phaseCoherenceFit;
results.polynomialCoefficients = polynomialCoefficients;
results.fitR2 = fitR2;
results.bestFrequencyHz = bestFrequency;
results.maximumCoherence = maximumCoherence;
save(fullfile(outputDir, [options.OutputStem, '_results.mat']), 'results');

fprintf('Analysis complete. Results saved to %s\n', outputDir);
end


function options = applyDefaults(options, defaults)
names = fieldnames(defaults);
for index = 1:numel(names)
    name = names{index};
    if ~isfield(options, name) || isempty(options.(name))
        options.(name) = defaults.(name);
    end
end
end
