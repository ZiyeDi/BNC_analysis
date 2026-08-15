%% Selected energy-envelope feature overlap analysis
clear; clc; close all;

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
biologicalDataFile = fullfile(projectRoot, 'data', 'raw', ...
    'biological_segments.xlsx');
hardwareDataFile = fullfile(projectRoot, 'data', 'raw', ...
    'hardware_segments.xlsx');
outputDirectory = fullfile(projectRoot, 'results', 'feature_overlap');

if ~isfile(biologicalDataFile)
    error('BNCAnalysis:MissingInput', ...
        'Biological input file not found: %s', biologicalDataFile);
end
if ~isfile(hardwareDataFile)
    error('BNCAnalysis:MissingInput', ...
        'Hardware input file not found: %s', hardwareDataFile);
end
if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end

fprintf('Loading biological and hardware signal segments...\n');
biologicalData = readmatrix(biologicalDataFile);
hardwareData = readmatrix(hardwareDataFile);

biologicalData = biologicalData(all(isfinite(biologicalData(:, 1:2)), 2), :);
hardwareData = hardwareData(all(isfinite(hardwareData(:, 1:2)), 2), :);

biologicalTime = biologicalData(:, 1);
hardwareTime = hardwareData(:, 1);
biologicalSamplingRate = 1 / median(diff(biologicalTime));
hardwareSamplingRate = 1 / median(diff(hardwareTime));

numberOfSamples = min([50, size(biologicalData, 2) - 1, ...
    size(hardwareData, 2) - 1]);
if numberOfSamples < 2
    error('BNCAnalysis:InsufficientData', ...
        'At least two signal columns are required in each input file.');
end

biologicalSignals = biologicalData(:, 2:numberOfSamples + 1);
hardwareSignalsOriginal = hardwareData(:, 2:numberOfSamples + 1);
hardwareSignals = zeros(size(biologicalSignals));

fprintf('Resampling hardware signals from %.3f Hz to %.3f Hz...\n', ...
    hardwareSamplingRate, biologicalSamplingRate);
for sampleIndex = 1:numberOfSamples
    hardwareSignals(:, sampleIndex) = interp1( ...
        hardwareTime, hardwareSignalsOriginal(:, sampleIndex), ...
        biologicalTime, 'spline');
end

samplingRate = biologicalSamplingRate;
upperFrequency = min(150, 0.95 * samplingRate / 2);
if upperFrequency <= 1
    error('BNCAnalysis:SamplingRateTooLow', ...
        'The sampling rate is too low for the configured bandpass filter.');
end

[filterB, filterA] = butter(4, ...
    [1, upperFrequency] / (samplingRate / 2), 'bandpass');
for sampleIndex = 1:numberOfSamples
    biologicalSignals(:, sampleIndex) = filtfilt( ...
        filterB, filterA, biologicalSignals(:, sampleIndex));
    hardwareSignals(:, sampleIndex) = filtfilt( ...
        filterB, filterA, hardwareSignals(:, sampleIndex));
end

biologicalSignalsNormalized = zscore(biologicalSignals, 0, 1);
hardwareSignalsNormalized = zscore(hardwareSignals, 0, 1);
biologicalSignalsNormalized(~isfinite(biologicalSignalsNormalized)) = 0;
hardwareSignalsNormalized(~isfinite(hardwareSignalsNormalized)) = 0;

featureNames = {'IntegratedEnvelopeSkewness', ...
    'EnvelopePSDSlope', 'BurstFanoFactor'};
biologicalFeatures = zeros(numberOfSamples, numel(featureNames));
hardwareFeatures = zeros(numberOfSamples, numel(featureNames));

fprintf('Extracting energy-envelope features...\n');
for sampleIndex = 1:numberOfSamples
    biologicalFeatures(sampleIndex, :) = extractEnvelopeFeatures( ...
        biologicalSignalsNormalized(:, sampleIndex), samplingRate);
    hardwareFeatures(sampleIndex, :) = extractEnvelopeFeatures( ...
        hardwareSignalsNormalized(:, sampleIndex), samplingRate);
end

biologicalFeatures(~isfinite(biologicalFeatures)) = 0;
hardwareFeatures(~isfinite(hardwareFeatures)) = 0;
features = [biologicalFeatures; hardwareFeatures];
labels = [ones(numberOfSamples, 1); zeros(numberOfSamples, 1)];
groupLabels = [repmat({'Biological'}, numberOfSamples, 1); ...
    repmat({'Hardware'}, numberOfSamples, 1)];

standardizedFeatures = zscore(features, 0, 1);
standardizedFeatures(~isfinite(standardizedFeatures)) = 0;

fprintf('Training a cross-validated linear SVM...\n');
svmModel = fitcsvm(features, labels, ...
    'Standardize', true, 'KernelFunction', 'linear');
crossValidatedModel = crossval(svmModel, 'KFold', 10);
[~, predictionScores] = kfoldPredict(crossValidatedModel);
positiveClassIndex = find(svmModel.ClassNames == 1, 1);
[falsePositiveRate, truePositiveRate, ~, areaUnderCurve] = ...
    perfcurve(labels, predictionScores(:, positiveClassIndex), 1);

classificationFigure = figure('Name', 'Feature Space and Classification', ...
    'Color', 'w', 'Position', [100, 100, 1100, 480]);
subplot(1, 2, 1);
scatter3(standardizedFeatures(labels == 1, 1), ...
    standardizedFeatures(labels == 1, 2), ...
    standardizedFeatures(labels == 1, 3), 40, 'b', 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6); hold on;
scatter3(standardizedFeatures(labels == 0, 1), ...
    standardizedFeatures(labels == 0, 2), ...
    standardizedFeatures(labels == 0, 3), 40, 'r', '^', 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6);
xlabel('Standardized Skewness');
ylabel('Standardized PSD Slope');
zlabel('Standardized Fano Factor');
title('A. Standardized Feature Space');
legend('Biological', 'Hardware', 'Location', 'northeast');
view(35, 25); grid on; axis tight;

subplot(1, 2, 2);
plot(falsePositiveRate, truePositiveRate, 'LineWidth', 3, ...
    'Color', [0.85, 0.33, 0.10]); hold on;
plot([0, 1], [0, 1], 'k--', 'LineWidth', 1.5);
xlabel('False Positive Rate'); ylabel('True Positive Rate');
title('B. Linear SVM Classification');
text(0.5, 0.3, sprintf('AUC = %.3f', areaUnderCurve), ...
    'FontSize', 14, 'FontWeight', 'bold', ...
    'Color', [0.85, 0.33, 0.10]);
legend('Cross-validated ROC', 'Chance reference', 'Location', 'southeast');
grid on;

pValues = zeros(1, numel(featureNames));
for featureIndex = 1:numel(featureNames)
    pValues(featureIndex) = ranksum( ...
        biologicalFeatures(:, featureIndex), ...
        hardwareFeatures(:, featureIndex));
end

distributionFigure = figure('Name', 'Feature Distributions', ...
    'Color', 'w', 'Position', [150, 150, 1200, 450]);
plotLabels = {'Integrated Envelope Skewness', ...
    'Envelope PSD Slope', 'Burst Fano Factor'};
plotColors = {[0.49, 0.18, 0.56], [0.00, 0.45, 0.74], ...
    [0.47, 0.67, 0.19]};

for featureIndex = 1:numel(featureNames)
    subplot(1, 3, featureIndex);
    featureValues = features(:, featureIndex);
    categoricalGroups = categorical(groupLabels);
    boxchart(categoricalGroups, featureValues, ...
        'BoxFaceColor', plotColors{featureIndex}, 'MarkerStyle', 'none');
    hold on;
    swarmchart(categoricalGroups, featureValues, 20, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.4, 'XJitterWidth', 0.5);
    ylabel(plotLabels{featureIndex});
    title(sprintf('%s\nP = %.4g', plotLabels{featureIndex}, ...
        pValues(featureIndex)));
    grid on;
end

saveas(classificationFigure, ...
    fullfile(outputDirectory, 'feature_space_and_roc.png'));
savefig(classificationFigure, ...
    fullfile(outputDirectory, 'feature_space_and_roc.fig'));
saveas(distributionFigure, ...
    fullfile(outputDirectory, 'feature_distributions.png'));
savefig(distributionFigure, ...
    fullfile(outputDirectory, 'feature_distributions.fig'));

rocTable = table(falsePositiveRate, truePositiveRate, ...
    'VariableNames', {'FalsePositiveRate', 'TruePositiveRate'});
writetable(rocTable, fullfile(outputDirectory, 'roc_curve.csv'));

standardizedFeatureTable = table(categorical(groupLabels), ...
    standardizedFeatures(:, 1), standardizedFeatures(:, 2), ...
    standardizedFeatures(:, 3), ...
    'VariableNames', {'Group', 'StandardizedSkewness', ...
    'StandardizedPSDSlope', 'StandardizedFanoFactor'});
writetable(standardizedFeatureTable, ...
    fullfile(outputDirectory, 'standardized_features.csv'));

rawVariableNames = [{'Group'}, featureNames];
rawFeatureTable = table(categorical(groupLabels), ...
    features(:, 1), features(:, 2), features(:, 3), ...
    'VariableNames', rawVariableNames);
writetable(rawFeatureTable, ...
    fullfile(outputDirectory, 'raw_features.csv'));

biologicalVariableNames = cellstr( ...
    ["Time_s", "Biological_" + string(1:numberOfSamples)]);
hardwareVariableNames = cellstr( ...
    ["Time_s", "Hardware_" + string(1:numberOfSamples)]);
biologicalNormalizedTable = array2table( ...
    [biologicalTime, biologicalSignalsNormalized], ...
    'VariableNames', biologicalVariableNames);
hardwareNormalizedTable = array2table( ...
    [biologicalTime, hardwareSignalsNormalized], ...
    'VariableNames', hardwareVariableNames);
writetable(biologicalNormalizedTable, ...
    fullfile(outputDirectory, 'biological_signals_normalized.xlsx'));
writetable(hardwareNormalizedTable, ...
    fullfile(outputDirectory, 'hardware_signals_normalized.xlsx'));

results = struct();
results.numberOfSamplesPerGroup = numberOfSamples;
results.biologicalSamplingRateHz = biologicalSamplingRate;
results.hardwareSamplingRateHz = hardwareSamplingRate;
results.featureNames = featureNames;
results.biologicalFeatures = biologicalFeatures;
results.hardwareFeatures = hardwareFeatures;
results.pValues = pValues;
results.linearSVMAUC = areaUnderCurve;
save(fullfile(outputDirectory, 'feature_overlap_results.mat'), 'results');

fprintf('Feature analysis complete. Results saved to %s\n', outputDirectory);


function features = extractEnvelopeFeatures(signal, samplingRate)
features = zeros(1, 3);
rawEnvelope = abs(hilbert(signal));

membraneWindow = max(1, round(samplingRate * 0.02));
integratedEnvelope = smoothdata(rawEnvelope, 'gaussian', membraneWindow);
features(1) = skewness(integratedEnvelope);

macroWindow = max(1, round(samplingRate * 0.05));
macroEnvelope = smoothdata(rawEnvelope, 'gaussian', macroWindow);
zeroMeanEnvelope = macroEnvelope - mean(macroEnvelope);
[powerSpectrum, frequencies] = pwelch( ...
    zeroMeanEnvelope, hamming(numel(zeroMeanEnvelope)), [], ...
    numel(zeroMeanEnvelope), samplingRate);

fitMask = frequencies >= 2 & frequencies <= 18 & powerSpectrum > 0;
if nnz(fitMask) > 3
    coefficients = robustLineFit(log10(frequencies(fitMask)), ...
        log10(powerSpectrum(fitMask)));
    features(2) = coefficients(2);
else
    features(2) = 0;
end

burstThreshold = prctile(macroEnvelope, 80);
minimumPeakDistance = max(1, round(samplingRate * 0.05));
[~, peakLocations] = findpeaks(macroEnvelope, ...
    'MinPeakHeight', burstThreshold, ...
    'MinPeakDistance', minimumPeakDistance);

fanoWindow = max(1, round(0.25 * samplingRate));
numberOfWindows = floor(numel(macroEnvelope) / fanoWindow);
eventCounts = zeros(numberOfWindows, 1);
for windowIndex = 1:numberOfWindows
    windowStart = (windowIndex - 1) * fanoWindow + 1;
    windowEnd = windowIndex * fanoWindow;
    eventCounts(windowIndex) = sum( ...
        peakLocations >= windowStart & peakLocations <= windowEnd);
end

if ~isempty(eventCounts) && mean(eventCounts) > 0
    features(3) = var(eventCounts) / mean(eventCounts);
else
    features(3) = 0;
end
end


function coefficients = robustLineFit(predictor, response)
predictor = predictor(:);
response = response(:);
designMatrix = [ones(size(predictor)), predictor];
coefficients = designMatrix \ response;

for iteration = 1:200
    residuals = response - designMatrix * coefficients;
    robustScale = median(abs(residuals - median(residuals))) / 0.6745;
    if ~isfinite(robustScale) || robustScale <= eps
        break;
    end

    scaledResiduals = residuals / (4.685 * robustScale);
    weights = zeros(size(scaledResiduals));
    included = abs(scaledResiduals) < 1;
    weights(included) = (1 - scaledResiduals(included) .^ 2) .^ 2;
    if nnz(weights) < 2
        break;
    end

    weightedDesign = designMatrix .* sqrt(weights);
    weightedResponse = response .* sqrt(weights);
    updatedCoefficients = weightedDesign \ weightedResponse;
    if norm(updatedCoefficients - coefficients) <= ...
            1e-8 * max(1, norm(coefficients))
        coefficients = updatedCoefficients;
        break;
    end
    coefficients = updatedCoefficients;
end
end
