%% Burst analysis for a BNC hardware output signal
clear; clc; close all;

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(scriptDirectory);

dataFile = fullfile(projectRoot, 'data', 'raw', ...
    'bnc_hardware_bursting.xlsx');
outputDirectory = fullfile(projectRoot, 'results', 'bnc_hardware');

config = struct();
config.SignalColumn = 2;
config.SignalScale = 1;
config.DetectValleys = false;
config.ThresholdFactor = 1.9;
config.EnableAdaptiveThreshold = true;
config.AdaptiveThresholdFactors = [2.0, 1.5, 1.2, 1.0, 0.8, 0.6];
config.MinimumAdaptiveSpikes = 50;
config.MaximumAdaptiveSpikes = 1000;
config.MaxISIWithinBurstSeconds = 0.04;
config.MinimumSpikesPerBurst = 3;
config.BurstingMeanSizeThreshold = 3.0;
config.DatasetLabel = 'BNC hardware output signal';
config.OutputStem = 'bnc_hardware_burst_analysis';

metrics = run_burst_analysis(dataFile, outputDirectory, config);
disp(metrics);
