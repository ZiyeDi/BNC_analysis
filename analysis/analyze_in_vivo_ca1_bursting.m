%% Burst analysis for an in vivo CA1 extracellular recording
clear; clc; close all;

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(scriptDirectory);

dataFile = fullfile(projectRoot, 'data', 'raw', ...
    'ca1_in_vivo_recording.csv');
outputDirectory = fullfile(projectRoot, 'results', 'in_vivo_ca1');

config = struct();
config.SignalColumn = 3;
config.SignalScale = 1;
config.DetectValleys = true;
config.ThresholdFactor = 1.5;
config.EnableAdaptiveThreshold = false;
config.MaxISIWithinBurstSeconds = 0.04;
config.MinimumSpikesPerBurst = 3;
config.BurstingMeanSizeThreshold = 3.0;
config.DatasetLabel = 'In vivo CA1 extracellular signal';
config.OutputStem = 'in_vivo_ca1_burst_analysis';

metrics = run_burst_analysis(dataFile, outputDirectory, config);
disp(metrics);
