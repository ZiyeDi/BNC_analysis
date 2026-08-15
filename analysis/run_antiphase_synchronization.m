%% Anti-phase synchronization analysis for paired BNC signals
clear; clc; close all;

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(scriptDirectory);

dataFile = fullfile(projectRoot, 'data', 'raw', ...
    'phase_locking_antiphase.xlsx');
outputDirectory = fullfile(projectRoot, 'results', ...
    'phase_locking_antiphase');

options = struct();
options.FrequencyRange = 10:2:200;
options.RelativeBandwidth = 0.20;
options.PhaseBins = 100;
options.LoopLabel = 'Mutual inhibition';
options.ExpectedPhaseRadians = pi;
options.OutputStem = 'antiphase_synchronization';

results = compute_phase_frequency_map(dataFile, outputDirectory, options);
disp(results);
