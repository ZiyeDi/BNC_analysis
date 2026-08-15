%% In-phase synchronization analysis for paired BNC signals
clear; clc; close all;

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(scriptDirectory);

dataFile = fullfile(projectRoot, 'data', 'raw', ...
    'phase_locking_inphase.xlsx');
outputDirectory = fullfile(projectRoot, 'results', ...
    'phase_locking_inphase');

options = struct();
options.FrequencyRange = 2:2:200;
options.RelativeBandwidth = 0.05;
options.PhaseBins = 100;
options.LoopLabel = 'Mutual excitation';
options.ExpectedPhaseRadians = 0;
options.OutputStem = 'inphase_synchronization';

results = compute_phase_frequency_map(dataFile, outputDirectory, options);
disp(results);
