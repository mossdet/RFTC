clear all; close all; clc;

paths = getFilesPaths();
workspace = paths.workspacePath;
run('convertMicromedToMatlab.m');
run('savePatientData_RFTC.m');
run('detectAndCharacterize.m');
run('artefactCorrectedAvgChannelFeatures.m');
run('completeDataForZoneAssignment.m')


% Paper
run('generatePatientsTable.m');
run('biomarkerZones_AnalysisZones_Correlation.m');
run('plot_Zone_Zone_Correlation_Circles.m');
run('plot_biomarkerZones_AnalysisZones_Correlation.m');
run('plotSignificantDecreaseCircles.m');
run('newGroupAnalysis_PrePost_ImprNonImp.m');

%Prediction
run('pdfsFeaturesZones_DiffOutcomes_HFO.m');
run('newGroupAnalysis_PrePost_ImprNonImp.m');
run('getPredictionZoneCharacterization.m'); % predict outcome with k-means
run('getAverageEventsFromPatsAndZones.m');

