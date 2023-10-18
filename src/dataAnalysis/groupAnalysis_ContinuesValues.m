clc; clear all; close all;
paths = getFilesPaths();
%files  = getAllFiles();%getAllFiles, getPreFiles, getPostFiles
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
