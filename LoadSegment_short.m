% function [ D, L, B ] = LoadSegment( fname, wavelength, Btotal )
function [ D, L] = LoadSegment_short( fname, wavelength )
%function to load and segment the data stored in fname
%phase data should be stored in fname.Phase

% wavelength = 660;
% filename = 'data\QPM40X_pos54_frame_92.mat';
% load(filename);
% I = Phase;
% D = I.*wavelength;

% Load segment loads the current frame and subtracts polynomial fit
   % This is usually used when running the FindCommonBack_Phasics
   
% ERP modified 8/26/20 to work with a new data type. Previously, this code
% loaded images that were single precision, however now we are loading
% images that are short integers (16-bit signed integers). Therefore, this
% code needs to convert it from being in the range [-65536,-65536] to the
% range [-2pi, 2pi]. To do this we just divide out 2^16 and multiply by 2pi

% TM modified 10/15/20 (1) Made it only load Phase matrix not the whole file to
% save ram space (2) Changed the reduction from 50 to 1 to remove
% randmization

load(fname,'Phase');
floating = double(Phase)*(2*pi)/65536;
D = floating.*wavelength;

DF = medfilt2(D);
L = imagesegment_aggressive_v4(DF);      %segment image (detect distinct cell regions and disconnect connected regions)


end

