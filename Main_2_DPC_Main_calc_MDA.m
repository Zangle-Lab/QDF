% last edited by Lei Tian (lei_tian@alum.mit.edu), 04/29/2015

% reference: 
% Lei Tian and Laura Waller, "Quantitative differential phase contrast
% imaging in an LED array microscope," Opt. Express 23, 11394-11403 (2015)

clear all; clc; close all;

global F Ft
% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
froot = 'Data/Cells/MDAMB231/';
savefolder = 'Data/Cells/ProcessedData/';
fstart = 'QPM10x_MDAMB231_';
numPos = 1;
numFrame =16; % 144
startframe = 7;

 

newSaveFolder = [savefolder,froot];
mkdir(newSaveFolder)
for ii = startframe:numFrame
    frame = ii;
    disp(['computing frame ', num2str(ii)])
    tic
    for jj = 1:numPos
    try
        disp([froot, 'pos', num2str(jj),'\frame',num2str(ii)])
            trickingParfor_v3([froot, 'pos', num2str(jj),'\frame',num2str(ii),'\'],newSaveFolder, fstart, jj, ii)
    end
    end

    
    toc
end