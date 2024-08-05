clear, close all


% Copyright Tarek Moustafa (Zangle lab, University of Utah) 2024
% Code to commute Empty background correct for DF Images
%% File location and location + Scaling
tic 
froot_logs = 'Data/Cells/EmptyDarkfieldImages/';
froot = 'Data/Cells/ProcessedData/';
numPos_QDF = [1:9];
numFrames_QDF = [1];

AverageEmptyTL = zeros(1200,1920);
AverageEmptyTR = zeros(1200,1920);
AverageEmptyBL = zeros(1200,1920);
AverageEmptyBR = zeros(1200,1920);

c =0; % counter

%%
for pp = numPos_QDF
    tic
    display(['Computing pos ',num2str(pp)])

for ff = numFrames_QDF
    
            FRoot = [froot_logs,'pos',num2str(pp),'/frame',num2str(ff),'/'];
            TL = floor(double(imread([FRoot, 'QrtImage_',num2str(ff),'TL.tiff'])))./16;
            TR = floor(double(imread([FRoot, 'QrtImage_',num2str(ff),'TR.tiff'])))./16;
            BL = floor(double(imread([FRoot, 'QrtImage_',num2str(ff),'BL.tiff'])))./16;
            BR = floor(double(imread([FRoot, 'QrtImage_',num2str(ff),'BR.tiff'])))./16;

            AverageEmptyTL = (AverageEmptyTL.*c + TL)./(c+1);
            AverageEmptyTR = (AverageEmptyTR.*c + TR)./(c+1);
            AverageEmptyBL = (AverageEmptyBL.*c + BL)./(c+1);
            AverageEmptyBR = (AverageEmptyBR.*c + BR)./(c+1);

            c=c+1; 
            
end
end

save([froot, 'AverageEmpty'],'AverageEmptyTL','AverageEmptyTR','AverageEmptyBL','AverageEmptyBR')
disp('Empty Darkfield is saved')
toc
