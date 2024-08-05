% function [ D, L, B ] = LoadSegment_Btotal_short_rollingball( fname, wavelength, Btotal )
function [ D2, L, L_tight, M,DF,CorrDF] = LoadSegment_CFit_mdamb231_darkfield( fname, fname2, wavelength, Btotal, AverageEmptyTL, AverageEmptyTR,AverageEmptyBR, AverageEmptyBL,order, polyfitReductionParams)

%function [ D, L, B ] = LoadSegment( fname, wavelength )
%function to load and segment the data stored in fname
%phase data should be stored in fname.Phase

% Btotal version loads Btotal to subtract it from the loaded image
%   This is usually used when running the tracking code

%TM 10/12/2020 Edited the load line to only load phase and use it directly
%This to save memory needed by half. Combined the steps to make D into 1
%step

load(fname,'Phase');
if isa(Phase,'single')
    D = double(Phase).*wavelength;
elseif isa(Phase,'int16')
    D = double(Phase)*(4*pi).*wavelength/65536;
else
    error('Error. \nInput must be a int16 or single, not a %s.',class(Phase))
end
B = regexp(strtrim(fname2(end-3:end)),'[0-9]','match');
Num = char(1);
for ii= 1:length(B)
      Num(ii)=B{ii};

end
        frame = str2double(Num);
        TL = double(imread([fname2, '\QrtImage_',num2str(frame),'TL.tiff']));
        TR = double(imread([fname2, '\QrtImage_',num2str(frame),'TR.tiff']));
        BL = double(imread([fname2, '\QrtImage_',num2str(frame),'BL.tiff']));
        BR = double(imread([fname2, '\QrtImage_',num2str(frame),'BR.tiff']));
  
        TL = floor(TL./16)-AverageEmptyTL;
        TR = floor(TR./16)-AverageEmptyTR;
        BR = floor(BR./16)-AverageEmptyBR;
        BL = floor(BL./16)-AverageEmptyBL;
% polyfitReductionParams defines the size reduction used for the polyfit
% in order to decrease the time required to compute polyfit surface
[B1, M, DF,CorrDF] = imagebackground_polyn_reduced2_DF_v1(D-Btotal,TL,TR,BR,BL,order,polyfitReductionParams(1),polyfitReductionParams(2));    %compute image background
D1 = D-B1-Btotal;           %subtract background from data files
D2 = D1 - mean(D1(M)); % reset background to zero

DF1 = medfilt2(D2); % median filter background corrected phase image
[L,L_tight] = imagesegment_aggressive_MDA_DF(DF1); %segment image (detect distinct cell regions and disconnect connected regions)

end