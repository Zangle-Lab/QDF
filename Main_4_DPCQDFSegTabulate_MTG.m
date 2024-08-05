

% script to run cell tracking codel-ko
% TAZ 5/18/10

clear all; close all;
tic
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set up file names of images to be analyzed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

%%%%%%%%%%removed anything to do with wells replaced these 2 lines%%%%%%%%%
% Where to find Phase images and commonBackground
froot = [ ...
        'Data/Cells/ProcessedData\MTG021_MTG084\'; ...
        ];
% Where to find QRT Images
froot2 = [ ...
        'Data/Cells/MTG021_MTG084\'; ...
        ];
cellType = [ 'MTG'];

fstart  = ['QPM10x_MTG_pos'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavelength = 624; %nm
pxlsize = 5.3619e-4; %mm/pixel, for 10x
savefolder = froot;

%%% First, define which files the script will work on
fext = '.mat'; %file extension
overwrite = 1; %set to 1 to enable overwrite of pre-stored data files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second, define image analysis parameters

%%% define min and max area and mean intensity (MI) thresholds
%%% only objects which fall between these values will be counted as "cells"
%%% These parameters should be adjusted for each sample to capture the
%%% objects of interest. See Figure 11 to evaluate where these values fall
%%% relative to the properties of the image. See figures 12 and 13 to see
%%% which objects in the first and last frames are counted as "cells"
% minAreathresh = 500;
% 10um dataset

minAreathresh = 500; 
maxAreathresh = 100000; 
minMIthresh = 80;
maxMIthresh = 800;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Third, define tracking software parameters
scalingFactor   = .8; %This is a scaling factor for DF vs QDF due to them being unequal in signal
minpathlength = 20; %min path length to use in plotting results. only paths
%                    of this length or longer will be displayed. this does
%                    not affect the tracking software (tracks shorter than
%                    minpathlength will still be computed and stored)
% .5 went well5


%%% tracking parameters below affect the tracking software itself
max_disp = 120; %max displacement for particle tracking
%               max_disp is an estimate of the maximum distance that a
%               particle would move in a single time interval. It should be
%               set to a value somewhat less than the mean spacing between
%               the particles
massfact = 1;   %factor to multiply mass by in tracking step. use this to
%account for differences in how much the cell moves vs. how
%much mass changes over time (was set to 2.5 previously, 2/25/15)
param.mem = 0;  %this is the number of time steps that a particle can be
%               'lost' and then recovered again.  If the particle reappears
%               after this number of frames has elapsed, it will be
%               tracked as a new particle. The default setting is zero.
%               this is useful if particles occasionally 'drop out' of
%               the data.
param.dim = 3; %number of dimensions of coordinate data to track. If you
%               set this value to 2, then it will track based on position
%               in x and y. If param.dim = 3, then the software will track
%               based on mass as well.
param.good = 0; %set this keyword to eliminate all trajectories with
%                fewer than param.good valid positions.  This is useful
%                for eliminating very short, mostly 'lost' trajectories
%                due to blinking 'noise' particles in the data stream.
param.quiet = 1; %set this keyword to 1 if you don't want any text
%                 displayed while the tracking algorithm is running


%%

%%


[LocList, numLoc] = getloclist_v2(froot, fstart, fext, cellType);
%pre-processing and variable initialization before loop begins:
Loc = 1;
filelist = dir([froot, fstart, char(LocList(Loc)), '_*', fext]);
fileNames = char(sort_nat({filelist.name}'));
filelist2 = dir([froot2,'pos','*']);
foldersOnly = [filelist2.isdir];
fileNamesDF = char(sort_nat({filelist2(foldersOnly).name}'));

posNames = fileNamesDF;

%%

numf = length(fileNames);

fname = strtrim([froot, fileNames(1,:)]);

% reductionParams defines the size reduction used for the polyfit
% in order to decrease the time required to compute polyfit surface
% use caution when changing these parameters
polyfitReductionParams = [10, 10];
% order is the "order of polynomial fit"
polyfitOrder = 8;


 load([froot2, 'Btotal_allpos_MTG'],'Btotal');
load(['ProcessedData\AverageEmpty.mat']);

load(fname,'Phase')

time0 = LoadTime(fname);


clearvars D1
%%
for Loc = 1:numLoc %parfor                                                                                                                                                                                  
    if ~exist([savefolder, 'data_', num2str(Loc), '.mat'],'file') || overwrite
        
        filelist = dir([froot, fstart, char(LocList(Loc)), '_*', fext]);
        fileNames = char(sort_nat({filelist.name}'));
        filelist2 = dir([froot2,strtrim(fileNamesDF(Loc,:)),'\darkfield\frame','*']);
        foldersOnly = [filelist2.isdir];
        fileNames2 = char(sort_nat({filelist2(foldersOnly).name}'));
    

        numf = length(fileNames(:,1)); % track all frames at position
        CFit = cell([1,numf]);
        %%%%%%%%%%%%%tfnh%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% grab first frame for analysis and detection of the correct cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fname = strtrim([froot, fileNames(1,:)]);

        polyfitParams = struct('polyfitOrder', {polyfitOrder}, ...
                               'polyfitReductionParams', {polyfitReductionParams}, ...
                               'CFit', {CFit});
        load(fname,'Phase')

        %preallocate variables for speed
        yshift_store = zeros(numf);
        xshift_store = zeros(numf);
        t_stored = zeros(numf);
        % these can be commented out to save memory, if you do that then
        % make sure to comment out the lines that use them belows
        D_stored = zeros([size(Phase),numf], 'int16');
        L_stored = zeros([size(Phase),numf], 'uint16');
        DF_stored = zeros([size(Phase),numf], 'int16');
        CorrDF_stored = zeros([size(Phase),numf], 'int16');


        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% loop through first numf file names stored in fnum and store analysis
        %%% results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tt = 1; %initialize tt, the index of the tracking array
        T_array = [];

        xshift_old = 0;
        yshift_old = 0;
        segID = [];
        labelID = [];

        FROOT2 = strtrim([froot2,'pos',LocList{Loc}(3)]);
        frootDF = [FROOT2,'\darkfield\'];
        for jj =1:numf % numf
%             try
            fname = strtrim([froot, fileNames(jj,:)]);
            disp(fname)
            fname2 = strtrim([frootDF,fileNames2(jj,:)]);
            disp(fname2)
            [D, L, L_tight,~,DF,CorrDF] = LoadSegment_CFit_mdamb231_darkfield(fname,fname2,wavelength, Btotal, AverageEmptyTL, AverageEmptyTR,AverageEmptyBR, AverageEmptyBL,...
                    polyfitOrder, polyfitReductionParams);

            timen = LoadTime(fname);
            time = (datenum(timen)-datenum(time0)).*24; %store time in hours
            QDF = scalingFactor.*DF - CorrDF;
            [V, M, A, MI, P, SF] = imageprops_SF(L, D, pxlsize); %compute image properties based on the regions stored in L
            
            % outputs of imageprops are in same order as labels in label
            % matrix, L. Therefore, we just need labelID to be 1:max(L(:))
            [~, ~, ~, MIDF, ~, ~] = imageprops_SF(L, DF, pxlsize);
            [~, ~, AQDF, MIQDF, ~, ~] = imageprops_SF(L, QDF, pxlsize);

%             [~, ~, ~, MIDF_tight, ~, ~] = imageprops_SF(L_tight, DF, pxlsize);


                % save D and L into D_stored and L_stored
                D_stored(:,:,jj) = int16(D(:,:)/4/pi*65536/wavelength);
                L_stored(:,:,jj) = uint16(L(:,:));
                t_stored(jj,Loc) = time;
                DF_stored(:,:,jj) =int16(DF);
                CorrDF_stored(:,:,jj) = int16(CorrDF);
                for ii = 1:length(V)
                        T_array(tt,1)  = AQDF(ii); %Area of puncta
                        T_array(tt,2)  = MIQDF(ii); %Mean intensity of puncta
                    tt = tt +1;
                end
        end
    end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %save D_stored and L_stored in a separate .mat file
        parsave([savefolder 'data', num2str(Loc), '.mat'], D_stored, 'D_stored', 0)
        parsave([savefolder 'data', num2str(Loc), '.mat'], L_stored, 'L_stored', 1)
        parsave([savefolder 'data', num2str(Loc), '.mat'], DF_stored, 'DF_stored', 1)
        parsave([savefolder 'data', num2str(Loc), '.mat'], CorrDF_stored, 'CorrDF_stored', 1)
        parsave([savefolder 'Tabledata', num2str(Loc), '.mat'], T_array, 'T_array', 0)

end 
%%
clearvars -except savefolder numLoc
T_arrayComplete = [];
for Loc = 1:numLoc
    load([savefolder 'Tabledata', num2str(Loc), '.mat'])
    T_arrayComplete = [T_arrayComplete;T_array];
end

clearvars -except savefolder T_arrayComplete

T_array = T_arrayComplete;
parsave([savefolder 'TabledataAll.mat'], T_array, 'T_array', 0)
