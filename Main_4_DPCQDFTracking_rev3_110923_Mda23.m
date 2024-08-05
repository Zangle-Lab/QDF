

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
        'Data/Cells/ProcessedData\MDAMB231\'; ...
        ];
% Where to find QRT Images
froot2 = [ ...
        'Data/Cells/MDAMB231\'; ...
        ];
cellType = [ 'MDAMB231'];

fstart  = ['QPM10x_MDAMB231_pos'];
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
max_disp = 60; %max displacement for particle tracking
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


 load([froot2, 'Btotal_allpos_MDA'],'Btotal');
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
        fname = strtrim([froot, fileNames(1,:)]);
%             disp(fname)
        fname2 = strtrim([frootDF,fileNames2(1,:)]);
%             disp(fname2)

        [D, ~, ~,~,~,~] = LoadSegment_CFit_mdamb231_darkfield_v2(fname,fname2,wavelength, Btotal, AverageEmptyTL, AverageEmptyTR,AverageEmptyBR, AverageEmptyBL,...
                polyfitOrder, polyfitReductionParams);
        D_old = D;
        for jj =1:numf % numf
%             try
            fname = strtrim([froot, fileNames(jj,:)]);
            disp(fname)
            fname2 = strtrim([frootDF,fileNames2(jj,:)]);
            disp(fname2)
            [D, L, L_tight,~,DF,CorrDF] = LoadSegment_CFit_mdamb231_darkfield_v2(fname,fname2,wavelength, Btotal, AverageEmptyTL, AverageEmptyTR,AverageEmptyBR, AverageEmptyBL,...
                    polyfitOrder, polyfitReductionParams);


            timen = LoadTime(fname);
            time = (datenum(timen)-datenum(time0)).*24; %store time in hours
            QDF = scalingFactor.*DF - CorrDF;
            [V, M, A, MI, P, SF] = imageprops_SF(L, D, pxlsize); %compute image properties based on the regions stored in L
            [~, M_tight, A_tight, ~, ~, SF_tight] = imageprops_SF(L_tight, D, pxlsize);
            % outputs of imageprops are in same order as labels in label
            % matrix, L. Therefore, we just need labelID to be 1:max(L(:))
            [~, ~, ~, MIDF, ~, ~] = imageprops_SF(L, DF, pxlsize);
            [~, ~, ~, MIQDF, ~, ~] = imageprops_SF(L, QDF, pxlsize);

            [~, ~, ~, MIDF_tight, ~, ~] = imageprops_SF(L_tight, DF, pxlsize);
            [~, ~, ~, MIQDF_tight, ~, ~] = imageprops_SF(L_tight, QDF, pxlsize);
            if std(D(:))~=0 %skip if blank image
                [yshift, xshift] = CorrShift(D_old,D); %find average shift between current frame and first frame
                yshift = yshift+yshift_old;
                xshift = xshift+xshift_old;
                D_old = D;
                xshift_old = xshift;
                yshift_old = yshift;

                yshift_store(jj,Loc) = yshift;
                xshift_store(jj,Loc) = xshift;
                % save D and L into D_stored and L_stored
                D_stored(:,:,jj) = int16(D(:,:)/4/pi*65536/wavelength);
                L_stored(:,:,jj) = uint16(L(:,:));
                t_stored(jj,Loc) = time;
                DF_stored(:,:,jj) =int16(DF);
                CorrDF_stored(:,:,jj) = int16(CorrDF);
                %next, loop through all items identified in V and find only the ones
                %which meet area and mean intensity requirements
                for ii = 1:length(V)
                    %first, check that 1) there is something at index ii, 2) that
                    if(~isnan(P(ii).Centroid(1)) && A(ii) > minAreathresh && A(ii) < maxAreathresh && MI(ii) > minMIthresh && MI(ii) < maxMIthresh)
                        T_array(tt,1:2) = P(ii).Centroid; %store position in first two columns of T_array
                        T_array(tt,1:2) = T_array(tt,1:2) - [xshift, yshift]; %remove shift due to movement of the entire frame
                        T_array(tt,3)   = M(ii);          %store mass in third column
                        T_array(tt,4)   = time;           %store time from first frame in seconds in 4th column
                        T_array(tt,5)   = A(ii); %store area in fifth column
                        T_array(tt,6)   = SF(ii); %store shape factor in sixth column
                        T_array(tt,7)   = MIDF(ii);
                        T_array(tt,8)   = MIQDF(ii);
                        % Stores the new tighter Label variables
                        T_array(tt,9)   = A_tight(ii);
                        T_array(tt,10)   = M_tight(ii);
                        T_array(tt,11)   = MIDF_tight(ii);
                        T_array(tt,12)   = MIQDF_tight(ii);
                        T_array(tt,13)   = SF_tight(ii);
                        labelID(tt) = ii;
                        tt = tt+1;                        %increment T_array index
                    end
                end
            end
%             end
            %clear D B L V M A MI %erase data from this step
        end
        segID = [segID; labelID];

        if ~isempty(T_array) && sum(T_array(:,4) ~= T_array(1,4))>0
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Cell tracking using the track function
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %T_array starts as [x, y, m, t, A, SF]
            T_array(:,3) = T_array(:,3).*massfact; %change mass weighting in T_array
            T_array = sortrows(T_array, 4); %sort T_array based on time vectors
            minTx =  min(T_array(:,1));
            T_array(:,1) = T_array(:,1) -minTx +1; %make sure all x positions are positive\
            minTy =  min(T_array(:,2));
            T_array(:,2) = T_array(:,2) -minTy +1; %make sure all y positions are positive, new with rev6
            %move time to last column, T_array will now be [x, y, m, A, SF, t]
%             T_array_temp = [T_array(:,1:3), T_array(:,5:8), segID', T_array(:,4)];
            T_array_temp = [T_array(:,1:3), T_array(:,5:13), segID', T_array(:,4)];

            tracks = track(T_array_temp,max_disp,param);
            tracks0 = tracks;
            
            %move area back to 5th column, T_array will now be [x, y, m, t, A, SF] and
            %tracks will now be [x, y, m, t, cellnum, A, SF]
            % originally tracks comes out of the track function with 7
            % columns, I inserted a 6th column into T_array (the input for
            % the track function) so now it comes out with 8 columns

            % the way the track function seems to work is it takes in an
            % array (T_array) and then adds a column onto the end, which is
            % the cellnum. This is what the tracking does, it labels each
            % row with a cellnum so and then reorganizes the rows so that
            % the cellnum is constant in order by frame until the next cell
            % num comes up?

%             tracks_temp = [tracks(:,1:3), tracks(:,6:7), tracks(:,4:5)];
            % Temp Tracks for no
%             tracks_temp = [tracks(:,1:3), tracks(:,9:10), tracks(:,4:7)];
            tracks_temp = [tracks(:,1:3), tracks(:,14:15), tracks(:,4:12)];
            sortedSeg = tracks(:,13);
            tracks = tracks_temp;

            T_array(:,3) = T_array(:,3)./massfact; %adjust mass weighting back to the way it was
            T_array(:,1) = T_array(:,1) + minTx -1; %set all x positions back to the way they were
            T_array(:,1) = T_array(:,2) + minTy -1; %set all y positions back to the way they were
            tracks(:,3) = tracks(:,3)./massfact; %adjust for mass weighting in the tracks array as well
            tracks(:,1) = tracks(:,1) +minTx -1; %set all x positions back to the way they were
            tracks(:,2) = tracks(:,2) +minTy -1; %set all y positions back to the way they were
        else
            % initialize tracks and sortedSeg in case condition not met
            % above
            tracks = [];
            sortedSeg = [];
            tracks_temp = [];
        end
        %%%%%%%%%This is the only addition from previous analysis%%%%%%%%%%
        % for each location, insert each location number
        tracks(:,8) = Loc;
        times = tracks(:,4);

        % if you turn off parfor, the way you access t_stored might be an
        % error (something like t_stored(:,2)
        % finds the unique time for each frame, should have a length equal
        % to the number of frames in the experiment
        uniqueTimes = t_stored(:,Loc);

        % numf is also the length of uniqueTimes, so for each unique time:
        for ii = 1:numf
            % associate frame number with timestamp on the image
            currentTime = tracks(:,4) == uniqueTimes(ii);
            % insert frame number
            tracks(currentTime,9) = ii;
        end

        % find drug ID and concentration of drug based on its location on
        % the 96 well plate
        % TODO: need to define plate structure out here, and send it to 
        % function, NEED to save plate structure in data_allframes
%         if(panel == 1)
%             [drugID, conc, ~, ~, ~] = getCondition_PanelA_HCI27BS(Loc);
%         elseif(panel == 2)
%             [drugID, conc, ~, ~, ~] = getCondition_PanelB_HCI27BS(Loc);
%         end
        drugID = 0;
        conc   = 0;

        sizeTracks = size(tracks);
        numRows = sizeTracks(1);
        % check to make sure that there is something in segID 
        if isempty(segID) || numRows == 1
            % if empty there should only be one row in tracks, set 10
            % element to zero
            tracks(:,10) = 0;
            tracks(:,11) = 0;
            tracks(:,12) = 0;
            % Store MIDF and MIQDF
            tracks(:,13) = 0;
            tracks(:,14) = 0;
            % Sores Tight variables (Area, Mass, DF, QDF,SF)
            tracks(:,15) = 0;
            tracks(:,16) = 0;
            tracks(:,17) = 0;
            tracks(:,18) = 0;
            tracks(:,19) = 0;
        else
            % if it isn't empty then the size of these arrays should match5
            tracks(:,10) = sortedSeg;
    
            tracks(:,11) = drugID;
            tracks(:,12) = conc;
            % Store MIDF and MIQDF
            tracks(:,13) = tracks_temp(:,8);
            tracks(:,14) = tracks_temp(:,9);
            % Sores Tight variables (Area, Mass, DF, QDF,SF)
            tracks(:,15) = tracks_temp(:,10);
            tracks(:,16) = tracks_temp(:,11);
            tracks(:,17) = tracks_temp(:,12);
            tracks(:,18) = tracks_temp(:,13);
            tracks(:,19) = tracks_temp(:,14);
        end

        %tracks will now be [x, y, m, t, cellnum, A, SF]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %save D_stored and L_stored in a separate .mat file to save memory

        %save(['data', row, col, '.mat'], 'D_stored', 'L_stored') %
        % line above is an alternative to the two lines below (I think the
        % parsave() is for saving in parallel to speed up processing.
        parsave([savefolder 'datav3_', num2str(Loc), '.mat'], D_stored, 'D_stored', 0)
        parsave([savefolder 'datav3_', num2str(Loc), '.mat'], L_stored, 'L_stored', 1)
        parsave([savefolder 'datav3_', num2str(Loc), '.mat'], DF_stored, 'DF_stored', 1)
        parsave([savefolder 'datav3_', num2str(Loc), '.mat'], CorrDF_stored, 'CorrDF_stored', 1)
        %clear D_stored L_stored

        tracksColHeaders = {'X', 'Y', 'Mass (pg)', 'Time (h)', 'id', 'Area', 'shape factor', ...
                            'Location ID', 'Frame ID', 'segmentLabel', 'condition_drugID', ...
                            'condition_concentration (um)','Darkfield Mean Intensity',...
                            'QDF Mean Intensity (au)','Area tight','Mass tight (pg)','Darkfield Mean Intensity tight (au)',...
                            'QDF Mean Intensity tight (au)','Shape Factor tight'};


        %save tracks data and clear vector so that code can be parallelized
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], tracks, 'tracks', 0)
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], T_array, 'T_array', 1)
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], xshift_store, 'xshift_store', 1)
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], yshift_store, 'yshift_store', 1)
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], t_stored, 't_stored', 1)
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], tracksColHeaders, 'tracksColHeaders', 1)
        %clear T_array xshift_store yshift_store t_stored row_stored col_stored

    end %close if statement looking for stored data
end %close rows for loop
% 
% clear tracks T_array xshift_store yshift_store t_stored D_stored L_stored T_array_temp tracks_temp
% toc
% disp('End of parfor loop')
% tic
% datetime('now')
% % matlabpool close
% %%
% tic
% %reconstitute tracks vector
% tracks_comp = [];
% xshift_store_c = [];
% yshift_store_c = [];
% t_stored_c = [];
% ii_stored = [];
% maxindex = 0;
% Loc_stored_c = [];
% for Loc = 46:63
% % for Loc = 61
% %     Loc = 1;
%     % change 54 to Loc in the next line for multiple positions
%     load([savefolder, 'Tdata', num2str(Loc), '.mat']);
%     if ~isempty(tracks)
%         tracks(:,5) = tracks(:,5)+maxindex;
%     end
% 
%     tracks_comp = [tracks_comp; tracks];
%     xshift_store_c = [xshift_store_c; xshift_store(:,Loc)];
%     yshift_store_c = [yshift_store_c; yshift_store(:,Loc)];
%     t_stored_c = [t_stored_c; t_stored(:,Loc)];
%     ii_stored = [ii_stored, 1:length(t_stored(:,Loc))'];
%     Loc_stored_c = [Loc_stored_c, (1:length(t_stored(:,Loc))').*0 + Loc];
% 
%     if ~isempty(tracks_comp)
%         maxindex = max(tracks_comp(:,5));
%     end
% 
%     clear tracks xshift_store yshift_store t_stored
% 
% end
% 
% Loc_stored = Loc_stored_c;
% tracks = tracks_comp;
% xshift_store = xshift_store_c;
% yshift_store= yshift_store_c;
% t_stored = t_stored_c;
% ii_stored = ii_stored';
% clear tracks_comp xshift_store_c yshift_store_c t_stored_c maxindex Loc_stored_c
% if length(tracks) > 10
%     T0=min(tracks(:,4)); %find time of first image in the set
%     tracks(:,4) = (tracks(:,4)-T0);
%     t_stored = t_stored - T0;
% end
% toc
% %%
% % tic
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%Code to extract growth rates
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % [num, indices] = track_numpart(tracks,minpathlength); %find tracks >= minpathlength
% % for ii = 1:num
% %     currentnum = tracks(indices(ii),5);
% %     [x,y,z,t] = track_partn_SF(tracks,currentnum);
% %     BF= polyfit(t, z, 1); %find best fit line to the data
% %     GRate(ii) = BF(1); %store growth rate
% %     Mass0(ii) = z(1); %store initial masss
% %     MassF(ii) = z(end); %store final mass
% %     Length(ii) = length(x); % store length of track?
% % end
% % 
% % toc
% %%
% % figure(2)
% % % figure
% % %plot growth over time
% % minpathlength = 30;
% % [num, indices] = track_numpart(tracks,minpathlength); %find tracks >= minpathlength
% % 
% % %plot those tracks
% % for ii = 1:num
% %     currentnum = tracks(indices(ii),5);
% %     [x,y,z,t] = track_partn_SF(tracks,currentnum);
% %     plot(t, z, '.-')
% %     hold on
% %     Pf = polyfit(t, z,1);
% %     text(t(end)*1.02, z(end), [num2str(currentnum), ', ', num2str(z(1), '%0.2f'),', ', num2str(Pf(1), '%0.4f'), ', ', num2str(Pf(2))])
% %     plot([t(1) t(end)], [Pf(1)*t(1) Pf(1)*t(end)]+Pf(2), 'g')
% %     %text(t(end)*1.02, z(end)*0.99, [num2str(Pf(1), '%0.4f'), ', ', num2str(Pf(2))])
% %     %text(t(end)*1.02, z(end)*0.98, num2str(z(1)))
% % end
% % hold off
% % 
% % ylabel('Mass (pg)', 'FontSize', 14)
% % xlabel('time (h)', 'FontSize', 14)
% % set(gca, 'FontSize', 14)
% % title('labels: "cell #, initial mass, best fit slope, best fit intercept"', 'FontSize', 10)
% % 
% % 
% % %find and store data from this plot
% % h = findobj(gca,'Type','line');
% % fig2x=get(h,'Xdata');
% % fig2y=get(h,'Ydata');
% 
% %%
% % figure(3)
% % %plot normalized growth over time
% % [num, indices] = track_numpart(tracks,minpathlength); %find tracks >= minpathlength
% % 
% % for ii = 1:num
% %     currentnum = tracks(indices(ii),5);
% %     [x,y,z,t] = track_partn_SF(tracks,currentnum);
% %     z_norm = z./z(1); %normalize by first mass values
% %     plot(t, z_norm, '.-')
% %     hold on
% % end
% % hold off
% % ylabel('Normalized Mass', 'FontSize', 14)
% % xlabel('time (h)', 'FontSize', 14)
% % set(gca, 'FontSize', 14)
% 
% 
% %%
% %%% added these two lines to remove superfluous images and save this data%%%
% tic
% clear D D1 D1s D_old L L1 Btotal
% if(panel == 1)
%     [drugID, conc, ~, ~, ~] = getCondition_PanelA_HCI27BS(Loc);
% elseif(panel == 2)
%     [drugID, conc, ~, ~, ~] = getCondition_PanelB_HCI27BS(Loc);
% end 

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %save data
% if ~exist([savefolder 'data_allframes']) || overwrite
%     save([savefolder 'data_allframes'])
% end
% 
% % I modified this line to work with new data output
% % f3_AnalyzeTData_v3(savefolder)
% toc
% disp('data_allframes is created')
% datetime('now')
% 
% %%
% % load([savefolder 'data1'])
% %% Check mass/area/meanintensity of cells in a particular frame
% % 
% % fignum = 2;
% % figure(fignum)
% % % loc is position/location number
% % % loc =54;
% % % nn is the frame number
% % for Loc = [46, 55, 82, 91, 118, 127, 154, 163]
% % 
% % for nn =10
% % 
% %     [LocList, numLoc] = getloclist_v2(froot, fstart, fext, date);
% %     %pre-processing and variable initialization before loop begins:
% %     filelist = dir([froot, fstart, char(LocList(Loc)), '_*', fext]);
% %     fileNames = char(sort_nat({filelist.name}'));
% %     filelist2 = dir([froot2,'pos','*']);
% %     foldersOnly = [filelist2.isdir];
% %     fileNames2 = char(sort_nat({filelist2(foldersOnly).name}'));
% % %     if numLoc~=length(fileNames2)
% % %         error('Unequal number of phase and darkfield positions')
% % %     end
% %     posNames = fileNames2;
% %     fname = strtrim([froot, fileNames(nn,:)]);
% % 
% %     load([froot, 'Btotal'], 'Btotal');
% % %     Btotal = B3D;
% %     polyfitReductionParams = [10, 10];
% %     % order is the "order of polynomial fit"
% %     polyfitOrder = 8;
% %     % L1 = L_stored(:,:,nn);
% %     % [D1, ~] = LoadSegment_Btotal_short(fname, wavelength, Btotal);
% %             [D1,L1, ~,~] =LoadSegment_CFit_QDF_060222(fname,wavelength, Btotal, ...
% %                     polyfitOrder, polyfitReductionParams);                     % [D1] = LoadSegment(fname,wavelength);
% %     time0 = LoadTime(fname);
% % 
% % %     L1(190:220,635:670) = L1(252,637);
% % 
% %     D1s = zeros([size(D1),numLoc], 'single');
% % %     se90 = strel('line', 5, 90);
% % %     se0 = strel('line', 5, 0);
% % %     L1 = imdilate(L1,[se90 se0]);
% %     scale = [-300 1200];
% %     % D1F = imfilter(abs(D_stored(:,:,nn)), fspecial('gaussian', [10 10], 2));
% %     % D1F = medfilt2(D_stored(:,:,nn));
% %     % D1F = D_stored(:,:,nn);
% %     % L2 = imagesegment_aggressive(D1F);
% %     % L2 = L_stored(:,:,nn);
% %     % [V1, M1, A1, MI1, P1, SF1] = imageprops_SF(L_stored(:,:,nn), D1F, pxlsize); %compute image properties
% %     [V1, M1, A1, MI1, P1, SF1] = imageprops_SF(L1, D1, pxlsize); %compute image properties
% %     % imagesc(D_stored(:,:,nn))
% %     % plotBWoutlines(D_stored(:,:,nn),L_stored(:,:,nn),12);
% %     plotBWoutlines(D1, L1,fignum);
% %     caxis(scale)
% %     colormap parula
% %     axis image
% %     title(['areas in ',num2str(nn),'th frame that are counted as "cells" labeled with masses'])
% %     lowA1 = find(A1 > minAreathresh & A1 < maxAreathresh & MI1 > minMIthresh & MI1 < maxMIthresh); %find indices of objects which are counted as "cells"
% %     % lowA1 = find(A1); %uncomment this to select all the cells in the image
% % 
% %     for k = 1 : length(lowA1)
% %         rectangle('Position', P1(lowA1(k)).BoundingBox, 'EdgeColor','y'); %plot yellow box around each cell
% % 
% %         %comment out the next line if you want to label with something other than mass
% %         text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(M1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame    text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(M1(lowA1(k))), 'Color', [1 1 0]) %label each cell with its mass in the first frame
% % 
% %         %uncomment the next line if you want to label with area instead of mass
% %     %     text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(A1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame
% % 
% %         %uncomment the next line if you want to label with mean intensity instead of mass
% %     %     text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(MI1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame
% % 
% %         %uncomment the next line if you want to label with shape factor instead of mass
% %         %text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(SF1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame
% % 
% %     end
% %     hold off;
% %     pause(0.5)
% % end
% % end
% % %%
% % 
% % fignum =4;
% % figure(fignum)
% % % nn is the frame number
% % nn =5;
% % scale = [0 300];
% % % D1F = imfilter(abs(D_stored(:,:,nn)), fspecial('gaussian', [10 10], 2));
% % % D1F = double(D_stored(:,:,nn))*2*pi./65536*wavelength;
% % % D1F = double(DF_stored(:,:,nn));
% % 
% % % D1F = D_stored(:,:,nn);
% % % L2 = imagesegment_aggressive(D1F);
% % L2 = L_stored(:,:,nn);
% % % [V1, M1, A1, MI1, P1, SF1] = imageprops_SF(L_stored(:,:,nn), D1F, pxlsize); %compute image properties
% % % [V1, M1, A1, MI1, P1, SF1] = imageprops_SF(L2, D1F, pxlsize); %compute image properties
% % % imagesc(D_stored(:,:,nn))
% % % plotBWoutlines(D_stored(:,:,nn),L_stored(:,:,nn),12);
% % plotBWoutlines(0.8*double(DF_stored(:,:,nn))-double(CorrDF_stored(:,:,nn)),L2,fignum);
% % % plotBWoutlines(double(DF_stored(:,:,nn)),L2,fignum);
% % 
% % caxis(scale)
% % colormap gray
% % % title('areas in nth frame that are counted as "cells" labeled with masses')
% % % lowA1 = find(A1 > minAreathresh & A1 < maxAreathresh & MI1 > minMIthresh & MI1 < maxMIthresh); %find indices of objects which are counted as "cells"
% % % % lowA1 = find(A1); %uncomment this to selectfdytg all the cells in the image
% % % 
% % % for k = 1 : length(lowA1)
% % %     rectangle('Position', P1(lowA1(k)).BoundingBox, 'EdgeColor','y'); %plot yellow box around each cell
% % %     
% % %     %comment out the next line if you want to label with something other than mass
% % %     text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(M1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame    text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(M1(lowA1(k))), 'Color', [1 1 0]) %label each cell with its mass in the first frame
% % %     
% % %     %uncomment the next line if you want to label with area instead of mass
% % % %     text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(A1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame
% % %     
% % %     %uncomment the next line if you want to label with mean intensity instead of mass
% % % %     text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(MI1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame
% % %     
% % %     %uncomment the next line if you want to label with shape factor instead of mass
% % %     %text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(SF1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame
% % %     
% % % end
% % % hold off;
% % 
% % % toc
% % % %% %%
% % % figure(13)
% % % for nn = 1:1
% % % %     D1F = imfilter(abs(D_stored(:,:,nn)), fspecial('gaussian', [10 10], 2));
% % %     D1F = D_stored(:,:,nn);
% % % %     L2 = imagesegment_aggressive(D1F);
% % %     L2 = L_stored(:,:,nn);    
% % %     scale = [-300 1200];
% % %     % [V1, M1, A1, MI1, P1, SF1] = imageprops_SF(L_stored(:,:,nn), D1F, pxlsize); %compute image properties
% % %     [V1, M1, A1, MI1, P1, SF1] = imageprops_SF(L2, D1F, pxlsize); %compute image properties
% % %     %imagesc(D_stored(:,:,nn))
% % %     % plotBWoutlines(D_stored(:,:,nn),L_stored(:,:,nn),12);
% % %     plotBWoutlines(D_stored(:,:,nn), scale,L2,13);
% % %     colormap parula
% % %     title('areas in nth frame that are counted as "cells" labeled with masses')
% % %     lowA1 = find(A1 > minAreathresh & A1 < maxAreathresh & MI1 > minMIthresh & MI1 < maxMIthresh); %find indices of objects which are counted as "cells"
% % % %     lowA1 = find(A1); %uncomment this to select all the cells in the image
% % % 
% % %     for k = 1 : length(lowA1)
% % %         rectangle('Position', P1(lowA1(k)).BoundingBox, 'EdgeColor','y'); %plot yellow box around each cell
% % % 
% % %         %comment out the next line if you want to label with something other than mass
% % %             text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(M1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame    text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(M1(lowA1(k))), 'Color', [1 1 0]) %label each cell with its mass in the first frame
% % % 
% % %         %uncomment the next line if you want to label with area instead of mass
% % % %         text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(A1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame
% % % 
% % %         %uncomment the next line if you want to label with mean intensity instead of mass
% % %     %     text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(MI1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame
% % % 
% % %         %uncomment the next line if you want to label with shape factor instead of mass
% % %         %text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(SF1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame
% % % 
% % %     end
% % %     hold off;
% % %     pause(.1);
% % % end
% % %% 
