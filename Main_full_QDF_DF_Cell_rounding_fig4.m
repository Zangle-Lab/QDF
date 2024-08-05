clearvars
close all
load('Data/Cells/FullDataSet/data_allframes_MDA.mat')

%% Specify min path length, cell ID location and average min mass
minp          = 5; % min number of frames a cell has to be considered a cell
cellsloc      = 5; % location of cell ID array in the tracks
minmass       = 150; %pg % minimum mass to be considered a cell
pxlsize       = 5.3619e-1; %um
%% Find IDs that are above or equal  minp
cellsID       = unique(tracks(:,cellsloc));
Ncount        = histc(tracks(:,cellsloc),cellsID); 
targetcells   = cellsID(Ncount>= minp); 
NewCount      = Ncount(Ncount>= minp);
%% Check min average mass, QDF and DF of each cell

remove= zeros(1,length(targetcells)); % Empty array to track how many cells were removed
mass  =  tracks(:,3);
DF    = tracks(:,13);
QDF   = tracks(:,14);
IDs   = tracks(:,5);
Area  = tracks(:,6);
CorrDF = 0.8*tracks(:,13) -tracks(:,14);
% SF    = tracks(:,7);
% TightSF= tracks(:,19);
% can be parallelized with parfor
for i = 1:length(targetcells)
    y = IDs==targetcells(i);
    m = mean(mass(y));
    minDF = min(DF(y));
    minQDF = min(QDF(y));
    maxArea = max(Area(y));
    minCorrDF = min(CorrDF(y));
    % requirement to remove false data and bad segmentations
    if m <= minmass || minDF<0 || minQDF <0 || maxArea >8000 || minCorrDF <0
        remove(i) = 1;
    end
end

remove = remove >0;
%%

targetcells(remove) = [];
NewCount(remove)    = [];
%% Reconstruct tracks based on valid cells

newtracks = nan(sum(NewCount),19);
newtracks(1:NewCount(1),:) = tracks(tracks(:,5) ==targetcells(1),:);

for z = 1:(length(targetcells)-1)
    NewIndex = z +1;
    startingRow = sum(NewCount(1:(NewIndex-1)));
    endingRow   = sum(NewCount(1:NewIndex));
    newtracks((startingRow+1):endingRow,:)  = tracks(tracks(:,5) ==targetcells(NewIndex),:);

end
newtracks(:,6) = newtracks(:,6).*(pxlsize.^2);
%% variables for fitting
cleanMass = newtracks(:,3);
cleanDF   = newtracks(:,13).*newtracks(:,6);
cleanQDF  = newtracks(:,14).*newtracks(:,6);
cleanArea = newtracks(:,6);
cleanMeanMass =cleanMass./cleanArea;
cleanDeltaDF = cleanDF -cleanQDF;
cleanCorrDF = 0.8*cleanDF - cleanQDF;
cleanCorrMeanDF = cleanCorrDF./cleanArea;
cleanDFPerMass = cleanDF./cleanMass;
cleanQDFPerMass = cleanQDF./cleanMass;
cleanCorrPerMass = cleanCorrDF./cleanMass;
cleanSF    = newtracks(:,7);
cleanMeanDF = newtracks(:,13);
cleanMeanQDF = newtracks(:,14);



%% Creating figure

%% Fig 1 MeanMass (x-axis) vs total DF and total QDF
% Darkfield
[xData1a, yData1a] = prepareCurveData( cleanMeanMass, cleanDF );

% Set up fittype and options.
ft1a = fittype( 'poly1' );

% Fit model to data.
[fitresult1a, gof1a] = fit( xData1a, yData1a, ft1a );

% Plot fit with data.
f1 = figure(1);
f1.WindowState = 'maximized';

% h = plot( fitresult1a, xData1a, yData1a );
scatter(xData1a,yData1a,'MarkerEdgeColor','none',...
              'MarkerFaceColor',[0 0.4470 0.7410])
hold on
% legend( h, 'cleanDF vs. cleanMeanMass', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Mass per Area (pg/\mum^{2})', 'Interpreter', 'tex' );
% ylabel( 'cleanDF', 'Interpreter', 'none' );
% grid on

% QDF
[xData1b, yData1b] = prepareCurveData( cleanMeanMass, cleanQDF );

% Set up fittype and options.
ft1b = fittype( 'poly1' );

% Fit model to data.
[fitresult1b, gof1b] = fit( xData1b, yData1b, ft1b );

% Plot fit with data.
scatter(xData1b,yData1b,'MarkerEdgeColor','none',...
              'MarkerFaceColor',[0.9290 0.6940 0.1250])
hold off
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult1b, xData1b, yData1b );
% legend( h, 'cleanQDF vs. cleanMeanMass', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
% xlabel( 'cleanMeanMass', 'Interpreter', 'none' );
ylabel( 'Darkfield signal (au)', 'Interpreter', 'none' );
% grid on


%% figure 2 Mass per Area vs DF/QDF per mass

[xData2a, yData2a] = prepareCurveData( cleanMeanMass, cleanDFPerMass );

% Set up fittype and options.
ft2a = fittype( 'poly1' );

% Fit model to data.
[fitresult2a, gof2a] = fit( xData2a, yData2a, ft2a );

% Plot fit with data.
f2 =figure(2);

    f2.WindowState = 'maximized';
scatter(xData2a,yData2a,'MarkerEdgeColor','none',...
              'MarkerFaceColor',[0 0.4470 0.7410])
hold on

xlabel( 'Mass per Area (pg/\mum^{2})', 'Interpreter', 'tex' );

% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult2a, xData2a, yData2a );
% legend( h, 'cleanDFPerMass vs. cleanMeanMass', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
% xlabel( 'cleanMeanMass', 'Interpreter', 'none' );
% ylabel( 'cleanDFPerMass', 'Interpreter', 'none' );
% grid on

[xData2b, yData2b] = prepareCurveData( cleanMeanMass, cleanQDFPerMass );

% Set up fittype and options.
ft2b = fittype( 'poly1' );

% Fit model to data.
[fitresult2b, gof2b] = fit( xData2b, yData2b, ft2b );
scatter(xData2b,yData2b,'MarkerEdgeColor','none',...
              'MarkerFaceColor',[0.9290 0.6940 0.1250])
hold off

ylabel( 'Darkfield signal per mass (/pg)', 'Interpreter', 'tex' );

% Plot fit with data.
% figure(2);
% h = plot( fitresult2a, xData2a, yData2a );
% legend( h, 'cleanQDFPerMass vs. cleanMeanMass', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
% xlabel( 'cleanMeanMass', 'Interpreter', 'none' );
% % ylabel( 'cleanQDFPerMass', 'Interpreter', 'none' );
% grid on
%     saveas(f2, ['E:\Work\MdaMB231_NegativeQDFExperiment\DPCImages\MassPerAreaVsDFPerMass.svg'])


%% Figure 3 mass per area vs DF/QDF per area

[xData3a, yData3a] = prepareCurveData( cleanMeanMass, cleanMeanDF );

% Set up fittype and options.
ft3a = fittype( 'poly1' );

% Fit model to data.
[fitresult3a, gof3a] = fit( xData3a, yData3a, ft3a );
f3 = figure(3);
    f3.WindowState = 'maximized';

scatter(xData3a,yData3a,'MarkerEdgeColor','none',...
              'MarkerFaceColor',[0 0.4470 0.7410])
hold on

xlabel( 'Mass per Area (pg/\mum^{2})', 'Interpreter', 'tex' );

[xData3b, yData3b] = prepareCurveData( cleanMeanMass, cleanMeanQDF );

% Set up fittype and options.
ft3b = fittype( 'poly1' );

% Fit model to data.
[fitresult3b, gof3b] = fit( xData3b, yData3b, ft3b );


scatter(xData3b,yData3b,'MarkerEdgeColor','none',...
              'MarkerFaceColor',[0.9290 0.6940 0.1250])
hold off

ylabel( 'Darkfield signal per area (/\mum^{2})', 'Interpreter', 'tex' );

R_coefficient = [sqrt(gof1a.rsquare) sqrt(gof1b.rsquare);...
                sqrt(gof2a.rsquare) sqrt(gof2b.rsquare);...
                sqrt(gof3a.rsquare) sqrt(gof3b.rsquare);...
                ]

