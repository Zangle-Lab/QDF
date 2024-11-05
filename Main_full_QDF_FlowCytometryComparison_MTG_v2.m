clear all
load('Data/Cells/FullDataSet/Tabledata_allframes_MTG.mat')

%%


%%

y = T_array(:,3) >350 & T_array(:,7) >0 & T_array (:,8) >0 & T_array(:,5) <20000 ;

%% MTG21 variables
AllowedMTG21Pos= [46,49,50,51,52,53,54,55,60,61,62];
% filtered position based focus stability, debris, and having cells
x = ismember(T_array(:,13),AllowedMTG21Pos);

QDF21 = T_array(x&y,8);
DF21 = T_array(x&y,7);
Area21 = T_array(x&y,5);
Mass21 = T_array(x&y,3);
MeanMass21 = Mass21./Area21;
TQDF21 = QDF21.*Area21;
TDF21   = DF21.*Area21;
QDFPerMass21 = TQDF21./Mass21;
DFPerMass21  = TDF21./Mass21;
AreaPuncta21 = T_array(x&y,11);
QDFPunctaperArea21 = T_array(x&y,12);

%% MTG84 variables
AllowedMTG84Pos = [154,156,157,159,160,162,164,1655,166,168,169,171];
z = ismember(T_array(:,13),AllowedMTG84Pos);
z = z & T_array(:,5)>2000;
QDF84 = T_array(z&y,8);
DF84 = T_array(z&y,7);
Area84 = T_array(z&y,5);
Mass84 = T_array(z&y,3);
MeanMass84 = Mass84./Area84;
TQDF84 = QDF84.*Area84;
TDF84   = DF84.*Area84;

QDFPerMass84 = TQDF84./Mass84;
DFPerMass84  = TDF84./Mass84;
AreaPuncta84 = T_array(z&y,11);
QDFPunctaperArea84 = T_array(z&y,12);



%%
fh = figure(5111);
Areacomp = [Area21; Area84];
TotalQDF21 = QDF21.*Area21;
TotalQDF84 = QDF84.*Area84;
TotalDF21 = DF21.*Area21;
TotalDF84 = DF84.*Area84;


TotalDF = [TotalDF21;TotalDF84];
TotalQDF = [TotalQDF21;TotalQDF84];
grp = [21*ones(size(TotalQDF21)); 84.*ones(size(TotalQDF84))];

I = grp == 84;
ZQDF = relativeEntropy(TotalQDF,I)
ZDF = relativeEntropy(TotalDF,I)

violinplot(TotalQDF,grp);
title('Total QDF')
% ylim([-1e5 13e5])
fh.WindowState= 'maximized';

fh2 = figure(5112);
TotalDF = [TotalDF21;TotalDF84];
grp = [21*ones(size(TotalDF21)); 84.*ones(size(TotalDF84))];
violinplot(TotalDF,grp);
fh2.WindowState= 'maximized';