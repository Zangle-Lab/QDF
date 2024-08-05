clear all
load('Data/Cells/ProcessedData\MTG021_MTG084\Tabledata1.mat')
TotalQDF21 = T_array(:,1).*T_array(:,2);
load('Data/Cells/ProcessedData\MTG021_MTG084\Tabledata2.mat')
TotalQDF21 = [TotalQDF21;T_array(:,1).*T_array(:,2)];
load('Data/Cells/ProcessedData\MTG021_MTG084\Tabledata3.mat')
TotalQDF84 = [T_array(:,1).*T_array(:,2)];

%%

%%


TotalQDF = [TotalQDF21;TotalQDF84];
grp = [21*ones(size(TotalQDF21)); 84.*ones(size(TotalQDF84))];

violinplot(TotalQDF,grp);
title('Total QDF')
fh.WindowState= 'maximized';
