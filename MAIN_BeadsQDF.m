clear, close all


% Copyright Tarek Moustafa (Zangle lab, University of Utah) 2024
%% load files
BL = double(imread('Data/Beads/QrtDF_BL.tiff'));
BR = double(imread('Data/Beads/QrtDF_BR.tiff'));
TL = double(imread('Data/Beads/QrtDF_TL.tiff'));
TR = double(imread('Data/Beads/QrtDF_TR.tiff'));


%% Scale to 12 bit
BL = ceil(BL./16);
BR = ceil(BR./16);
TL = ceil(TL./16);
TR = ceil(TR./16);

%%

pxlsize=5.316e-4; %mm


%% Make a mask for background
DF_Mask = TL + BR +BL +TR;
mask = DF_Mask>= 2e4/16;
maskfilled = imfill(mask, 'holes');

%% Correct Backgrounnd for all

BL = BL - mean(BL(~maskfilled));
BR = BR - mean(BR(~maskfilled));

TL = TL - mean(TL(~maskfilled));
TR = TR - mean(TR(~maskfilled));
%% Calculate DF, Edge and QDF
DF_calc = TL + TR + BR +BL;
Edge = abs(double(TL)-double(BR)) + abs(double(BL)-double(TR));
Edge = Edge - mean(Edge(~maskfilled));
ScalingFactor = .8;
QDF = ScalingFactor.*DF_calc - Edge;
%%

% clean bead image
figure(8)

% imagesc(double(BL(300:460,940:1100)))
imagesc(double(BL))

axis off image
colormap gray
caxis([0 5e4/16])
pxlsize=5.316e-4;

scalebar(140,140,10,pxlsize)

figure(9)

imagesc(double(BR(300:460,940:1100)))
axis off image
colormap gray
caxis([0 5e4/16])
pxlsize=5.316e-4;

figure(10)

imagesc(double(TR(300:460,940:1100)))
axis off image
colormap gray
caxis([0 5e4/16])
pxlsize=5.316e-4;


figure(11)

imagesc(double(TL(300:460,940:1100)))
axis off image
colormap gray
caxis([0 5e4/16])
pxlsize=5.316e-4;
%%
figure(12)
imagesc(double(DF_calc(300:460,940:1100)-mean(DF_calc(~maskfilled),'all')))
imagesc(double(DF_calc(300:460,940:1100)))

axis off image
colormap gray
caxis([0 5e4/16])
%%
figure(13)

imagesc(double(Edge(300:460,940:1100)))
axis off image
colormap gray
caxis([0 5e4/16])
%%
figure(14)
imagesc(QDF(300:460,940:1100))
axis off image
colormap gray
caxis([0 2.5e4/16])

%% Not so clean bead
figure(15)

imagesc(double(BL(930:1090,1435:1595)))
axis off image
colormap gray
caxis([0 5e4/16])
pxlsize=5.316e-4;

scalebar(140,140,10,pxlsize)

figure(16)

imagesc(double(BR(930:1090,1435:1595)))
axis off image
colormap gray
caxis([0 5e4/16])
pxlsize=5.316e-4;

figure(17)

imagesc(double(TR(930:1090,1435:1595)))
axis off image
colormap gray
caxis([0 5e4/16])
pxlsize=5.316e-4;


figure(18)

imagesc(double(TL(930:1090,1435:1595)))
axis off image
colormap gray
caxis([0 5e4/16])
pxlsize=5.316e-4;
%%
figure(19)
imagesc(double(DF_calc(930:1090,1435:1595)))
axis off image
colormap gray
caxis([0 5e4/16])
%%
figure(20)
imagesc(double(Edge(930:1090,1435:1595)))
axis off image
colormap gray
caxis([0 5e4/16])
%%
figure(21)
imagesc(QDF(930:1090,1435:1595))
axis off image
colormap gray
caxis([0 2.5e4/16])