function [B,M,D2_DF,D2_CorrDF] = imagebackground_polyn_reduced2_DF_v1(I, TL,TR,BR,BL,order, r1, r2)
%function to find the background of an image, I, using a "orderth" order
%polynomial fit to the background pixels
%_reduced2 also reduces the number of points at which to evaluate polyfitn
%by a factor of 'r1', as well as the polyvaln at the end by factor of r2
%reasonable values from initial testing on DPC: r1 = 50, r2 =10
%_v2 modified TAZ 10/16/20 to include cells near edges and tighten borders
%around cells
%input: I, the grayscale image to find the background of
%output: B, the background of I; M, the mask used to find B
%method: find 'objects' in I, mask them from the image, paint the remaining
%area using the inpaint_nans function

%%setup variables for later
sz = size(I);
[X,Y] = meshgrid(1:sz(2), 1:sz(1));
%reduce the list of points for polyvaln
Xfit = imresize(X, 1./r2);
Yfit = imresize(Y, 1./r2);
szfit = size(Xfit);

%%fit surface to entire image using r1 and r2 parameters, then make
%initial 'profile' image from iterative polyfit
RandMask = logical(0.*I+1);
RandMask(randi([0 50-1], size(I))~=0) = 0; %reduce # points by factor of 50

XList = X(RandMask(:));
YList = Y(RandMask(:));    
IList = I(RandMask(:));
CFit = polyfitn([XList,YList], IList, order);

PolyProfile =(reshape(polyvaln(CFit, [Xfit(:), Yfit(:)]), szfit(1), szfit(2)));
PolyProfile = imresize(PolyProfile,size(I));

D1 = I - min(I,PolyProfile);

TL_med = medfilt2(TL,[5 5]);
TR_med = medfilt2(TR,[5 5]);
BL_med = medfilt2(BL,[5 5]);
BR_med = medfilt2(BR,[5 5]);


C_TR = TR_med>100;
C_TL = TL_med>100;
C_BL = BL_med>100;
C_BR = BR_med>100;
%%run edge detection on image after this initial correction
BWs = edge(D1, 'prewitt', 20);

%dilate edge detection image slightly




%subplot(1,2,1)
%imagesc(BWs)
%axis image

% create structuring element for erosion, to seperate large regions
% from the edges
ErStrel0 = strel('line', 2, 0);
ErStrel90 = strel('line', 2, 90);

% erode mask using structuring element
BWdil1 = imdilate(BWs, [ErStrel0 ErStrel90]);
% ErStrel0 = strel('line', 2, 0);
% ErStrel90 = strel('line', 2, 90);
% BWdil2 = imdilate(BWdil1, [ErStrel0 ErStrel90]);
BWdil2 = BWdil1;

% BWclear = imclearborder(BWdil2);
%manually clear border pixels instead of removing entire regions (cells)
BWclear = BWdil2;
BWclear(1:7,:) = 0;
BWclear(end-6:end,:) = 0;
BWclear(:,1:7) = 0;
BWclear(:,end-6:end) = 0;
BWopen = bwareaopen(BWclear, 500);

% dilate cells
seCirc = strel('disk', 6);
BWdil3 = imdilate(BWopen, seCirc);
% BWdil4 = imdilate(BWdil3, seCirc);
% ErStrel0 = strel('line', 4, 0);
% ErStrel90 = strel('line', 4, 90);
% BWdil5 = imdilate(BWdil3, [ErStrel0 ErStrel90]);
BWfilled = imfill(BWdil3, 'holes');

BWfinal = BWfilled;

%reduce BWfinal with random mask (reduce by factor of reduce)
% reduce = 4; 
BWfinal(randi([0 r1-1], size(I))~=0) = 1;

IList = I(~BWfinal);

TLList= TL_med(~BWfinal);
TRList= TR_med(~BWfinal);
BLList= BL_med(~BWfinal);
BRList= BR_med(~BWfinal);


sz = size(I);
[X,Y] = meshgrid(1:sz(2), 1:sz(1));
XList = X(~BWfinal);
YList = Y(~BWfinal);

%reduce the list of points for polyvaln
Xfit = imresize(X, 1./r2);
Yfit = imresize(Y, 1./r2);
szfit = size(Xfit);

if sum(~isnan(BWfinal)) ~= 0
    CFit = polyfitn([XList,YList], IList, order);
    CFit_TL = polyfitn([XList,YList], TLList, order);
    CFit_TR = polyfitn([XList,YList], TRList, order);
    CFit_BL = polyfitn([XList,YList], BLList, order);
    CFit_BR = polyfitn([XList,YList], BRList, order);

%     B2 =(reshape(polyvaln(CFit, [X(:), Y(:)]), sz(1), sz(2)));
    B =(reshape(polyvaln(CFit, [Xfit(:), Yfit(:)]), szfit(1), szfit(2)));
    B_TL = (reshape(polyvaln(CFit_TL, [Xfit(:), Yfit(:)]), szfit(1), szfit(2)));
    B_TR = (reshape(polyvaln(CFit_TR, [Xfit(:), Yfit(:)]), szfit(1), szfit(2)));
    B_BL = (reshape(polyvaln(CFit_BL, [Xfit(:), Yfit(:)]), szfit(1), szfit(2)));
    B_BR = (reshape(polyvaln(CFit_BR, [Xfit(:), Yfit(:)]), szfit(1), szfit(2)));
    B = imresize(B,size(I));
    B_TL = imresize(B_TL,size(I));
    B_TR = imresize(B_TR,size(I));
    B_BL = imresize(B_BL,size(I));
    B_BR = imresize(B_BR,size(I));

else
    B = I;
    B_TL =TL_med;
    B_TR =TR_med;
    B_BL =BL_med;
    B_BR =BR_med;
end

M = ~BWfilled; %return the mask used for processing


D1_TL = TL-B_TL;
D1_TR = TR-B_TR;
D1_BL = BL-B_BL;
D1_BR = BR-B_BR;

D2_TL = D1_TL - mean(D1_TL(M));
D2_TR = D1_TR - mean(D1_TR(M));
D2_BL = D1_BL - mean(D1_BL(M));
D2_BR = D1_BR - mean(D1_BR(M));

DF = D2_TL + D2_TR + D2_BL + D2_BR;
CorrDF = abs(D2_TL - D2_BR) + abs(D2_TR - D2_BL);

CorrDF_med= medfilt2(CorrDF,[5 5]);
DF_med = medfilt2(DF,[5 5]);



CorrDFList= CorrDF_med(~BWfinal);
DFList    = DF(~BWfinal);

if sum(~isnan(BWfinal)) ~= 0
    CFit_CorrDF = polyfitn([XList,YList], CorrDFList, order);
    CFit_DF     = polyfitn([XList,YList], DFList, order);

%     B =(reshape(polyvaln(CFit, [X(:), Y(:)]), sz(1), sz(2)));
    B_CorrDF =(reshape(polyvaln(CFit_CorrDF, [Xfit(:), Yfit(:)]), szfit(1), szfit(2)));
    B_DF =(reshape(polyvaln(CFit_DF, [Xfit(:), Yfit(:)]), szfit(1), szfit(2)));

    B_CorrDF = imresize(B_CorrDF,size(I));
    B_DF = imresize(B_DF,size(I));

else
    B_CorrDF = CorrDF_med;
    B_DF = DF_med;
end
D1_CorrDF = CorrDF-B_CorrDF;
D1_DF     = DF- B_DF;


D2_CorrDF = D1_CorrDF - mean(D1_CorrDF(M));
D2_DF     = D1_DF     - mean(D1_DF(M));