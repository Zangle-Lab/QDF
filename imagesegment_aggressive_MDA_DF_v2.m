function [L,L_tight] = imagesegment_aggressive_MDA_DF_v2(I)
%function [L] = imagesegment_aggressive_smaller2(I)
%function to segment an image of cells using the matlab help file
%procedure followed by breaking up connected cells using a watershed image
%transform
%input: I, the grayscale image to segment
%output: L, the labeled, segmented image
%smaller2 uses a gaussian filter before watershed

%step 2, detect entire cell
% [junk threshold] = edge(I, 'sobel');
% fudgeFactor = .5;
% BWs = edge(I,'sobel', threshold * fudgeFactor);

% BWs = adaptivethreshold(abs(I),[400 400], 0);

If1 = imfilter(I, fspecial('sobel'));
If2 = imfilter(I, fspecial('sobel')');
t= 60;
BWs = abs(If1)>t | abs(If2)>t;
%fix problem with borders of image
BWs(1,:) = 0;
BWs(end,:) = 0;
BWs(:,1) = 0;
BWs(:,end) = 0;


% figure(2)
% imshow(BWs)
% title('binary gradient mask');

%step 2, Erode to remove wave lines in the image
% seD = strel('diamond',3);
% BWsrod = imerode(BWs, seD);
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWsrod= imerode(BWs,[se90 se0]);

% se90 = strel('line', 2, 90);
% se0 = strel('line', 2, 0);
% BWsrod2= imerode(BWsrod,[se90 se0]);
%step 3, dilate the image
se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);

BWsdil = imdilate(BWsrod, [se90 se0]);

se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);
BWsdil2 = imdilate(BWsdil, [se90 se0]);
% figure(3)
% imshow(BWsdil)
% title('dilated gradient mask')

%step 4, fill gaps
BWdfill = imfill(BWsdil2, 'holes');
% figure(4)
% imshow(BWdfill);
% title('binary image with filled holes');



%step 6, smooth image
seD = strel('diamond',3);
%BWfinal = imerode(BWnobord,seD);
BWfiller = imerode(BWdfill,seD);
BWcleaner = bwareaopen(BWfiller, 300); %remove regions smaller than 300 pixels
seD = strel('diamond',3);
BWfinal = imdilate(BWcleaner, seD);

%to segment connected cells, follow this page:
%http://blogs.mathworks.com/steve/2006/06/02/cell-segmentation/
%new make a mask for the watershed transform
Igr = mat2gray(imfilter(I, fspecial('gaussian', [20 20], 2)));
thr = .05;

%%%%%%%%% adjustment in next line
mask_em_1 = imextendedmax(Igr, thr); 
%%%%%%%%%
mask_em = imfill(mask_em_1, 'holes');
%Next step: complement the image so that the peaks become valleys. We do
%this because we are about to apply the watershed transform, which
%identifies low points, not high points.
I_c = imcomplement(I);
%Next: modify the image so that the background pixels and the extended
%maxima pixels are forced to be the only local minima in the image.
I_mod = imimposemin(I_c, ~BWfinal | mask_em);
%now, compute watershed transform
L = watershed(I_mod);

%step 5, remove connected images on border
L = imclearborder(L, 4);

seD = strel('diamond',2);

L_tight = imerode(L,seD);
