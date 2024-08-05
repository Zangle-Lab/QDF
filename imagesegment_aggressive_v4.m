function [L] = imagesegment_aggressive_v4(I)
%function [L] = imagesegment_aggressive_v3(I)
%function to segment an image of cells using a watershed image
%transform
%input: I, the grayscale image to segment;
%output: L, the labeled, segmented image
%modified heavily for DPC image analysis

% I = medfilt2(Phase*625-Btotal);
%step 2, detect cell

If1 = imfilter(I, fspecial('sobel'));
If2 = imfilter(I, fspecial('sobel')');
% BWs = adaptivethreshold(abs(If1),[400 400], 0) | adaptivethreshold(abs(If2),[400 400], 0);
BWs = abs(If1)>120 | abs(If2)>120;

%step 3, dilate the image
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);

BWsdil = imdilate(BWs, [se90 se0]);
% figure(3)
% imshow(BWsdil)
% title('dilated gradient mask')

%step 4, fill gaps
BWdfill = imfill(BWsdil, 'holes');
% figure(4)
% imshow(BWdfill);
% title('binary image with filled holes');



%step 6, smooth image
seD = strel('diamond',1);
%BWfinal = imerode(BWnobord,seD);
BWdfErode = imerode(BWdfill,seD);

%step 6: dilate image

se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);

BWsdil2 = imdilate(BWdfErode, [se90 se0]);

BWclear = bwareaopen(BWsdil2, 300); %remove regions smaller than 300 pixels


se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);

BWsdil3 = imdilate(BWclear, [se90 se0]);
BWclear = imclearborder(BWsdil3);

% figure(6)
% imshow(BWfinal)
% title('segmented image')

%to segment connected cells, follow this page:
%http://blogs.mathworks.com/steve/2006/06/02/cell-segmentation/
%new make a mask for the watershed transform
% % Igr = mat2gray(I);
% thr = graythresh(Igr);
% % thr = 0.1; % originally .1

%%%%%%%%% adjustment in next line
% % mask_em = imextendedmax(Igr, thr); %change this to 0.8 or so
%%%%%%%%%
% imagesc(mask_em)

% mask_em = imclose(mask_em, ones(4,4)); %this was ones(5,5) in original imagesegment
% mask_em = imerode(mask_em, strel('diamond', 1));
% % mask_em = imfill(mask_em, 'holes');
% mask_em = bwareaopen(mask_em, 1);%this parameter was originally 40, depends on mag used
%Next step: complement the image so that the peaks become valleys. We do
%this because we are about to apply the watershed transform, which
%identifies low points, not high points. 
% % I_c = imcomplement(I);
%Next: modify the image so that the background pixels and the extended
%maxima pixels are forced to be the only local minima in the image. 
% % I_mod = imimposemin(I_c, ~BWfinal | mask_em);
%now, compute watershed transform
% % L = watershed(I_mod);
% imshow(label2rgb(L))

% L1 = imclearborder(BWfinal);

L = bwlabel(BWclear);
