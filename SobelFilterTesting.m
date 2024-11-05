clear, close all

% Copyright Tarek Moustafa (Zangle lab, University of Utah) 2024

%% Set to 0 to display the full image, or 1 to display only the ROIs
use_ROI = 1;

%% Load files
BL = double(imread('Data/Beads/QrtDF_BL.tiff'));
BR = double(imread('Data/Beads/QrtDF_BR.tiff'));
TL = double(imread('Data/Beads/QrtDF_TL.tiff'));
TR = double(imread('Data/Beads/QrtDF_TR.tiff'));

%% Scale to 12 bit
BL = ceil(BL ./ 16);
BR = ceil(BR ./ 16);
TL = ceil(TL ./ 16);
TR = ceil(TR ./ 16);

%% Pixel size
pxlsize = 5.316e-4; % mm

%% Make a mask for background
DF_Mask = TL + BR + BL + TR;
mask = DF_Mask >= 2e4 / 16;
maskfilled = imfill(mask, 'holes');

%% Correct Background for all
BL = BL - mean(BL(~maskfilled));
BR = BR - mean(BR(~maskfilled));
TL = TL - mean(TL(~maskfilled));
TR = TR - mean(TR(~maskfilled));

%% Test
%% Calculate DF, Edge and QDF
DF_calc = TL + TR + BR + BL;
Edge = abs(double(TL) - double(BR)) + abs(double(BL) - double(TR));
Edge = Edge - mean(Edge(~maskfilled));
ScalingFactor = 0.8;
QDF = ScalingFactor .* DF_calc - Edge;

%% Define Regions of Interest (ROIs)
ROIs = {
    DF_calc(300:460, 940:1100), ...
    DF_calc(930:1090, 1435:1595)
};
Edge_ROIs = {
    Edge(300:460, 940:1100), ...
    Edge(930:1090, 1435:1595)
};
QDF_ROIs = {
    QDF(300:460, 940:1100), ...
    QDF(930:1090, 1435:1595)
};

% Names for the two ROIs
roi_names = {"Clean Bead", "Bead with Imperfections"};

% Best Sobel thresholds
sobel_thresholds = [200, 430, 600];

% Structuring element sizes for dilation
dilation_sizes = [3,5];  % Test sizes from 2 to 5

%% Sobel Edge Detection Figures
if use_ROI == 1
    % Use ROIs for processing
    data_sets = {ROIs, Edge_ROIs, QDF_ROIs};
    data_names = roi_names;
else
    % Use the full images for processing
    data_sets = {{DF_calc}, {Edge}, {QDF}};
    data_names = {"Full Image"};
end

% Loop over each dataset (either full image or ROIs)
for data_idx = 1:length(data_sets{1})
    DF = data_sets{1}{data_idx};
    Edge_ROI = data_sets{2}{data_idx};
    QDF_ROI = data_sets{3}{data_idx};
    
    % Apply Gaussian filter to smooth the image for Sobel detection
    I_smooth = imgaussfilt(DF, 2);

    % Loop over each Sobel threshold
    for i = 1:length(sobel_thresholds)
        % Loop over each structuring element size for dilation
        for se_size = dilation_sizes
            % Create structuring element
            se = strel('disk', se_size);
            
            % Sobel edge detection
            BW_sobel = edge(I_smooth, 'sobel', sobel_thresholds(i));
            BW_sobel_dilated = imdilate(BW_sobel, se);  % Apply dilation
            masked_image_sobel = DF;  % Original image for Sobel
            masked_image_sobel(~BW_sobel_dilated) = 0;  % Masked image after dilation
            sobel_difference = DF - masked_image_sobel;  % Difference image

            % Create a new figure for the current threshold and dilation size, and maximize it
            fig = figure;
            set(fig, 'WindowState', 'maximized');

            % Set a common color scale limit
            max_val = prctile([DF(:); Edge_ROI(:); QDF_ROI(:); masked_image_sobel(:); sobel_difference(:)], 99); % 99th percentile
            min_val = 0;  % Set minimum value to 0
            
            % First Row: DF, Edge, QDF
            % Column 1: DF
            subplot(2, 3, 1);
            imagesc(DF, [min_val max_val]);  % Use same color scale
            colormap(gray);
            axis image off;
            title('DF');

            % Column 2: Edge
            subplot(2, 3, 2);
            imagesc(Edge_ROI, [min_val max_val]);  % Use same color scale
            colormap(gray);
            axis image off;
            title('Edge');

            % Column 3: QDF
            subplot(2, 3, 3);
            imagesc(QDF_ROI, [min_val max_val]);  % Use same color scale
            colormap(gray);
            axis image off;
            title('QDF');

            % Second Row: Original Image, Masked after Dilation, Difference
            % Column 1: Original Image (DF)
            subplot(2, 3, 4);
            imagesc(DF, [min_val max_val]);  % Use same color scale
            colormap(gray);
            axis image off;
            title('DF');

            % Column 2: Masked Image after Dilation (Edge Detection)
            subplot(2, 3, 5);
            imagesc(masked_image_sobel, [min_val max_val]);  % Use same color scale
            colormap(gray);
            axis image off;
            title('Edge (Digital)');

            % Column 3: Difference Image
            subplot(2, 3, 6);
            imagesc(sobel_difference, [min_val max_val]);  % Use same color scale
            colormap(gray);
            axis image off;
            title('DF without Edge (Digital)');
            
            % Set a common title for the figure
            sgtitle([data_names{data_idx}, ' - Sobel Threshold: ', num2str(sobel_thresholds(i)), ...
                     ', Dilation Size: ', num2str(se_size)]);
        end
    end
end
