clear all, close all

load('data61.mat','DF_stored')
DF= double(DF_stored(200:260,610:670,12));

%200:260,610:670
I_smooth = imgaussfilt(DF, 2);

% Best Sobel thresholds
sobel_thresholds = 1;

% Structuring element sizes for dilation
dilation_sizes = [5];  % Test sizes from 2 to 5

output_dir = 'Results_cell_MDAMB231';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

    for i = 1:length(sobel_thresholds)
        for se_size = dilation_sizes
            se = strel('disk', se_size);
            figure
            % Sobel edge detection
            BW_sobel = edge(I_smooth, 'sobel', sobel_thresholds(i));
            BW_sobel_dilated = imdilate(BW_sobel, se);
            masked_image_sobel = DF;
            masked_image_sobel(~BW_sobel_dilated) = 0;
            sobel_difference = DF - masked_image_sobel;

            % Set color scale limits
            min_val = 0;
            max_val = 400;
            
                        % First Row: DF, Edge, QDF
            % Column 1: DF
            subplot(2, 3, 1);
            imagesc(DF, [min_val max_val]);  % Use same color scale
            colormap(gray);
            axis image off;
            title('DF');
% 
% %             % Column 2: Edge
% %             subplot(2, 3, 2);
% %             imagesc(Edge_ROI, [min_val max_val]);  % Use same color scale
% %             colormap(gray);
% %             axis image off;
% %             title('Edge');
% % 
% %             % Column 3: QDF
% %             subplot(2, 3, 3);
% %             imagesc(QDF_ROI, [min_val max_val]);  % Use same color scale
% %             colormap(gray);
% %             axis image off;
% %             title('QDF');

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
            sgtitle([' - Sobel Threshold: ', num2str(sobel_thresholds(i)), ...
                     ', Dilation Size: ', num2str(se_size)]);
            % Construct the full filename with the directory path
            base_filename = sprintf('Threshold_%d_Dilation_%d', sobel_thresholds(i), se_size);

            % Save each file in the specified output directory
%             imwrite(uint16(mat2gray(DF, [min_val max_val]) * 65535), fullfile(output_dir, sprintf('%s_DF.tiff', base_filename)));
%             imwrite(uint16(mat2gray(Edge_ROI, [min_val max_val]) * 65535), fullfile(output_dir, sprintf('%s_Edge.tiff', base_filename)));
%             imwrite(uint16(mat2gray(QDF_ROI, [min_val max_val]) * 65535), fullfile(output_dir, sprintf('%s_QDF.tiff', base_filename)));
%             imwrite(uint16(mat2gray(masked_image_sobel, [min_val max_val]) * 65535), fullfile(output_dir, sprintf('%s_Edge_Masked.tiff', base_filename)));
%             imwrite(uint16(mat2gray(sobel_difference, [min_val max_val]) * 65535), fullfile(output_dir, sprintf('%s_Difference.tiff', base_filename)));
        end
    end
% end
