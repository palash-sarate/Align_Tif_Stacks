classdef Image_viewer < handle
    properties
        stack_info
        utils
    end
    methods
        % constructor
        function obj = Image_viewer(stack_info, utils)
            obj.stack_info = stack_info;
            obj.utils = utils;    
        end
        function setFrame(k)
            persistent counter;
            % initialize the counter
            if isempty(counter)
                counter = 0;
            end
            % reset imgs after 100 frames
            if counter > 100
                % display_warning("Resetting images");
                counter = 0;
                stack_info.img_data.imgs = cell(1, stack_info.img_data.num_imgs);
            end
            image_idx = stack_info.start_index + k - 1;
            % fprintf('Frame number: %d\n image idx %d', k, image_idx);
            % Get the current slider value
            set(ui.frame_number, 'String', num2str(k)); % Update the frame number display
            set(ui.slider, 'Value', k);  % update the slider value
            % load the image if it's not already loaded
            if isempty(stack_info.img_data.imgs{image_idx})            
                % fprintf('Loading image %d as its empty', image_idx);
                stack_info.img_data.imgs{image_idx} = mat2gray(imread(fullfile(stack_info.img_data.img_files(image_idx).folder, stack_info.img_data.img_files(image_idx).name)));
                counter = counter + 1;
            end
            % if aligned displacements exist, apply them to the image
            if stack_info.aligned == true
                displaced_img = imtranslate(stack_info.img_data.imgs{image_idx}, - stack_info.displacements(image_idx, :));
                if isfield(stack_info, 'masked') && stack_info.masked == true && isfield(stack_info, 'mask')
                    displaced_img = imoverlay(displaced_img, stack_info.mask, 'r');
                end 
                imshow(displaced_img, 'Parent', ui.ax1);
            else
                imshow(stack_info.img_data.imgs{image_idx}, 'Parent', ui.ax1);
            end
            if isfield(stack_info, 'timestamps')
                timestamp = stack_info.timestamps{image_idx};
                if ~isempty(timestamp)
                    % draw the timestamp on the image
                    % disp(image_idx);
                    text(ui.ax1, 'Units', 'normalized', 'Position', [0.99, 0.03], ...
                        'String', sprintf("%.2f Sec", ...
                        time_2_sec(timestamp)-time_2_sec(stack_info.timestamps{1})), ...
                        'Color', 'white', 'FontSize', 18, 'HorizontalAlignment', 'right');
                end
            end
        end
    end 
end