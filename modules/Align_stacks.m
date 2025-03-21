classdef Align_stacks < handle
    properties
        stack_info
        utils
        Gr_module
    end
    methods
        % constructor
        function obj = Align_stacks(stack_info, utils)
            obj.stack_info = stack_info;
            obj.utils = utils;
            obj.Gr_module = Gr_module(stack_info);
        end
        function displacements = matchTemplate(template, templatePosition)
            displacements = zeros(stack_info.img_data.num_imgs, 2);
            % f = waitbar(0,'Please wait...','Name','Aligning stack...');
            WaitMessage = parfor_wait(stack_info.img_data.num_imgs, 'Waitbar', true);
            img_data = stack_info.img_data;
            parfor k = 1:stack_info.img_data.num_imgs
                if skip_alignment
                    error('Alignment skipped');
                end
                % display_warning(num2str(k));
                baseFileName = img_data.img_files(k).name;
                fullFileName = fullfile(img_data.img_files(k).folder, baseFileName);
                bwImageGPU = gpuArray(imread(fullFileName));
                templateGPU = gpuArray(template); % Convert the template to a GPU array
                % Extract the region of interest (ROI) with `n` pixel padding
                x1 = max(1, round(templatePosition(1)) - searchWindow); % X start
                y1 = max(1, round(templatePosition(2)) - searchWindow); % Y start
                x2 = min(size(bwImageGPU, 2), round(templatePosition(1) + templatePosition(3)) + searchWindow);
                y2 = min(size(bwImageGPU, 1), round(templatePosition(2) + templatePosition(4)) + searchWindow);
    
                % Crop the padded region
                croppedImageGPU = bwImageGPU(y1:y2, x1:x2);
                % find the template in the image
                c = normxcorr2(templateGPU, croppedImageGPU);
                [ypeak, xpeak] = find(gather(c) == max(c(:)), 1);
                yoffSet = ypeak-size(template,1);
                xoffSet = xpeak-size(template,2);
                displacements(k, :) = [xoffSet, yoffSet];
                WaitMessage.Send;
            end
            WaitMessage.Destroy
        end
        % GPU accelerated function to compare first n images with previous image
        function maxDiffIndex = matchImages(n)
            differences = zeros(n, 1);
            % f = waitbar(0,'Please wait...','Name','Aligning stack...');
            WaitMessage = parfor_wait(n, 'Waitbar', true);
            img_data = stack_info.img_data;
            roi = [200, 750, 400, 100]; % [x, y, width, height]
            for k = 2:n
                % display_warning(num2str(k));
                baseFileName = img_data.img_files(k-1).name;
                fullFileName = fullfile(img_data.img_files(k-1).folder, baseFileName);
                template = imread(fullFileName);
                template = imcrop(mat2gray(template), roi);
                displaced_img1 = imtranslate(template, - stack_info.displacements(k-1, :));
                % imshow(displaced_img1, 'Parent', ax1);
                baseFileName = img_data.img_files(k).name;
                fullFileName = fullfile(img_data.img_files(k).folder, baseFileName);
                bwImage = imread(fullFileName);
                bwImage = imcrop(mat2gray(bwImage), roi);
                displaced_img2 = imtranslate(bwImage, - stack_info.displacements(k, :));
    
                difference = imabsdiff(displaced_img1, displaced_img2);
                differences(k) = sum(difference(:));
                WaitMessage.Send;
            end
            WaitMessage.Destroy
            maxDiffIndex = find(differences == max(differences));
            stack_info.maxDiffIndex = maxDiffIndex;
            save_stack_callback();
        end
        function [template,position] = getTemplate(mode, frame_num)
            % set slider to first image
            image_idx = stack_info.start_index + frame_num - 1;
            display_warning(sprintf("Getting template from image %d frame %d", image_idx, frame_num));
            setFrame(frame_num);
            windowSize = 100;
            x_offset = 5;
            y_offset = 10;
            h = drawrectangle('Parent', ui.ax1,'Position',[800-windowSize-x_offset,800-windowSize-y_offset,windowSize,windowSize]);
            if (mode == "manual")
                % display_warning("select a template by drawing a rectangle");
                display_warning("Press enter to confirm the template");
                wait(h);
            end
            % Wait for the user to press the Enter key to confirm the template
            position = round(h.Position);
            template = imcrop(mat2gray(stack_info.img_data.imgs{image_idx}), position);
        end
        function skip_alignment_callback(~, ~)
            % skip the current stack
            skip_alignment = true;
        end
        % load_images_callback();
        function plot_displacements()
            displacements = stack_info.displacements;
            % clear axis
            cla(ui.ax2);
            % plot the displacements on ui.ax2
            plot(ui.ax2, displacements(:,1), 'r');
            hold(ui.ax2, 'on');
            plot(ui.ax2, displacements(:,2), 'b');
            % draw vertical lines at start and end indices
            if stack_info.shortened == true
                % disp(stack_info.start_index);
                plot(ui.ax2, stack_info.start_index, [-5, 5], 'g');
                plot(ui.ax2, stack_info.end_index, [-5, 5], 'g');
            end
            title('Displacements');
            xlabel('Image number');
            ylabel('Displacement');
            legend('x', 'y');
            axis(ui.ax2, 'tight'); 
        end
        function align_stack_callback(varargin)
            skip_alignment = false;
            skip_button.Enable = 'on';
            % Check if the function is called as a callback
            if nargin > 0 && isa(varargin{1},'matlab.ui.control.UIControl')
                % Called as a callback, set default mode
                mode = 'manual';
            else
                % Called with arguments, use inputParser to parse the arguments
                p = inputParser;
                addParameter(p, 'mode', 'manual', @(x) ischar(x) || isstring(x));
                parse(p, varargin{:});
                mode = char(p.Results.mode); % Ensure mode is a character array
            end
            if(stack_info.img_data.num_imgs == 0)
                display_warning("load some images first");
            else
                % Get the template
                frame_num = 1;
                [template,templatePosition] = getTemplate(mode, frame_num);
                % if the template is all zeros, get the next frame
                while all(template(:) == 0)
                    display_warning("Template is all zeros, getting next frame");
                    frame_num = frame_num + 1;
                    shorten_stack_callback(0, 0, frame_num, stack_info.end_index);
                    [template, templatePosition] = getTemplate(mode, frame_num);
                end
                % match the template with each image in the stack
                displacements = matchTemplate(template, templatePosition);
                if ~skip_alignment
                    % subtract the displacements of the first frame from all the displacements
                    displacements(:,1) = displacements(:,1)-displacements(stack_info.start_index,1);
                    displacements(:,2) = displacements(:,2)-displacements(stack_info.start_index,2);
                    stack_info.displacements = displacements;
                    % Display the displacements
                    plot_displacements();
                    % save displacements for later use
                    stack_info.aligned = true;
                    toggle_indicator(ui.aligned_indicator, true);
                    save_stack_callback();
                end
            end
            % disable skip button
            skip_button.Enable = 'off';
            fprintf('Aligned stack %s\n', path);
        end
        function remove_alignment_callback(~, ~)
            % remove the displacements file
            stack_info.aligned = false;
            stack_info.displacements = [];
            toggle_indicator(ui.aligned_indicator, false);
            % get current frame
            frame = get(ui.slider, 'Value');
            setFrame(frame);
            % add star to save button
            save_button.String = 'Save *';
        end
        function align_all_stacks_callback(~,~)
            % iterate over all the stacks and align them
            for i = 1:length(stack_paths)
                set(ui.stack_dropdown, 'Value', i);
                load_images_callback();
                if stack_info.aligned == true
                    logs{end+1} = sprintf("Trial %s already aligned",stack_paths(i));
                    continue;
                end
                if contains(path, 'time_control')
                    stack_info.aligned = true;
                    save_stack_callback();
                    continue;
                end
                align_stack_callback('mode', 'auto');          
            end
        end
    end
end