classdef Aligner < handle
    properties
        app
    end
    methods
        function obj = Aligner(app)
            obj.app = app;
        end
        function remove_alignment_callback(obj, ~, ~)
            % remove the displacements file
            obj.app.stack_info.aligned = false;
            obj.app.stack_info.displacements = [];
            obj.app.utils.toggle_indicator(obj.app.ui.info.alignedIndicator, false);
            % get current frame
            frame = get(obj.app.ui.controls.slider, 'Value');
            obj.app.utils.setFrame(frame);
            % add star to save button
            obj.app.ui.controls.saveButton.String = 'Save *';
        end
        function skip_alignment_callback(obj, ~, ~)
            % skip the current stack
            obj.app.skip_alignment = true;
        end
        function align_stack_callback(obj, varargin)
            obj.app.skip_alignment = false;
            obj.app.ui.controls.skipButton.Enable = 'on';
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
            if(obj.app.stack_info.img_data.num_imgs == 0)
                obj.app.utils.display_warning("load some images first");
            else
                % Get the template
                frame_num = 1;
                [template,templatePosition] = getTemplate(mode, frame_num);
                % if the template is all zeros, get the next frame
                while all(template(:) == 0)
                    obj.app.utils.display_warning("Template is all zeros, getting next frame");
                    frame_num = frame_num + 1;
                    obj.app.shortener.shorten_stack_callback(0, 0, frame_num, obj.app.stack_info.end_index);
                    [template, templatePosition] = getTemplate(mode, frame_num);
                end
                % match the template with each image in the stack
                displacements = matchTemplate(template, templatePosition);
                if ~obj.app.skip_alignment
                    % subtract the displacements of the first frame from all the displacements
                    displacements(:,1) = displacements(:,1)-displacements(obj.app.stack_info.start_index,1);
                    displacements(:,2) = displacements(:,2)-displacements(obj.app.stack_info.start_index,2);
                    obj.app.stack_info.displacements = displacements;
                    % Display the displacements
                    obj.plot_displacements();
                    % save displacements for later use
                    obj.app.stack_info.aligned = true;
                    obj.app.utils.toggle_indicator(obj.app.ui.info.alignedIndicator, true);
                    obj.app.utils.save_stack_callback();
                end
            end
            % disable skip button
            obj.app.ui.controls.skipButton.Enable = 'off';
            fprintf('Aligned stack %s\n', obj.app.path);
            function displacements = matchTemplate(template, templatePosition)
                displacements = zeros(obj.app.stack_info.img_data.num_imgs, 2);
                % f = waitbar(0,'Please wait...','Name','Aligning stack...');
                WaitMessage = parfor_wait(obj.app.stack_info.img_data.num_imgs, 'Waitbar', true);
                img_data = obj.app.stack_info.img_data;
                parfor k = 1:obj.app.stack_info.img_data.num_imgs
                    if obj.app.skip_alignment
                        error('Alignment skipped');
                    end
                    % obj.app.utils.display_warning(num2str(k));
                    baseFileName = img_data.img_files(k).name;
                    fullFileName = fullfile(img_data.img_files(k).folder, baseFileName);
                    bwImageGPU = gpuArray(imread(fullFileName));
                    templateGPU = gpuArray(template); % Convert the template to a GPU array
                    % Extract the region of interest (ROI) with `n` pixel padding
                    x1 = max(1, round(templatePosition(1)) - obj.app.searchWindow); % X start
                    y1 = max(1, round(templatePosition(2)) - obj.app.searchWindow); % Y start
                    x2 = min(size(bwImageGPU, 2), round(templatePosition(1) + templatePosition(3)) + obj.app.searchWindow);
                    y2 = min(size(bwImageGPU, 1), round(templatePosition(2) + templatePosition(4)) + obj.app.searchWindow);
        
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
            function [template,position] = getTemplate(mode, frame_num)
                % set slider to first image
                image_idx = obj.app.stack_info.start_index + frame_num - 1;
                obj.app.utils.display_warning(sprintf("Getting template from image %d frame %d", image_idx, frame_num));
                obj.app.utils.setFrame(frame_num);
                windowSize = 100;
                x_offset = 5;
                y_offset = 10;
                h = drawrectangle('Parent', obj.app.ui.controls.ax1,'Position',[800-windowSize-x_offset,800-windowSize-y_offset,windowSize,windowSize]);
                if (mode == "manual")
                    % obj.app.utils.display_warning("select a template by drawing a rectangle");
                    obj.app.utils.display_warning("Press enter to confirm the template");
                    wait(h);
                end
                % Wait for the user to press the Enter key to confirm the template
                position = round(h.Position);
                template = imcrop(mat2gray(obj.app.stack_info.img_data.imgs{image_idx}), position);
            end
        end
        function plot_displacements(obj)
            displacements = obj.app.stack_info.displacements;
            % clear axis
            cla(obj.app.ui.controls.ax2);
            % plot the displacements on obj.app.ui.controls.ax2
            plot(obj.app.ui.controls.ax2, displacements(:,1), 'r');
            hold(obj.app.ui.controls.ax2, 'on');
            plot(obj.app.ui.controls.ax2, displacements(:,2), 'b');
            % draw vertical lines at start and end indices
            if obj.app.stack_info.shortened == true
                % disp(obj.app.stack_info.start_index);
                plot(obj.app.ui.controls.ax2, obj.app.stack_info.start_index, [-5, 5], 'g');
                plot(obj.app.ui.controls.ax2, obj.app.stack_info.end_index, [-5, 5], 'g');
            end
            title('Displacements');
            xlabel('Image number');
            ylabel('Displacement');
            legend('x', 'y');
            axis(obj.app.ui.controls.ax2, 'tight'); 
        end
        function align_all_stacks_callback(obj, ~,~)
            % iterate over all the stacks and align them
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.load_images_callback();
                if obj.app.stack_info.aligned == true
                    obj.app.logs{end+1} = sprintf("Trial %s already aligned",obj.app.stack_paths(i));
                    continue;
                end
                if contains(obj.app.path, 'time_control')
                    obj.app.stack_info.aligned = true;
                    obj.app.utils.save_stack_callback();
                    continue;
                end
                obj.align_stack_callback('mode', 'auto');          
            end
        end
    end
end