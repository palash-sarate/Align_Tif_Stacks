classdef Shortener < handle
    properties
        app
    end
    methods
        function obj = Shortener(app)
            obj.app = app;
        end
        function [start_index, end_index] = shorten_stack_callback(obj, ~, ~, start_, end_)
            obj.remove_shortening_callback();
            % check if stack_info has maxDiffIndex
            if ~isfield(obj.app.stack_info, 'maxDiffIndex')
                maxDiffIndex = matchImages(300);
                obj.app.utils.display_warning(sprintf("Max difference index: %d", maxDiffIndex));
            end
            obj.app.utils.setFrame(obj.app.stack_info.maxDiffIndex);
            if nargin > 2
                start_index = start_;
                end_index = end_;
            else
                obj.app.utils.display_warning("go to start frame and press n");
                obj.app.utils.wait_for_keypress("n");
                start_index = obj.app.stack_info.start_index + round(get(obj.app.ui.controls.slider, 'Value')) - 1;
                obj.app.utils.display_warning("go to end frame and press n");
                obj.app.utils.setFrame(obj.app.stack_info.end_index-obj.app.stack_info.start_index + 1);
                obj.app.utils.wait_for_keypress("n");
                end_index = obj.app.stack_info.start_index + round(get(obj.app.ui.controls.slider, 'Value')) - 1;
            end
            obj.app.masking.shorten_slider(start_index, end_index);
            obj.app.stack_info.shortened = true;
            obj.app.utils.toggle_indicator(obj.app.ui.info.shortenedIndicator, true);
            % save the start and end indices in a obj.app.stack_info.mat file
            obj.app.utils.save_stack_callback();
            % GPU accelerated function to compare first n images with previous image
            function maxDiffIndex = matchImages(n)
                differences = zeros(n, 1);
                % f = waitbar(0,'Please wait...','Name','Aligning stack...');
                WaitMessage = parfor_wait(n, 'Waitbar', true);
                img_data = obj.app.stack_info.img_data;
                roi = [200, 750, 400, 100]; % [x, y, width, height]
                for k = 2:n
                    % obj.app.utils.display_warning(num2str(k));
                    baseFileName = img_data.img_files(k-1).name;
                    fullFileName = fullfile(img_data.img_files(k-1).folder, baseFileName);
                    template = imread(fullFileName);
                    template = imcrop(mat2gray(template), roi);
                    displaced_img1 = imtranslate(template, - obj.app.stack_info.displacements(k-1, :));
                    % imshow(displaced_img1, 'Parent', obj.app.ui.controls.ax1);
                    baseFileName = img_data.img_files(k).name;
                    fullFileName = fullfile(img_data.img_files(k).folder, baseFileName);
                    bwImage = imread(fullFileName);
                    bwImage = imcrop(mat2gray(bwImage), roi);
                    displaced_img2 = imtranslate(bwImage, - obj.app.stack_info.displacements(k, :));

                    difference = imabsdiff(displaced_img1, displaced_img2);
                    differences(k) = sum(difference(:));
                    WaitMessage.Send;
                end
                WaitMessage.Destroy
                maxDiffIndex = find(differences == max(differences));
                obj.app.stack_info.maxDiffIndex = maxDiffIndex;
                obj.app.utils.save_stack_callback();
            end
        end
        function shorten_all_stack_callback(obj, ~,~)
            % iterate over all the stacks and align them
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.load_images_callback();
                if isfield(obj.app.stack_info, 'shortened')
                    if obj.app.stack_info.shortened == true
                        obj.app.logs{end+1} = sprintf("Trial %s already shortened",obj.app.stack_paths(i));
                        continue;
                    end
                end
                if contains(obj.app.path, 'time_control')
                    obj.app.stack_info.shortened = true;
                    obj.app.utils.save_stack_callback();
                    continue;
                end
                obj.shorten_stack_callback();
            end
        end
        function remove_shortening_callback(obj, ~, ~)
            % remove the displacements file
            obj.app.stack_info.shortened = false;
            obj.app.stack_info.start_index = 1;
            obj.app.stack_info.end_index = numel(obj.app.stack_info.img_data.img_files);
            obj.app.masking.shorten_slider(1, numel(obj.app.stack_info.img_data.img_files));
            obj.app.utils.toggle_indicator(obj.app.ui.info.shortenedIndicator, false);
            % get current frame
            frame = get(obj.app.ui.controls.slider, 'Value');
            obj.app.utils.setFrame(frame);
            % add star to save button
            obj.app.ui.controls.saveButton.String = 'Save *';
        end
    end
end