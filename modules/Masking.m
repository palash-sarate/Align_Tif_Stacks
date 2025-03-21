classdef Masking < handle
    properties
        stack_info
        utils
    end
    methods
        % constructor
        function obj = Masking(stack_info, utils)
            obj.stack_info = stack_info;
            obj.utils = utils;    
        end
        function [start_index, end_index] = shorten_stack_callback(~, ~, start_, end_)
            remove_shortening_callback();
            % check if stack_info has maxDiffIndex
            if ~isfield(stack_info, 'maxDiffIndex')
                maxDiffIndex = matchImages(300);
                display_warning(sprintf("Max difference index: %d", maxDiffIndex));
            end
            setFrame(stack_info.maxDiffIndex);
            if nargin > 2
                start_index = start_;
                end_index = end_;
            else
                display_warning("go to start frame and press n");
                wait_for_keypress("n");
                start_index = stack_info.start_index + round(get(ui.slider, 'Value')) - 1;
                display_warning("go to end frame and press n");
                setFrame(stack_info.end_index-stack_info.start_index + 1);
                wait_for_keypress("n");
                end_index = stack_info.start_index + round(get(ui.slider, 'Value')) - 1;
            end
            shorten_slider(start_index, end_index);
            stack_info.shortened = true;
            toggle_indicator(ui.shortened_indicator, true);
            % save the start and end indices in a stack_info.mat file
            save_stack_callback();
        end
        function shorten_all_stack_callback(~,~)
            % iterate over all the stacks and align them
            for i = 1:length(stack_paths)
                set(ui.stack_dropdown, 'Value', i);
                load_images_callback();
                if isfield(stack_info, 'shortened')
                    if stack_info.shortened == true
                        logs{end+1} = sprintf("Trial %s already shortened",stack_paths(i));
                        continue;
                    end
                end
                if contains(path, 'time_control')
                    stack_info.shortened = true;
                    save_stack_callback();
                    continue;
                end
                shorten_stack_callback();
            end
        end
        function draw_masks_callback(~,~)
            % draw two right angle triangles on the image
            if isfield(stack_info, 'masked')
                if stack_info.masked == true
                    left_vertices = stack_info.mask_vertices(1:3, :);
                    right_vertices = stack_info.mask_vertices(4:6, :);
                end
            else
                left_vertices = [0 400;0 800;250 800];
                right_vertices = [800 400;800 800;550 800];
            end 

            left_triangle = drawpolygon(ui.ax1, 'Position', left_vertices);
            right_triangle = drawpolygon(ui.ax1, 'Position', right_vertices);
            display_warning("Press n to save the mask");
            wait_for_keypress("n");
            % get the mask from the drawn triangles
            mask = createMask(left_triangle) | createMask(right_triangle);
            % save the mask
            stack_info.mask = mask;
            stack_info.mask_vertices = [left_triangle.Position; right_triangle.Position];
            stack_info.masked = true;
            save_stack_callback();
            % clear the drawn triangles
            delete(left_triangle);
            delete(right_triangle);
            setFrame(get(ui.slider, 'Value'));
            fprintf('Masked stack %s\n', path);
        end
        function draw_all_masks_callback(~,~)
            % iterate over all the stacks and align them
            for i = 1:length(stack_paths)
                set(ui.stack_dropdown, 'Value', i);
                load_images_callback();
                if isfield(stack_info, 'masked')
                    if stack_info.masked == true
                        logs{end+1} = sprintf("Trial %s already masked",stack_paths(i));
                        continue;
                    end
                end
                if contains(path, 'time_control')
                    stack_info.masked = true;
                    save_stack_callback();
                    continue;
                end
                draw_masks_callback();
            end
        end
        function shorten_slider(start_index, end_index)
            % display_warning(sprintf("Shortening stack from %d to %d", start_index, end_index));
            % shorten stack
            stack_info.start_index = start_index;
            stack_info.end_index = end_index;
            set(ui.slider, 'Min', 1, 'Max', end_index - start_index + 1, 'Value', 1, ...
                'SliderStep', [1/(end_index-start_index) , 1/(end_index-start_index)], 'Visible', 'On');
            set(ui.to_frame, 'String', num2str(end_index - start_index + 1));
            setFrame(1);
        end
        function remove_shortening_callback(~, ~)
            % remove the displacements file
            stack_info.shortened = false;
            stack_info.start_index = 1;
            stack_info.end_index = numel(stack_info.img_data.img_files);
            shorten_slider(1, numel(stack_info.img_data.img_files));
            toggle_indicator(ui.shortened_indicator, false);
            % get current frame
            frame = get(ui.slider, 'Value');
            setFrame(frame);
            % add star to save button
            save_button.String = 'Save *';
        end
    end
end