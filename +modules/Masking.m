classdef Masking < handle
    properties
        app
    end
    methods
        function obj = Masking(app)
            obj.app = app;
        end
        %%%%%%%%%%%%%%%%%%%%%% MASKING %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function draw_masks_callback(obj, ~,~)
            left_vertices = [0 400; 0 800; 250 800];
            right_vertices = [800 400; 800 800; 550 800];
            % draw two right angle triangles on the image
            if isfield(obj.app.stack_info, 'masked')
                if obj.app.stack_info.masked == true
                    left_vertices = obj.app.stack_info.mask_vertices(1:3, :);
                    right_vertices = obj.app.stack_info.mask_vertices(4:6, :);
                end
            end 

            left_triangle = drawpolygon(obj.app.ui.controls.ax1, 'Position', left_vertices);
            right_triangle = drawpolygon(obj.app.ui.controls.ax1, 'Position', right_vertices);
            obj.app.utils.display_warning("Press n to save the mask");
            obj.app.utils.wait_for_keypress("n");
            % get the mask from the drawn triangles
            mask = createMask(left_triangle) | createMask(right_triangle);
            % save the mask
            obj.app.stack_info.mask = mask;
            obj.app.stack_info.mask_vertices = [left_triangle.Position; right_triangle.Position];
            obj.app.stack_info.masked = true;
            obj.app.utils.save_stack_callback();
            % clear the drawn triangles
            delete(left_triangle);
            delete(right_triangle);
            obj.app.utils.setFrame(get(obj.app.ui.controls.slider, 'Value'));
            fprintf('Masked stack %s\n', obj.app.path);
        end
        function draw_all_masks_callback(obj, ~,~)
            % iterate over all the stacks and align them
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.load_images_callback();
                if isfield(obj.app.stack_info, 'masked')
                    if obj.app.stack_info.masked == true
                        obj.app.logs{end+1} = sprintf("Trial %s already masked",obj.app.stack_paths(i));
                        continue;
                    end
                end
                if contains(obj.app.path, 'time_control')
                    obj.app.stack_info.masked = true;
                    obj.app.utils.save_stack_callback();
                    continue;
                end
                obj.draw_masks_callback();
            end
        end
        function shorten_slider(obj, start_index, end_index)
            % obj.app.utils.display_warning(sprintf("Shortening stack from %d to %d", start_index, end_index));
            % shorten stack
            obj.app.stack_info.start_index = start_index;
            obj.app.stack_info.end_index = end_index;
            set(obj.app.ui.controls.slider, 'Min', 1, 'Max', end_index - start_index + 1, 'Value', 1, ...
                'SliderStep', [1/(end_index-start_index) , 1/(end_index-start_index)], 'Visible', 'On');
            set(obj.app.ui.info.toFrame, 'String', num2str(end_index - start_index + 1));
            obj.app.utils.setFrame(1);
        end
        function change_drive_callback(obj, ~, ~)
            % change the drive letter
            current_drive = 'E:';
            new_drive = 'F:';
            WaitMessage = parfor_wait(length(obj.app.stack_paths), 'Waitbar', true);
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                current_idx = get(obj.app.ui.controls.stackDropdown, 'Value');
                obj.app.path = obj.app.stack_paths{current_idx};
                obj.app.utils.update_info(obj.app.path);
                [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
                if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                    obj.app.stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
                    while isfield(obj.app.stack_info, 'stack_info')
                        obj.app.stack_info = obj.app.stack_info.stack_info;
                    end
                    % replace the parentDir with the new drive letter
                    obj.app.stack_info.parentDir = strrep(obj.app.stack_info.parentDir, current_drive, new_drive);
                    % replace the drive letter in each stack_info.img_data.img_files.folder
                    for j = 1:numel(obj.app.stack_info.img_data.img_files)
                        obj.app.stack_info.img_data.img_files(j).folder = strrep(obj.app.stack_info.img_data.img_files(j).folder, current_drive, new_drive);
                    end
                    obj.app.utils.save_stack_callback();
                    % assignin('base', 'new_stack_info', stack_info);
                else
                    continue;
                end
                WaitMessage.Send;
            end
            WaitMessage.Destroy
        end
    end
end