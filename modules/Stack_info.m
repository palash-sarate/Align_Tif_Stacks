classdef StackInfo < handle
    properties
        img_data
        particle_locations
        distances
        gr
        gr_bins
        aligned = false
        shortened = false
    end
    
    methods
        function save(obj, parentDir, iteration)
            % Save the stack_info to a file
            stack_info = obj; % Create a copy with the variable name expected in the file
            save(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'stack_info');
        end
        function save_stack_callback(~,~)
            % save the current stack
            if stack_info.img_data.num_imgs == 0
                display_warning("load some images first");
            else
                while isfield(stack_info, 'stack_info')
                    stack_info = stack_info.stack_info;
                end
                temp_stack_info = stack_info;
                % get the current stack path
                path = stack_paths{get(ui.stack_dropdown, 'Value')};
                [iteration, parentDir] = getIteration(path);
                % remove the image data from the stack_info
                stack_info.img_data.imgs = cell(1, stack_info.img_data.num_imgs);
                % save the stack
                save(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'stack_info');
                display_warning(sprintf("Stack saved to %s//stack_info_%s.mat", parentDir, iteration));
                save_button.String = 'Save';
                stack_info = temp_stack_info;
                assignin('base', 'stack_info', stack_info);
            end
        end
        function stack_info = initialize_stack_info(img_data)
            % Image data
            path = stack_paths{get(ui.stack_dropdown, 'Value')};
            [iteration, parentDir] = getIteration(path);
            display_warning("Initializing stack info");
            stack_info = struct('start_index', 1, 'end_index', numel(img_data.img_files),...
                'parentDir', parentDir, 'iteration', iteration, ...
                'aligned', false, 'shortened', false, 'masked', false,...
                'displacements', zeros(numel(img_data.img_files), 2), 'img_data', img_data);
        end
    end
end