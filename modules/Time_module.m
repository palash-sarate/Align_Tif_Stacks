classdef Time_module < handle
    properties
        stack_info
    end
    methods
        % constructor
        function obj = Time_module(stack_info)
           obj.stack_info = stack_info;
        end
        function get_times(~, ~)
            % set get times button to cancel get times that can cancel the operation
            get_time.String = 'Get time *';
            % get the time from every nth image
            n = 100;
            % create empty cell array to store results if stack_info doesnt have it already
            if ~isfield(stack_info, 'ocr_results')
                stack_info.ocr_results = cell(1, stack_info.img_data.num_imgs);            
            end
            for i = 1:n:stack_info.img_data.num_imgs
                % disp(i)
                % skip if stack_info.ocr_results{i} already exists
                % if ~isempty(stack_info.ocr_results{i})
                %     continue;
                % end  
                py_results = get_time_ocr(i);
                % display_warning(sprintf("frame %d processed",i));
                % convert results from python list to cell array
                results = pyList2cell(py_results);
                % add results to stack info 
                stack_info.ocr_results{i} = results;
                % assignin('base', 'ocr_results', stack_info.ocr_results);
            end
            save_stack_callback();
            get_time.String = 'Get time';
        end
        function timestamp = get_binary_timestamp(image_idx)
            % disp(image_idx);
            img_path = fullfile(stack_info.img_data.img_files(image_idx).folder, ...
                                stack_info.img_data.img_files(image_idx).name);
            % Try to decode the timestamp
            try
                timestamp = py.pco.decode_timestamp(img_path);
                % convert timestamp from py.dict to matlab struct
                timestamp = structfun(@double, struct(timestamp), 'UniformOutput', false);
                valid_timestamp = isValidTimeStamp(timestamp);         
            catch
                display_warning('Failed to decode timestamp. Not a valid binary timestamp image.');
            end
            
            if ~valid_timestamp
                timestamp = [];
                fprintf('Invalid timestamp found in image %d\n', image_idx);
            end
        end
        function timestamp_stack(~,~)
            if ~isfield(stack_info, 'timestamps')
                stack_info.timestamps = cell(1, stack_info.img_data.num_imgs);
            end
            WaitMessage = parfor_wait(stack_info.img_data.num_imgs, 'Waitbar', true);
            for i = 1:stack_info.img_data.num_imgs
                timestamp = get_binary_timestamp(i);
                if ~isempty(timestamp)
                    stack_info.timestamps{i} = timestamp;
                    WaitMessage.Send;
                else
                    % if one time stamp is invalid, then binary stamps don't exist
                    break;
                end
            end
            save_stack_callback();
            WaitMessage.Destroy;
        end
        function timestamp_all_stacks(~,~)
            WaitMessage = parfor_wait(length(stack_paths), 'Waitbar', true);
            for i = 1:length(stack_paths)
                set(ui.stack_dropdown, 'Value', i);
                load_images_callback();
                timestamp_stack();
                WaitMessage.Send;
            end
            WaitMessage.Destroy;
        end
        function ocr_results = get_time_ocr(k)
            % get the time from the image
            slider_idx = k;
            image_idx = stack_info.start_index + slider_idx - 1;
    
            img_path = fullfile(stack_info.img_data.img_files(image_idx).folder, stack_info.img_data.img_files(image_idx).name);
            % h = drawrectangle('Parent', ax1);wait(h);
            % roi = round(h.Position);
            % roi_python = int32([roi(1), roi(2), roi(3), roi(4)]);
            % get the time from the cropped image
            ocr_results = py.EasyOcr.ocr_from_file(img_path, []);
            % assignin('base', 'ocr_results', char(ocr_results));
            % assignin('base', 'roi', roi);
            % display_warning(ocr_results.Text);
        end
        function plot_all_timestamps(~,~)
            Ns = [];
            % check if struct with gr of all stacks exists at F:\shake_table_data\Results
            if exist('F:\shake_table_data\Results\timestamp_all_stacks.mat', 'file')
                timestamp_all_stacks = load('F:\shake_table_data\Results\timestamp_all_stacks.mat');
                while isfield(timestamp_all_stacks, 'timestamp_all_stacks')
                    timestamp_all_stacks = timestamp_all_stacks.timestamp_all_stacks;
                end
            else
                timestamp_all_stacks = struct();
            end
            % clear axis
            cla(ui.ax2);hold on;
            WaitMessage = parfor_wait(length(stack_paths), 'Waitbar', true);
            % iterate over all the stacks
            for i = 1:length(stack_paths)
                set(ui.stack_dropdown, 'Value', i);
                current_idx = get(ui.stack_dropdown, 'Value');
                path = stack_paths{current_idx};
                [iteration, parentDir] = getIteration(path);
                [N, fs] = get_info(path);
                Ns = [Ns, N];
                if contains(path, 'time_control') || contains(path, 'temp')
                    fprintf('Skipping %s\n', path);           
                    continue;
                end
                if isfield(timestamp_all_stacks, sprintf('N%d', N)) && ...
                        isfield(timestamp_all_stacks.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                        isfield(timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration))
                    fprintf('Timestamps found in timestamp_all_stacks for %s\n', path);
                    timestamps = timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('timestamps');
                    plot(ui.ax2, cellfun(@(t) t.time_us, timestamps) - timestamps{1}.time_us, 'Color', get_color(N, [4,12,24,48]), 'DisplayName', 'None');
                else
                    fprintf('Timestamps not found in timestamp_all_stacks for %s\n now checking stack_info\n', path);
                    if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                        load_stack_info();
                        if isfield(stack_info, 'timestamps') && isValidTimeStamp(stack_info.timestamps{1})
                            timestamps = stack_info.timestamps;
                            timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('timestamps') = timestamps;
                            fprintf('Timestamps found in stack_info for %s\n plotting in colour %s', path,mat2str(get_color(N, [4,12,24,48])));
                            plot(ui.ax2, cellfun(@(t) t.time_us, timestamps) - timestamps{1}.time_us, 'Color', get_color(N, [4,12,24,48]), 'DisplayName', 'None');
                        end
                    else
                        fprintf('Timestamps not found in stack_info for %s\n', path);
                    end
                end
                WaitMessage.Send;
            end
            % save the gr_all_stacks
            save('F:\shake_table_data\Results\timestamp_all_stacks.mat', 'timestamp_all_stacks');
            WaitMessage.Destroy;
        end
    end
end