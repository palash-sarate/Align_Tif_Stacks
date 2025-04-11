classdef Timer < handle
    properties
        app
    end
    methods
        function obj = Timer(app)
            obj.app = app;
        end
        %%%%%%%%%%%%%%%%%%%%%% TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function validity = isValidTimeStamp(~, timestamp)
            valid_years = [2023, 2024, 2025]; % List of valid years
            validity = isfield(timestamp, 'YYYY') && isfield(timestamp, 'MM') && ...
                isfield(timestamp, 'DD') && isfield(timestamp, 'h') && ...
                isfield(timestamp, 'min') && isfield(timestamp, 's') && ...
                isfield(timestamp, 'us');
            if validity
                validity = ismember(timestamp.YYYY, valid_years) && ...
                    timestamp.MM >= 1 && timestamp.MM <= 12 && ...
                    timestamp.DD >= 1 && timestamp.DD <= 31 && ...
                    timestamp.h >= 0 && timestamp.h <= 23 && ...
                    timestamp.min >= 0 && timestamp.min <= 59 && ...
                    timestamp.s >= 0 && timestamp.s <= 59 && ...
                    timestamp.us >= 0 && timestamp.us <= 999999;
            end
        end
        function timestamp_stack(obj, ~,~)
            if ~isfield(obj.app.stack_info, 'timestamps')
                obj.app.stack_info.timestamps = cell(1, obj.app.stack_info.img_data.num_imgs);
            end
            binary_timestamp_available = true;
            WaitMessage = parfor_wait(obj.app.stack_info.img_data.num_imgs, 'Waitbar', true);
            for i = 1:obj.app.stack_info.img_data.num_imgs
                if binary_timestamp_available
                    timestamp = get_binary_timestamp(i);
                    if ~isempty(timestamp)
                        obj.app.stack_info.timestamps{i} = timestamp;
                    else
                        binary_timestamp_available = false;
                        fprintf('Binary timestamp not available for stack %s\n', obj.app.path);
                    end 
                end
                if ~binary_timestamp_available
                    % if one time stamp is invalid, then binary stamps don't exist
                    % predict timestamp from the averaged time graph
                    timestamp = obj.predict_timeStamp(i);
                    obj.app.stack_info.timestamps{i} = timestamp;
                end
                WaitMessage.Send;
            end
            obj.app.utils.save_stack_callback();
            WaitMessage.Destroy;
            function timestamp = get_binary_timestamp(image_idx)
                % disp(image_idx);
                img_path = fullfile(obj.app.stack_info.img_data.img_files(image_idx).folder, ...
                                    obj.app.stack_info.img_data.img_files(image_idx).name);
                % Try to decode the timestamp
                try
                    timestamp = py.pco.decode_timestamp(img_path);
                    % convert timestamp from py.dict to matlab struct
                    timestamp = structfun(@double, struct(timestamp), 'UniformOutput', false);
                    valid_timestamp = obj.app.isValidTimeStamp(timestamp);         
                catch
                    obj.app.utils.display_warning('Failed to decode timestamp. Not a valid binary timestamp image.');
                end
                
                if ~valid_timestamp
                    timestamp = [];
                    fprintf('Invalid timestamp found in image %d\n', image_idx);
                end
            end
        end
        function timestamp_all_stacks = average_timeStamps(~)
            % check if struct with gr of all stacks exists at F:\shake_table_data\Results
            if exist('F:\shake_table_data\Results\timestamp_all_stacks.mat', 'file')
                timestamp_all_stacks = load('F:\shake_table_data\Results\timestamp_all_stacks.mat');
                while isfield(timestamp_all_stacks, 'timestamp_all_stacks')
                    timestamp_all_stacks = timestamp_all_stacks.timestamp_all_stacks;
                end
            else
                fprintf('No timestamp_all_stacks found Please run plot_all_timestamps\n');
            end
            % ignore list (N4, F14, Iter2)(N12, F14, Iter3)
            % loop over all the fields in timestamp_all_stacks
            % max_length = get_max_length_timestamps(timestamp_all_stacks);
            % average_timeStamp = struct();
            % N = fieldnames(timestamp_all_stacks)';
            % for i = numel(N)
            %     F = fieldnames(timestamp_all_stacks.(N{i}))';
            %     for j = numel(F)
            %         Iter = fieldnames(timestamp_all_stacks.(N{i}).(F{j}))';
            %         for k = numel(Iter)
            %             % timestamps = timestamp_all_stacks.(N{i}).(F{j}).(Iter{k}).timestamps;                   
            %         end
            %     end
            % end
            timestamp_all_stacks.('avg') = timestamp_all_stacks.('N24').('F14').('Iter4').timestamps;
            function get_max_length_timestamps(timestamp_all_stacks)
                % get the maximum number of timestamps
                max_timestamps = 0;
                N = fieldnames(timestamp_all_stacks)';
                for i = 1:numel(N)
                    F = fieldnames(timestamp_all_stacks.(N{i}))';
                    for j = 1:numel(F)
                        Iter = fieldnames(timestamp_all_stacks.(N{i}).(F{j}))';
                        for k = 1:numel(Iter)
                            timestamps = timestamp_all_stacks.(N{i}).(F{j}).(Iter{k}).timestamps;
                            if numel(timestamps) > max_timestamps
                                max_timestamps = numel(timestamps);
                            end                    
                        end
                    end
                end
            end
        end
        function timestamp_all_stacks(obj, ~,~)
            WaitMessage = parfor_wait(length(obj.app.stack_paths), 'Waitbar', true);
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                if contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                obj.app.load_images_callback();
                obj.app.timestamp_stack();
                WaitMessage.Send;
            end
            WaitMessage.Destroy;
        end
        % function to start getting the times for every nth image and also has the option to stop the current operation
        function get_times(obj, ~, ~)
            % set get times button to cancel get times that can cancel the operation
            obj.app.ui.controls.get_time.String = 'Get time *';
            % get the time from every nth image
            n = 100;
            % create empty cell array to store results if stack_info doesn't have it already
            if ~isfield(obj.app.stack_info, 'ocr_results')
                obj.app.stack_info.ocr_results = cell(1, obj.app.stack_info.img_data.num_imgs);            
            end
            for i = 1:n:obj.app.stack_info.img_data.num_imgs
                % disp(i)
                % skip if stack_info.ocr_results{i} already exists
                % if ~isempty(stack_info.ocr_results{i})
                %     continue;
                % end  
                py_results = get_time_ocr(i);
                % obj.app.utils.display_warning(sprintf("frame %d processed",i));
                % convert results from python list to cell array
                results = pyList2cell(py_results);
                % add results to stack info 
                obj.app.stack_info.ocr_results{i} = results;
                % assignin('base', 'ocr_results', stack_info.ocr_results);
            end
            obj.app.utils.save_stack_callback();
            obj.app.ui.controls.get_time.String = 'Get time';
            function ocr_results = get_time_ocr(k)
                % get the time from the image
                slider_idx = k;
                image_idx = obj.app.stack_info.start_index + slider_idx - 1;
        
                img_path = fullfile(obj.app.stack_info.img_data.img_files(image_idx).folder, obj.app.stack_info.img_data.img_files(image_idx).name);
                % h = drawrectangle('Parent', obj.app.ui.controls.ax1);wait(h);
                % roi = round(h.Position);
                % roi_python = int32([roi(1), roi(2), roi(3), roi(4)]);
                % get the time from the cropped image
                ocr_results = py.EasyOcr.ocr_from_file(img_path, []);
                % assignin('base', 'ocr_results', char(ocr_results));
                % assignin('base', 'roi', roi);
                % obj.app.utils.display_warning(ocr_results.Text);
            end
            function out = pyList2cell(pyobj)
                % Converts a python list/tuple of lists/tuples into a nested cell array.
                out = cell(pyobj); % Convert top-level list/tuple to a cell array
                for i = 1:numel(out)
                    if isa(out{i}, 'py.list') || isa(out{i}, 'py.tuple')
                        out{i} = pyList2cell(out{i}); % Recursively handle nested lists
                    elseif isa(out{i}, 'py.str')
                        % Convert Python string to MATLAB string
                        matlab_str = char(out{i});
                        
                        % Try to convert to number if possible
                        num_val = str2double(matlab_str);
                        if ~isnan(num_val)
                            out{i} = num_val; % Use number if conversion succeeded
                        else
                            out{i} = matlab_str; % Keep as string if not a number
                        end
                    end
                end
            end
        end
        function value = get_jammed_or_flows(~, N, f, iter)
            persistent data;
            % load data if its not already loaded
            if isempty(data)
                    xl_path = "F:\shake_table_data\Results\Results.xlsx";
                    % load the excel file
                    opts = detectImportOptions(xl_path,VariableNamesRange="A1:L1",VariableNamingRule="preserve");
                    opts.VariableTypes = repmat("string",1,numel(opts.VariableTypes));
                    data = readtable(xl_path, opts);
                    % data = data(:,[2,4,6,7]);
            end
            % get the jammed or flows from the data
            % get the value of 7th column where value of 2nd column is N, 4th column is f and 6th column is iter
            % Ensure N, f, and iter are strings to match the table column types
            value = data{strcmp(data{:, 2}, N) & strcmp(data{:, 3}, f) & strcmp(data{:, 6}, iter), 7};
            if value == "jammed"
                value = 1;
            elseif value == "flows"
                value = 0;
            else
                value = -1;
            end
        end
        function plot_all_timestamps(obj, ~,~)
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
            cla(obj.app.ui.controls.ax2);hold on;
            WaitMessage = parfor_wait(length(obj.app.stack_paths), 'Waitbar', true);
            % iterate over all the stacks
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                current_idx = get(obj.app.ui.controls.stackDropdown, 'Value');
                obj.app.path = obj.app.stack_paths{current_idx};
                [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
                [N, fs] = obj.app.utils.get_info(obj.app.path);
                fprintf("jammed : %d\n",jammed);
                Ns = [Ns, N];
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                if isfield(timestamp_all_stacks, sprintf('N%d', N)) && ...
                        isfield(timestamp_all_stacks.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                        isfield(timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration))
                    fprintf('Timestamps found in timestamp_all_stacks for %s\n', obj.app.path);
                    timestamps = timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('timestamps');
                    plot(obj.app.ui.controls.ax2, cellfun(@(t) t.time_us, timestamps) - timestamps{1}.time_us, 'Color', obj.app.utils.get_color(N), 'DisplayName', 'None');
                else
                    fprintf('Timestamps not found in timestamp_all_stacks for %s\n now checking stack_info\n', obj.app.path);
                    if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                        obj.app.utils.load_stack_info();
                        if isfield(obj.app.stack_info, 'timestamps') && obj.app.isValidTimeStamp(obj.app.stack_info.timestamps{1})
                            timestamps = obj.app.stack_info.timestamps;
                            timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('timestamps') = timestamps;
                            % fprintf('Timestamps found in stack_info for %s\n plotting in colour %s', path,mat2str(obj.app.utils.get_color(N)));
                            plot(obj.app.ui.controls.ax2, cellfun(@(t) t.time_us, timestamps) - timestamps{1}.time_us, 'Color', obj.app.utils.get_color(N), 'DisplayName', 'None');
                        end
                    else
                        fprintf('Timestamps not found in stack_info for %s\n', obj.app.path);
                    end
                end
                WaitMessage.Send;
            end
            hold off;
            title('Timestamps');
            xlabel('Frame number');
            % average the timestamps
            timestamp_all_stacks = obj.app.average_timeStamps();
            % save the gr_all_stacks
            save('F:\shake_table_data\Results\timestamp_all_stacks.mat', 'timestamp_all_stacks');
            % save the plot
            saveas(obj.app.ui.controls.ax2, 'F:\shake_table_data\Results\timestamps.png');
            WaitMessage.Destroy;
        end
        function timestamp = predict_timeStamp(~, image_idx)
            % check if struct with gr of all stacks exists at F:\shake_table_data\Results
            persistent timestamp_all_stacks;
            if evalin('base', 'exist(''timestamp_all_stacks'', ''var'')')
                timestamp_all_stacks = evalin('base', 'timestamp_all_stacks');
            else
                if exist('F:\shake_table_data\Results\timestamp_all_stacks.mat', 'file')
                    timestamp_all_stacks = load('F:\shake_table_data\Results\timestamp_all_stacks.mat');
                    while isfield(timestamp_all_stacks, 'timestamp_all_stacks')
                        timestamp_all_stacks = timestamp_all_stacks.timestamp_all_stacks;
                    end
                    assignin('base', 'timestamp_all_stacks', timestamp_all_stacks);
                else
                    fprintf('No timestamp_all_stacks found. Please run plot_all_timestamps\n');
                    return;
                end
            end
            timestamp = timestamp_all_stacks.('avg'){image_idx};
        end
        function get_stack_durations(obj)
            % loop over all stacks and get duration by 
            % subtracting the start_index time stamp from the end_index timestamp
            stack_durations = struct();
            WaitMessage = parfor_wait(length(obj.app.stack_paths), 'Waitbar', true);
            
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                [N, fs] = obj.app.utils.get_info(obj.app.path);
                [iteration, ~] = obj.app.utils.getIteration(obj.app.path);
                if contains(obj.app.path, 'temp') || contains(obj.app.path, 'time_control')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                obj.app.utils.load_stack_info();
                obj.app.stack_info.duration = obj.time_2_sec(obj.app.stack_info.timestamps{obj.app.stack_info.end_index}) - ...
                    obj.time_2_sec(obj.app.stack_info.timestamps{obj.app.stack_info.start_index});
                stack_durations.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)) = obj.app.stack_info.duration;
                fprintf('Duration for stack %s is %d',obj.app.path,obj.app.stack_info.duration);
                % obj.app.utils.save_stack_callback();
                % assignin('base', 'stack_info', stack_info);
                WaitMessage.Send;
            end
            save('F:\shake_table_data\Results\stack_durations.mat', 'stack_durations');
            WaitMessage.Destroy;
        end
        function plot_stack_durations(obj)
            % check if struct with gr of all stacks exists at F:\shake_table_data\Results
            if exist('F:\shake_table_data\Results\stack_durations.mat', 'file')
                stack_durations = load('F:\shake_table_data\Results\stack_durations.mat');
                while isfield(stack_durations, 'stack_durations')
                    stack_durations = stack_durations.stack_durations;
                end
            else
                fprintf('No stack_durations found. Please run get_stack_durations\n');
                obj.get_stack_durations();
                obj.plot_stack_durations();
                return;
            end
            % clear axis
            cla(obj.app.ui.controls.ax2);hold on;
            % iterate over all the stacks
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                current_idx = get(obj.app.ui.controls.stackDropdown, 'Value');
                obj.app.path = obj.app.stack_paths{current_idx};
                [iteration, ~] = obj.app.utils.getIteration(obj.app.path);
                [N, fs] = obj.app.utils.get_info(obj.app.path);
                jammed = obj.get_jammed_or_flows(sprintf("%d",N),sprintf("%d",fs),iteration);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                if isfield(stack_durations, sprintf('N%d', N)) && ...
                        isfield(stack_durations.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                        isfield(stack_durations.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration))
                    fprintf('Duration found in stack_durations for %s\n', obj.app.path);
                    duration = stack_durations.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration));
                    scatter(obj.app.ui.controls.ax2, fs, duration, obj.app.utils.get_marker(N), 'MarkerFaceColor', obj.app.utils.get_color(N), 'MarkerEdgeColor', obj.app.utils.get_color(N),'HandleVisibility','off');
                    if jammed == 1
                        % add transparent yellow circle at same location
                        scatter(obj.app.ui.controls.ax2, fs, duration,40, obj.app.utils.get_marker(N), 'MarkerEdgeColor', 'yellow', 'LineWidth', 2, 'HandleVisibility', 'off');
                    end
                else
                    fprintf('Duration not found in stack_durations for %s\n', obj.app.path);
                end
            end
            unique_Ns = [4,12,24,48];
            legend_handles = zeros(1, numel(unique_Ns));
            legend_entries = cell(1, numel(unique_Ns));
            
            for j = 1:numel(unique_Ns)
                N_val = unique_Ns(j);
                % Create a "dummy" line just for the legend with the right color
                legend_handles(j) = scatter(obj.app.ui.controls.ax2, NaN, NaN, obj.app.utils.get_marker(N_val), 'MarkerFaceColor', obj.app.utils.get_color(N_val), 'MarkerEdgeColor',obj.app.utils.get_color(N_val));
                legend_entries{j} = sprintf('N = %d', N_val);
            end
            
            % Create legend using only our dummy lines
            legend(obj.app.ui.controls.ax2, legend_handles, legend_entries, 'Location', 'best');
            hold off;
            xlim([0, 22]);
            ylim([0, 1200]);
            title('Stack durations');
            xlabel('Frequency (Hz)');
            ylabel('Duration (s)');
            exportgraphics(obj.app.ui.controls.ax2, 'F:\shake_table_data\Results\stack_durations.png');
        end
        function secs = time_2_sec(~, timestamp)
            secs = timestamp.h * 3600 + timestamp.min * 60 + timestamp.s + timestamp.us / 1e6;
        end
    end
end
