% TODO: Time stamp
% TODO: Provision to combine stacks
% TODO: Time stamp from OCR
% TODO: plot phase space from the total time taken for the trial
% TODO:  

classdef App < handle
    properties
        % Add properties here
        ui
        utils
        path
        stack_info
        stack_paths
        searchWindow
        particle_locations_visible
        is_playing
        forced
        skip_alignment
        speed
        logs
        monitorChoice = 1;
        speeds = {1,2,4,8};
    end

    methods
        function run(app)
            % Erase all existing variables.
            if count(py.sys.path,pwd) == 0
                insert(py.sys.path,int32(0),pwd);
            end
            setenv('TCL_LIBRARY', 'C:/Python312/tcl/tcl8.6');
            setenv('TK_LIBRARY', 'C:/Python312/tcl/tk8.6');
            % workspace;  % Make sure the workspace panel is showing.
            % start parallel pool
            % parpool(4);
            app.particle_locations_visible = false;
            app.is_playing = false;
            app.forced = false;
            app.skip_alignment = false;
            app.speed = 1;
            app.logs = {}; % Initialize logs list

            app.searchWindow = 50;
            app.path = [];
            app.stack_info = struct();
            app.stack_paths = app.get_stack_paths();
            app.utils = modules.Utils();
            app.ui = modules.Ui(app);
        end

        function stack_paths = get_stack_paths(~)
            % directory of the tif stacks
            % folder = 'F:\shake_table_data\';
            % populate the list of paths to the tiff stacks
            Ns = [4,12,24,48];
            fs = [4,6,8,10,12,14,16,18,20];
            deg = 60;
            wd = 10;
            stack_paths = [];
            % fps = 47;

            for n = 1:length(Ns)
                for freq = 1:length(fs)
                    for w = 1:length(wd)
                        % Construct the path with escaped backslashes
                        c_path = sprintf("F:\\shake_table_data\\N%d\\%dhz_hopperflow\\%ddeg\\%dcm\\",Ns(n),fs(freq),deg,wd);
                        % find folders in the path directory
                        subFolders = dir(c_path);
                        subFolders = subFolders([subFolders.isdir]);  % Keep only directories
                        subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % Remove '.' and '..' directories

                        for k = 1:length(subFolders)
                            stack_paths = [stack_paths; fullfile(c_path, subFolders(k).name)];
                        end
                    end
                end
            end

            stack_paths = [stack_paths; "F:\shake_table_data\time_control\\1"];
            stack_paths = [stack_paths; "F:\shake_table_data\time_control\\2"];
            stack_paths = [stack_paths; "F:\shake_table_data\time_control\\3"];
        end

        function load_images_callback(obj, ~, ~)
            WaitMessage = parfor_wait(4, 'Waitbar', true);
            current_idx = get(obj.ui.controls.stackDropdown, 'Value');
            obj.path = obj.stack_paths{current_idx};
            % if path has time_control in it load the get_times button
            obj.toggle_get_time_ui();

            obj.update_info(obj.path);
            [iteration, parentDir] = obj.getIteration(obj.path);
            set(obj.ui.info.stackLabel, 'String', sprintf('Stack #%d of %d', current_idx, length(obj.stack_paths)));

            if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                obj.stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
                while isfield(obj.stack_info, 'stack_info')
                    obj.stack_info = obj.stack_info.stack_info;
                end
                WaitMessage.Send;
                assignin('base', 'stack_info', obj.stack_info);
            else
                img_data.img_files = dir(fullfile(obj.path, '*.tif'));
                obj.stack_info = obj.initialize_stack_info(img_data);
            end
            
            obj.stack_info.img_data.num_imgs = numel(obj.stack_info.img_data.img_files);
            obj.stack_info.img_data.imgs = cell(1, obj.stack_info.img_data.num_imgs);
            WaitMessage.Send;
            % if start and end indices are set, set the shortened indicator to green
            if obj.stack_info.shortened == true
                obj.toggle_indicator(obj.ui.info.shortenedIndicator, true);
            else
                obj.toggle_indicator(obj.ui.info.shortenedIndicator, false);
            end
            WaitMessage.Send;
            % if displacement_n.mat exists, set the aligned indicator to green
            if obj.stack_info.aligned == true && ~isempty(obj.stack_info.displacements)
                obj.toggle_indicator(obj.ui.info.alignedIndicator, true);
                % plot the displacements on obj.ui.controls.ax2
                obj.plot_displacements();
            else
                obj.toggle_indicator(obj.ui.info.alignedIndicator, false);
                % clear axis
                cla(obj.ui.controls.ax2);
            end
            obj.shorten_slider(obj.stack_info.start_index, obj.stack_info.end_index);
            WaitMessage.Send;
            fprintf('Loaded stack %s\n', obj.path);
            WaitMessage.Destroy;
        end
    %%%%%%%%%%%%%%%%%%%%%% SCALE %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function get_scales(obj, ~,~)
            num_particles_to_check = 5;
            for i = 1:length(obj.stack_paths)
                if contains(obj.stack_paths{i}, 'time_control') || contains(obj.stack_paths{i}, 'temp') || contains(obj.stack_paths{i}, 'cont')
                    continue;
                end
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.load_images_callback();
                if isfield(obj.stack_info, 'scale')
                    fprintf('Scales already exist for stack %s\n', obj.path);
                    if ~obj.forced
                        fprintf('Skipping stack %s\n', obj.path);
                        continue;
                    end
                end
                get_scale(num_particles_to_check);
                % average the obj.stack_info.radii
                radii_values = cellfun(@(x) x{2}, obj.stack_info.radii);
                avg_radius = mean(radii_values);
                std_radius = std(radii_values);
                obj.stack_info.scale = 2 * avg_radius;
                obj.stack_info.scale_std = std_radius;
                obj.save_stack_callback();
            end
            function get_scale(n)
                % set the viewer to first image
                obj.setFrame(1);
                fprintf('Getting scales for stack %s\n', obj.path);
                % get the particle locations
                particle_locations = obj.stack_info.particle_locations;
                % zoom onto 10 particles selected randomly from the image
                for i = 1:n
                    random_index = randi(size(particle_locations, 1));
                    % get the particle location
                    x = particle_locations.x(random_index);
                    y = particle_locations.y(random_index);
                    radius = particle_locations.size(random_index) * 3;
                    set(gcf, 'WindowScrollWheelFcn', {@zoom_callback,x,y});
                    % zoom onto the particle
                    xlim(obj.ui.controls.ax1, [x - 50, x + 50]);
                    ylim(obj.ui.controls.ax1, [y - 50, y + 50]);
        
                    % highlight the particle
                    hold(obj.ui.controls.ax1, 'on');
                    plot(obj.ui.controls.ax1, x, y, 'r*');
                    hold(obj.ui.controls.ax1, 'off');
                    % get the scale from the user by asking to draw a circle around the particle
                    h = drawcircle(obj.ui.controls.ax1, 'Center', [x, y], 'Radius', radius);
                    wait(h);
                    % save the circle to the stack_info
                    % add the scale to the stack_info
                    if ~isfield(obj.stack_info, 'radii')
                        obj.stack_info.radii = {};
                    end
                    obj.stack_info.radii{end+1} = {h.Center, h.Radius};
                end
                % delete the circle
                delete(h);
                % reset the zoom
                % zoom(obj.ui.controls.ax1,'reset');
                % enable the scroll callback
                set(gcf, 'WindowScrollWheelFcn', @obj.scrollWheelMoved);
                % assignin('base', 'stack_info', obj.stack_info);
            end
        end
        function plot_scales(obj)
            scales = zeros(1, length(obj.stack_paths));
            sizes = zeros(1, length(obj.stack_paths));
            scales_std = zeros(1, length(obj.stack_paths));
            sizes_std = zeros(1, length(obj.stack_paths));
            % check if csv already exists
            if isfile('F:\shake_table_data\Results\scales.csv')
                data = readmatrix('F:\shake_table_data\Results\scales.csv');
                scales = data(:,1);
                scales_std = data(:,2);
                sizes = data(:,3);
                sizes_std = data(:,4);
            else
                for i = 1:length(obj.stack_paths)
                    if contains(obj.stack_paths{i}, 'time_control') || contains(obj.stack_paths{i}, 'temp') || contains(obj.stack_paths{i}, 'cont')
                        continue;
                    end
                    set(obj.ui.controls.stackDropdown, 'Value', i);
                    obj.load_images_callback();
                    if isfield(obj.stack_info, 'scale')                
                        scales(i) = obj.stack_info.scale;
                        scales_std(i) = obj.stack_info.scale_std;
                    end
                    sizes(i) = mean(obj.stack_info.particle_locations.size);
                    sizes_std(i) = std(obj.stack_info.particle_locations.size);
                end
            end
            cla(obj.ui.controls.ax2);

            errorbar(obj.ui.controls.ax2, scales,scales_std, 'b', 'LineWidth', 2);
            hold(obj.ui.controls.ax2, 'on');
            errorbar(obj.ui.controls.ax2, sizes * 3 * 2,sizes_std, 'g', 'LineWidth', 2);
            hold(obj.ui.controls.ax2, 'off');
            title('Scales and particle sizes');
            xlabel('Stack number');
            ylabel('Scale');
            % save data to csv
            data = [scales; scales_std; sizes; sizes_std];
            data = data';
            writematrix(data, 'F:\shake_table_data\Results\scales.csv');
            saveas(gcf, 'F:\shake_table_data\Results\scales.png');
        end
    %%%%%%%%%%%%%%%%%%%%%% TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function validity = isValidTimeStamp(obj, timestamp)
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
            if ~isfield(obj.stack_info, 'timestamps')
                obj.stack_info.timestamps = cell(1, obj.stack_info.img_data.num_imgs);
            end
            binary_timestamp_available = true;
            WaitMessage = parfor_wait(obj.stack_info.img_data.num_imgs, 'Waitbar', true);
            for i = 1:obj.stack_info.img_data.num_imgs
                if binary_timestamp_available
                    timestamp = get_binary_timestamp(i);
                    if ~isempty(timestamp)
                        obj.stack_info.timestamps{i} = timestamp;
                    else
                        binary_timestamp_available = false;
                        fprintf('Binary timestamp not available for stack %s\n', obj.path);
                    end 
                end
                if ~binary_timestamp_available
                    % if one time stamp is invalid, then binary stamps don't exist
                    % predict timestamp from the averaged time graph
                    timestamp = obj.predict_timeStamp(i);
                    obj.stack_info.timestamps{i} = timestamp;
                end
                WaitMessage.Send;
            end
            obj.save_stack_callback();
            WaitMessage.Destroy;
            function timestamp = get_binary_timestamp(image_idx)
                % disp(image_idx);
                img_path = fullfile(obj.stack_info.img_data.img_files(image_idx).folder, ...
                                    obj.stack_info.img_data.img_files(image_idx).name);
                % Try to decode the timestamp
                try
                    timestamp = py.pco.decode_timestamp(img_path);
                    % convert timestamp from py.dict to matlab struct
                    timestamp = structfun(@double, struct(timestamp), 'UniformOutput', false);
                    valid_timestamp = obj.isValidTimeStamp(timestamp);         
                catch
                    obj.display_warning('Failed to decode timestamp. Not a valid binary timestamp image.');
                end
                
                if ~valid_timestamp
                    timestamp = [];
                    fprintf('Invalid timestamp found in image %d\n', image_idx);
                end
            end
        end
        function timestamp_all_stacks = average_timeStamps(obj)
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
            WaitMessage = parfor_wait(length(obj.stack_paths), 'Waitbar', true);
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.path = obj.stack_paths{i};
                if contains(obj.path, 'temp')
                    fprintf('Skipping %s\n', obj.path);           
                    continue;
                end
                obj.load_images_callback();
                obj.timestamp_stack();
                WaitMessage.Send;
            end
            WaitMessage.Destroy;
        end
        % function to start getting the times for every nth image and also has the option to stop the current operation
        function get_times(obj, ~, ~)
            % set get times button to cancel get times that can cancel the operation
            obj.ui.controls.get_time.String = 'Get time *';
            % get the time from every nth image
            n = 100;
            % create empty cell array to store results if stack_info doesn't have it already
            if ~isfield(obj.stack_info, 'ocr_results')
                obj.stack_info.ocr_results = cell(1, obj.stack_info.img_data.num_imgs);            
            end
            for i = 1:n:obj.stack_info.img_data.num_imgs
                % disp(i)
                % skip if stack_info.ocr_results{i} already exists
                % if ~isempty(stack_info.ocr_results{i})
                %     continue;
                % end  
                py_results = get_time_ocr(i);
                % obj.display_warning(sprintf("frame %d processed",i));
                % convert results from python list to cell array
                results = pyList2cell(py_results);
                % add results to stack info 
                obj.stack_info.ocr_results{i} = results;
                % assignin('base', 'ocr_results', stack_info.ocr_results);
            end
            obj.save_stack_callback();
            obj.ui.controls.get_time.String = 'Get time';
            function ocr_results = get_time_ocr(k)
                % get the time from the image
                slider_idx = k;
                image_idx = obj.stack_info.start_index + slider_idx - 1;
        
                img_path = fullfile(obj.stack_info.img_data.img_files(image_idx).folder, obj.stack_info.img_data.img_files(image_idx).name);
                % h = drawrectangle('Parent', obj.ui.controls.ax1);wait(h);
                % roi = round(h.Position);
                % roi_python = int32([roi(1), roi(2), roi(3), roi(4)]);
                % get the time from the cropped image
                ocr_results = py.EasyOcr.ocr_from_file(img_path, []);
                % assignin('base', 'ocr_results', char(ocr_results));
                % assignin('base', 'roi', roi);
                % obj.display_warning(ocr_results.Text);
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
        function value = get_jammed_or_flows(obj, N, f, iter)
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
            cla(obj.ui.controls.ax2);hold on;
            WaitMessage = parfor_wait(length(obj.stack_paths), 'Waitbar', true);
            % iterate over all the stacks
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                current_idx = get(obj.ui.controls.stackDropdown, 'Value');
                obj.path = obj.stack_paths{current_idx};
                [iteration, parentDir] = obj.getIteration(obj.path);
                [N, fs] = obj.get_info(obj.path);
                fprintf("jammed : %d\n",jammed);
                Ns = [Ns, N];
                if contains(obj.path, 'time_control') || contains(obj.path, 'temp')
                    fprintf('Skipping %s\n', obj.path);           
                    continue;
                end
                if isfield(timestamp_all_stacks, sprintf('N%d', N)) && ...
                        isfield(timestamp_all_stacks.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                        isfield(timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration))
                    fprintf('Timestamps found in timestamp_all_stacks for %s\n', obj.path);
                    timestamps = timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('timestamps');
                    plot(obj.ui.controls.ax2, cellfun(@(t) t.time_us, timestamps) - timestamps{1}.time_us, 'Color', obj.get_color(N), 'DisplayName', 'None');
                else
                    fprintf('Timestamps not found in timestamp_all_stacks for %s\n now checking stack_info\n', obj.path);
                    if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                        obj.load_stack_info();
                        if isfield(obj.stack_info, 'timestamps') && obj.isValidTimeStamp(obj.stack_info.timestamps{1})
                            timestamps = obj.stack_info.timestamps;
                            timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('timestamps') = timestamps;
                            % fprintf('Timestamps found in stack_info for %s\n plotting in colour %s', path,mat2str(obj.get_color(N)));
                            plot(obj.ui.controls.ax2, cellfun(@(t) t.time_us, timestamps) - timestamps{1}.time_us, 'Color', obj.get_color(N), 'DisplayName', 'None');
                        end
                    else
                        fprintf('Timestamps not found in stack_info for %s\n', obj.path);
                    end
                end
                WaitMessage.Send;
            end
            hold off;
            title('Timestamps');
            xlabel('Frame number');
            % average the timestamps
            timestamp_all_stacks = obj.average_timeStamps();
            % save the gr_all_stacks
            save('F:\shake_table_data\Results\timestamp_all_stacks.mat', 'timestamp_all_stacks');
            % save the plot
            saveas(obj.ui.controls.ax2, 'F:\shake_table_data\Results\timestamps.png');
            WaitMessage.Destroy;
        end
        function timestamp = predict_timeStamp(obj, image_idx)
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
            WaitMessage = parfor_wait(length(obj.stack_paths), 'Waitbar', true);
            
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.path = obj.stack_paths{i};
                [N, fs] = obj.get_info(obj.path);
                [iteration, ~] = obj.getIteration(obj.path);
                if contains(obj.path, 'temp') || contains(obj.path, 'time_control')
                    fprintf('Skipping %s\n', obj.path);           
                    continue;
                end
                obj.load_stack_info();
                obj.stack_info.duration = obj.time_2_sec(obj.stack_info.timestamps{obj.stack_info.end_index}) - ...
                    obj.time_2_sec(obj.stack_info.timestamps{obj.stack_info.start_index});
                stack_durations.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)) = obj.stack_info.duration;
                fprintf('Duration for stack %s is %d',obj.path,obj.stack_info.duration);
                % obj.save_stack_callback();
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
            cla(obj.ui.controls.ax2);hold on;
            % iterate over all the stacks
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                current_idx = get(obj.ui.controls.stackDropdown, 'Value');
                obj.path = obj.stack_paths{current_idx};
                [iteration, ~] = obj.getIteration(obj.path);
                [N, fs] = obj.get_info(obj.path);
                jammed = obj.get_jammed_or_flows(sprintf("%d",N),sprintf("%d",fs),iteration);
                if contains(obj.path, 'time_control') || contains(obj.path, 'temp')
                    fprintf('Skipping %s\n', obj.path);           
                    continue;
                end
                if isfield(stack_durations, sprintf('N%d', N)) && ...
                        isfield(stack_durations.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                        isfield(stack_durations.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration))
                    fprintf('Duration found in stack_durations for %s\n', obj.path);
                    duration = stack_durations.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration));
                    scatter(obj.ui.controls.ax2, fs, duration, obj.get_marker(N), 'MarkerFaceColor', obj.get_color(N), 'MarkerEdgeColor', obj.get_color(N),'HandleVisibility','off');
                    if jammed == 1
                        % add transparent yellow circle at same location
                        scatter(obj.ui.controls.ax2, fs, duration,40, obj.get_marker(N), 'MarkerEdgeColor', 'yellow', 'LineWidth', 2, 'HandleVisibility', 'off');
                    end
                else
                    fprintf('Duration not found in stack_durations for %s\n', obj.path);
                end
            end
            unique_Ns = [4,12,24,48];
            legend_handles = zeros(1, numel(unique_Ns));
            legend_entries = cell(1, numel(unique_Ns));
            
            for j = 1:numel(unique_Ns)
                N_val = unique_Ns(j);
                % Create a "dummy" line just for the legend with the right color
                legend_handles(j) = scatter(obj.ui.controls.ax2, NaN, NaN, obj.get_marker(N_val), 'MarkerFaceColor', obj.get_color(N_val), 'MarkerEdgeColor',obj.get_color(N_val));
                legend_entries{j} = sprintf('N = %d', N_val);
            end
            
            % Create legend using only our dummy lines
            legend(obj.ui.controls.ax2, legend_handles, legend_entries, 'Location', 'best');
            hold off;
            xlim([0, 22]);
            ylim([0, 1200]);
            title('Stack durations');
            xlabel('Frequency (Hz)');
            ylabel('Duration (s)');
            exportgraphics(obj.ui.controls.ax2, 'F:\shake_table_data\Results\stack_durations.png');
        end
        function secs = time_2_sec(~, timestamp)
            secs = timestamp.h * 3600 + timestamp.min * 60 + timestamp.s + timestamp.us / 1e6;
        end

    %%%%%%%%%%%%%%%%%%%%%% TRIAL CHARACTERIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
        function is_late = is_rec_start_late(obj, n,fs,iter)
            % ignore list (48,10,2)(12,14,3)(24,12,1) - first frame is false
            % return true is n,fs,iter is in the ignore list
            ignore_list = {[48,10,"2"],[12,14,"3"],[24,12,"1"]};
            is_late = false;
            % fprintf('Checking %d,%d,%s\n',n,fs,iter);
            if any(cellfun(@(x) isequal([n, fs, iter], x), ignore_list))
                % fprintf('Ignoring %d,%d,%s\n',n,fs,iter);
                is_late = true;
            end
        end
        function set_all_empty_or_not(obj)
            % function to loop over all stack and set if the stack emptied fully or not
            % should load the stack and go to the last frame and provide option to user to 
            % state whether the hopper has emptied or not
            % loop over all stacks and set the emptyornot field in stack_info
            % WaitMessage = parfor_wait(length(obj.stack_paths), 'Waitbar', true);
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.path = obj.stack_paths{i};
                % [N, fs] = obj.get_info(path);
                % [iteration, ~] = obj.getIteration(path);
                if contains(obj.path, 'time_control') || contains(obj.path, 'temp')
                    fprintf('Skipping %s\n', obj.path);           
                    continue;
                end
                obj.load_images_callback();
                obj.setFrame(obj.stack_info.end_index - obj.stack_info.start_index + 1);
                % show dialog box to user to set the emptyornot field
                options = {'Empty', 'Not Empty'};
                answer = questdlg('Is the hopper empty?', 'Empty or Not', options{1}, options{2}, options{1});
                if strcmp(answer, options{1})
                    obj.stack_info.empty = true;
                else
                    obj.stack_info.empty = false;
                end
                obj.save_stack_callback();
                % WaitMessage.Send;
            end
            % WaitMessage.Destroy;
        end
        function set_all_jam_or_not(obj)
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.path = obj.stack_paths{i};
                [N, ~] = obj.get_info(obj.path);
                % [iteration, ~] = obj.getIteration(path);
                if contains(obj.path, 'time_control') || contains(obj.path, 'temp')
                    fprintf('Skipping %s\n', obj.path);           
                    continue;
                end
                obj.load_images_callback();
                obj.setFrame(obj.stack_info.end_index - obj.stack_info.start_index + 1);
                if N == 4
                    obj.stack_info.jammed = true;
                else
                    % show dialog box to user to set the emptyornot field
                    options = {'Jammed', 'Not Jammed'};
                    answer = questdlg('Is the Jammed?', 'Jammed or Not', options{1}, options{2}, options{1});
                    if strcmp(answer, options{1})
                        obj.stack_info.jammed = true;
                    else
                        obj.stack_info.jammed = false;
                    end
                end
                obj.save_stack_callback();
                % WaitMessage.Send;
            end
            % WaitMessage.Destroy;
        end
    %%%%%%%%%%%%%%%%%%%%%% LOCAL ORDER PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%%
        function psi = localOrderParameter(obj, pos, m, r_cut)
            % obj.localOrderParameter computes the local bond-orientational order parameter.
            %
            %   psi = obj.localOrderParameter(pos, m) computes the order parameter for
            %   each particle given in pos (an N-by-2 array of [x,y] positions) using the
            %   symmetry index m. For instance, use m = 4 for square symmetry or m = 6 for
            %   triangular (hexagonal) symmetry.
            %
            %   psi = obj.localOrderParameter(pos, m, r_cut) only considers neighbors that are
            %   within the distance r_cut.
            %
            %   The local order parameter for particle i is defined as:
            %       psi_m(i) = | (1/N_i) * sum_{j in neighbors(i)} exp( i*m*theta_ij ) |
            %   where theta_ij is the angle between the vector (pos(j,:) - pos(i,:)) and the x-axis.
            %
            %   Example:
            %       % pos: N-by-2 array of particle coordinates
            %       psi_top = obj.localOrderParameter(pos(pos(:,2)>y_thresh,:), 4);
            %       psi_bottom = obj.localOrderParameter(pos(pos(:,2)<=y_thresh,:), 6);
            %
            %   See also delaunayTriangulation, atan2.
            %
            
            if nargin < 2
                error('You must provide at least the positions and symmetry parameter m.');
            end
            
            % Optional neighbor cutoff: if not provided, use all neighbors from the triangulation.
            if nargin < 3
                useCutoff = false;
            else
                useCutoff = true;
            end
            
            % Create a Delaunay triangulation from the particle positions.
            dt = delaunayTriangulation(pos(:,1), pos(:,2));
            
            % For each particle, find the indices of attached triangles (neighbors)
            attachList = vertexAttachments(dt);
            
            N = size(pos,1);
            psi = zeros(N,1);
            
            for i = 1:N
                % Extract the indices of triangles attached to particle i
                triIndices = attachList{i};
                % Get all vertices from these triangles
                nb = unique(dt.ConnectivityList(triIndices,:));
                % Remove the particle itself from its neighbor list
                nb(nb == i) = [];
                
                % If a cutoff distance is provided, filter the neighbors by distance.
                if useCutoff && ~isempty(nb)
                    distances = sqrt(sum((pos(nb,:) - pos(i,:)).^2, 2));
                    nb = nb(distances <= r_cut);
                end
                
                % If no neighbors are found, set order parameter to NaN.
                if isempty(nb)
                    psi(i) = NaN;
                    continue;
                end
                
                % Calculate the angle between particle i and each of its neighbors.
                angles = atan2(pos(nb,2) - pos(i,2), pos(nb,1) - pos(i,1));
                
                % Compute the local order parameter for particle i.
                psi(i) = abs(sum(exp(1i*m*angles)) / numel(angles));
            end 
        end
        function psi = freudLocalOrderParameter(obj, particles_path, bead_dia, l, nhood)
            psi = py.track.compute_ql(particles_path, bead_dia, int32(l), int32(nhood));
            % convert psi from python ndarray to matlab array
            psi = double(py.array.array('d', py.numpy.nditer(psi)));
        end
        % function to find the nhood where psi meand and std stabilise
        function find_nhood(obj)
            L = 6;
            nhood_values = [4, 6, 8, 10 , 12, 14, 16];
            colors = mat2cell(jet(length(nhood_values)), ones(1, length(nhood_values)), 3);
            cla(obj.ui.controls.ax2);
            for i = 1:length(nhood_values)
                nhood = nhood_values(i);
                % get the psi for the current nhood
                [iter, parentDir] = obj.getIteration(obj.path);
                particles_path = fullfile(parentDir, sprintf('particle_locations_%s.csv', iter));
                if ~isfile(particles_path)
                    fprintf('Particle locations not found for stack %s\n', obj.path);
                    continue;
                end
                psi = obj.freudLocalOrderParameter(particles_path, obj.stack_info.bd, L, nhood);
                % plot the mean and std of psi
                mean_psi = mean(psi);
                std_psi = std(psi);
                % check if mean and std are stable
                if i > 1 && abs(mean_psi - mean_psi_prev) < 0.01 && abs(std_psi - std_psi_prev) < 0.01
                    fprintf('Stable psi found for nhood = %d\n', nhood);
                    break;
                end
                mean_psi_prev = mean_psi;
                std_psi_prev = std_psi;
                % plot the psi histogram on obj.ui.controls.ax2
                histogram(obj.ui.controls.ax2, psi, 'Normalization', 'pdf', 'FaceColor', colors{i}, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                hold(obj.ui.controls.ax2, 'on');
            end
            obj.stack_info.nhood = nhood;
            obj.save_stack_callback();
            % save the obj.ui.controls.ax2
            exportgraphics(obj.ui.controls.ax2, fullfile(parentDir, sprintf('nhood_vs_psi_%s.png', iter)));
        end
        function find_all_nhood(obj)
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.path = obj.stack_paths{i};
                if contains(obj.path, 'time_control') || contains(obj.path, 'temp')
                    fprintf('Skipping %s\n', obj.path);           
                    continue;
                end
                obj.load_stack_info();
                obj.find_nhood();
            end
        end
        function find_all_psi(obj)
            % iterate over all the stacks
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.path = obj.stack_paths{i};
                if contains(obj.path, 'time_control') || contains(obj.path, 'temp')
                    fprintf('Skipping %s\n', obj.path);           
                    continue;
                end
                obj.load_stack_info();
                [iter, parentDir] = obj.getIteration(obj.path);
                particles_path = fullfile(parentDir, sprintf('particle_locations_%s.csv', iter));
                if ~isfile(particles_path)
                    fprintf('Particle locations not found for stack %s\n', obj.path);
                    continue;
                end
                bead_dia = obj.stack_info.bd;
                nhood = obj.stack_info.nhood;
                l = 6;
                psi = obj.freudLocalOrderParameter(particles_path, bead_dia, l, nhood);
                obj.stack_info.(sprintf("psi_%d",l)) = psi;
                % save the psi to the stack_info
                % assignin('base', 'stack_info', stack_info);
                obj.save_stack_callback();
                fprintf('Psi calculated for stack %s\n', obj.path);
            end
        end
    %%%%%%%%%%%%%%%%%%%%%% GR %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Gr_all_stacks_callback(obj, ~,~)
            % iterate over all the stacks
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.path = obj.stack_paths{i};
                if contains(obj.path, 'time_control') || contains(obj.path, 'temp') || contains(obj.path, 'cont')
                    continue;
                end
                obj.load_images_callback();
                if isfield(obj.stack_info, 'gr')
                    fprintf('Gr already exists for stack %s\n', obj.path);
                    if ~obj.forced
                        fprintf('Skipping stack %s\n', obj.path);
                        continue;
                    end
                end
                obj.get_Gr('mode', 'auto');
            end
        end
        function get_Gr(obj, ~,~)
            radius = 15 * 5;
            start = obj.stack_info.start_index;
            image_path = fullfile(obj.stack_info.img_data.img_files(start).folder, obj.stack_info.img_data.img_files(start).name);
            [iter, parentDir] = obj.getIteration(obj.path);
            save_path = fullfile(parentDir, sprintf('particle_locations_%s.csv', iter));

            % if save_path doesn't exist, get the particle locations from the first image
            if ~isfile(save_path)
                disp('getting particle locations');
                obj.get_particle_locations(image_path, save_path);
            end

            % check if gr exists in stack_info
            if ~isfield(obj.stack_info, 'gr') || obj.forced
                disp('calculating gr');
                % calculate gr for the first image
                bin_width = 3;
                obj.calculate_gr(radius, bin_width);
                % assignin('base', 'gr', gr);
            end
            % plot the gr
            obj.plot_gr(iter, parentDir);
            fprintf('Gr calculated for stack %s\n', obj.path);
        end
        function gr = calculate_gr(obj, r_max, dr)
            % calculate the radial distribution function
            % load the particle locations
            particle_locations = obj.stack_info.particle_locations;
            % bins
            bin_centers = dr:0.1:r_max-dr;
            % calculate the histogram
            disp('calculating radial distribution function');

            % Option 2: For precise window counting (if dr != bin width/2)
            counts = get_counts(bin_centers, dr, particle_locations);
            fprintf('Counts before normalizing: %d\n', counts);
            gr = normalize_count(counts, particle_locations);

            obj.stack_info.gr = gr;
            obj.stack_info.gr_bins = bin_centers;
            % assignin('base', 'stack_info', stack_info);
            obj.save_stack_callback();
            function counts = get_counts(bins_centers, dr, particle_locations)
                counts = zeros(1, numel(bins_centers));
                wait = waitbar(0, 'Calculating Count for each Particle');
                for i = 1:numel(particle_locations.x)            
                    % get the particle location
                    x = particle_locations.x(i);
                    y = particle_locations.y(i);
                    if obj.stack_info.mask(floor(y), floor(x)) == 1
                        continue;
                    end
                    % get the distances of the particle from all other particles
                    distances = sqrt((particle_locations.x - x).^2 + (particle_locations.y - y).^2);
                    for j = 1:numel(bins_centers)
                        bin_center = bins_centers(j);
                        % get the distances in the range [bin_center-dr, bin_center+dr]
                        in_range = (distances >= bin_center-dr) & (distances < bin_center+dr);
                        % normalize by the area of the bin with correction for border
                        in_range = in_range / get_sliceArea(x, y, bin_center, dr);
                        % add the count to the counts array
                        counts(j) = counts(j) + sum(in_range);                
                    end
                    waitbar(i/numel(particle_locations.x));
                end
                close(wait);
            end
            function area = get_sliceArea(x, y, bin_center, dr)
                % depending on the particle location, bin_center and dr, 
                % get the area of the bin considering the border and the mask
                % usually area = pi*(bin_center+dr)^2 - pi*(bin_center-dr)^2;
                % but if the bin is close to the border or the mask, the area will be less
                percent = calculate_unmasked_ring_percentage([y, x], bin_center, dr, obj.stack_info.mask);
                % if percent == 0
                %     fprintf('Percentage for x: %d, y: %d, bin: %d, percent: %d\n',x,y, bin_center, percent);
                % end
                area = (pi*(bin_center+dr)^2 - pi*(bin_center-dr)^2) * percent;
            end
            function percentage = calculate_unmasked_ring_percentage(center, radius, thickness, masked_image, mask_value)
                % Calculate the percentage of a complete circular ring that is both:
                % 1) Inside the image boundaries
                % 2) Not masked
                %
                % Parameters:
                % -----------
                % center : [row, col]
                %     Center coordinates of the circle in [row, column] format
                % radius : double
                %     Radius of the circle
                % thickness : double
                %     Thickness of the ring (dr), the ring extends from r-dr to r+dr
                % masked_image : 2D matrix
                %     Binary mask where mask_value indicates masked areas
                % mask_value : scalar, optional
                %     Value in masked_image that indicates masked areas (default: 0)
                %
                % Returns:
                % --------
                % percentage : double
                %     Percentage of the complete circular ring that is inside the image
                %     and not masked (0-100%)
                
                % Set default mask value if not provided
                if nargin < 5
                    mask_value = 1;
                end
                
                % Get image dimensions
                [height, width] = size(masked_image);
                
                % Extract center coordinates
                center_row = center(1);
                center_col = center(2);
                
                % Calculate inner and outer radius
                inner_radius = max(0, radius - thickness);
                outer_radius = radius + thickness;
                
                % Calculate the theoretical area of the complete circular ring
                % Area = π(R₂² - R₁²)
                % theoretical_ring_area = pi * (outer_radius^2 - inner_radius^2);
                
                % Create a large enough grid to cover the entire ring (ignoring image boundaries)
                row_min = floor(center_row - outer_radius);
                row_max = ceil(center_row + outer_radius);
                col_min = floor(center_col - outer_radius);
                col_max = ceil(center_col + outer_radius);
                
                % Create coordinate matrices for the full theoretical ring
                [cols, rows] = meshgrid(col_min:col_max, row_min:row_max);
                
                % Calculate distances from the center for all points
                dist_from_center = sqrt((rows - center_row).^2 + (cols - center_col).^2);
                
                % Create a mask for the full circular ring
                ring_mask = (dist_from_center >= inner_radius) & (dist_from_center <= outer_radius);
                
                % Create a mask for points that are inside the image boundaries
                valid_rows = (rows >= 1) & (rows <= height);
                valid_cols = (cols >= 1) & (cols <= width);
                inside_image_mask = valid_rows & valid_cols;
                
                % Count pixels in the theoretical complete ring
                total_ring_pixels = sum(ring_mask(:));
                
                % For pixels inside both the ring and image boundaries, check if they're masked
                valid_points = ring_mask & inside_image_mask;
                
                % Initialize a counter for unmasked pixels
                unmasked_count = 0;
                
                % Get linear indices of points inside both ring and image
                [valid_row_indices, valid_col_indices] = find(valid_points);
                
                % Check each valid point against the mask
                for i = 1:length(valid_row_indices)
                    img_row = valid_row_indices(i) + row_min - 1;  % Convert back to image coordinates
                    img_col = valid_col_indices(i) + col_min - 1;
                    
                    % Check if this point is not masked
                    if masked_image(img_row, img_col) ~= mask_value
                        unmasked_count = unmasked_count + 1;
                    end
                end
                
                % Calculate percentage: (unmasked pixels) / (total theoretical ring pixels) * 100
                if total_ring_pixels > 0
                    percentage = (unmasked_count / total_ring_pixels);
                else
                    percentage = 0;
                end
            end
            function norm_count = normalize_count(counts, particle_locs)
                % 1.Divide your total count by N, the number of reference particles you considered
                %  -- probably the total number of particles in your data.
                norm_count = counts/numel(particle_locs.x);
                % 2. Divide by the area of the bin
                %  -- probably the area of the ring between r and r+dr. 
                % NOTE: The area normalization is being done in the get_count function
                % norm_count = norm_count ./ get_slice_area(bin_centers, dr, particle_locs);
                % 3. Divide by the particle number density, i.e number/V
                if ~isfield(obj.stack_info, 'number_density')
                    calculate_number_density();
                end
                norm_count = norm_count ./ obj.stack_info.number_density;
            end
            function calculate_number_density()
                % calculate the number density of the particles
                % get the particle locations
                particle_locations = obj.stack_info.particle_locations;
                % get the box dimensions
                xmin = min(particle_locations.x);
                xmax = max(particle_locations.x);
                ymin = min(particle_locations.y);
                ymax = max(particle_locations.y);
                % get the mask area
                mask_area = sum(obj.stack_info.mask(:));
                % calculate the area of the box
                area = (xmax - xmin) * (ymax - ymin) - mask_area;
                % calculate the number density
                number_density = numel(particle_locations.x) / area;
                obj.stack_info.number_density = number_density;
            end
        end
        function plot_local_number_density(obj)
            % plot the local number density as a colormap only for 1st image
            if ~isfield(obj.stack_info, 'local_number_density')
                obj.calculate_local_number_density();
            end
            local_number_density = obj.stack_info.local_number_density;
            % obj.save_stack_callback();
            % plot the local number density
            cla(obj.ui.controls.ax2);
            imagesc(obj.ui.controls.ax2, local_number_density);
            set(obj.ui.controls.ax2, 'YDir', 'normal'); % Ensure the y-axis is not flipped
            title('Local number density');
            colorbar;
            clim([-1, 1]); % Set colorbar limits to 0 and 1
        end
        function calculate_local_number_density(obj)
            % get the particle locations
            particle_locations = obj.stack_info.particle_locations;
            bd = obj.stack_info.bd;
            % number density at each pixel is the...
            % number of particles in a box of size 20bd x 20bd
            % divided by the area of the box 20bd x 20bd - mask area
            obj.setFrame(1);
            [image_height,image_width] = size(obj.stack_info.img_data.imgs{obj.stack_info.start_index});
            % fprintf('height: %d, width: %d\n', image_height, image_width);
            fprintf('Calculating local number density\n');
            local_number_density = zeros(image_height, image_width);
            wait = waitbar(0, 'Calculating local number density');
            for i = 1:image_height
                for j = 1:image_width
                    % get the box corners
                    xmin = max(1, j - floor(10*bd));
                    xmax = min(image_width, j + floor(10*bd));
                    ymin = max(1, i - floor(10*bd));
                    ymax = min(image_height, i + floor(10*bd));
                    % fprintf('xmin: %d, xmax: %d, ymin: %d, ymax: %d\n', xmin, xmax, ymin, ymax);
                    % get the particles in the box
                    particles_in_box = particle_locations.x >= xmin & particle_locations.x <= xmax & ...
                        particle_locations.y >= ymin & particle_locations.y <= ymax;
                    % calculate the number density
                    mask_area = sum(obj.stack_info.mask(xmin:xmax, ymin:ymax), 'all');
                    area = (xmax - xmin) * (ymax - ymin) - mask_area;
                    local_number_density(i,j) = sum(particles_in_box) / area;
                    waitbar(i/image_height);
                end
            end
            wait.close();
            obj.stack_info.local_number_density = local_number_density;
            assignin('base', 'stack_info', obj.stack_info);
            % obj.save_stack_callback();
        end
        function plot_gr(obj, iter, parentDir)
            % clear axis
            cla(obj.ui.controls.ax2);
            % plot the gr on obj.ui.controls.ax2        
            plot(obj.ui.controls.ax2, obj.stack_info.gr_bins(1:end)/obj.stack_info.bd, obj.stack_info.gr);
            title('Radial distribution function');
            xlabel('r/bd');
            ylabel('g(r)');
            axis(obj.ui.controls.ax2, 'tight');
            % save the plot
            temp = figure(visible='off');
            plot(obj.stack_info.gr_bins(1:end)/obj.stack_info.bd, obj.stack_info.gr);
            title('Radial distribution function');
            xlabel('r/bd');
            ylabel('g(r)');
            saveas(temp, fullfile(parentDir, sprintf('gr_%s.png', iter)));
            close(temp);
        end
        function plot_all_gr(obj, ~,~)
            Ns = [];    
            % check if struct with gr of all stacks exists at F:\shake_table_data\Results
            if exist('F:\shake_table_data\Results\gr_all_stacks.mat', 'file')
                gr_all_stacks = load('F:\shake_table_data\Results\gr_all_stacks.mat');
                while isfield(gr_all_stacks, 'gr_all_stacks')
                    gr_all_stacks = gr_all_stacks.gr_all_stacks;
                end
            else
                gr_all_stacks = struct();
            end
            % clear axis
            cla(obj.ui.controls.ax2);hold on;
            % iterate over all the stacks
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                current_idx = get(obj.ui.controls.stackDropdown, 'Value');
                obj.path = obj.stack_paths{current_idx};
                [iteration, parentDir] = obj.getIteration(obj.path);
                [N, fs] = obj.get_info(obj.path);
                Ns = [Ns, N];
                if contains(obj.path, 'time_control') || contains(obj.path, 'temp') || contains(obj.path, 'cont')
                    fprintf('Skipping %s\n', obj.path);           
                    continue;
                end
                if obj.is_rec_start_late(N,fs,iteration)
                    fprintf('Ignoring %d,%d,%s\n as rec started late',N,fs,iteration);
                    continue;
                end
                % check if gr_all_stacks has gr, gr_bins data for N, fs, iteration
                if isfield(gr_all_stacks, sprintf('N%d', N)) && ...
                        isfield(gr_all_stacks.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                        isfield(gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration)) && ~obj.forced
                    fprintf('gr data found in gr_all_stacks for %s\n', obj.path);
                    gr = gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr');
                    gr_bins = gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr_bins');
                else
                    if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                        obj.load_stack_info();
                        if isfield(obj.stack_info, 'gr')
                            fprintf('gr data found in stack_info for %s\n', obj.path);
                            % add the data to gr_all_stacks
                            gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr') = obj.stack_info.gr;
                            gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr_bins') = obj.stack_info.gr_bins;
                            gr = obj.stack_info.gr;
                            gr_bins = obj.stack_info.gr_bins;
                        else
                            fprintf('gr data not found in stack_info for %s calculating now\n', obj.path);
                            obj.get_Gr('mode', 'auto');
                        end
                    end
                end
                % plot the gr on obj.ui.controls.ax2
                plot(obj.ui.controls.ax2, gr_bins(1:end)/obj.stack_info.bd, gr, "DisplayName", "None", "Color", obj.get_color(N));
            end
            % add legend with different colors for each N
            unique_Ns = unique(Ns);
            legend_handles = zeros(1, numel(unique_Ns));
            legend_entries = cell(1, numel(unique_Ns));
            
            for j = 1:numel(unique_Ns)
                N_val = unique_Ns(j);
                % Create a "dummy" line just for the legend with the right color
                legend_handles(j) = plot(obj.ui.controls.ax2, NaN, NaN, '-', 'Color', obj.get_color(N_val));
                legend_entries{j} = sprintf('N = %d', N_val);
            end
            
            % Create legend using only our dummy lines
            legend(obj.ui.controls.ax2, legend_handles, legend_entries, 'Location', 'best');
            hold off;
            axis(obj.ui.controls.ax2, 'tight');
            title('Radial distribution function');
            xlabel('r/bd');
            ylabel('g(r)');
            % save the gr_all_stacks
            save('F:\shake_table_data\Results\gr_all_stacks.mat', 'gr_all_stacks');
            % save the plot
            exportgraphics(obj.ui.controls.ax2, 'F:\shake_table_data\Results\gr_all_stacks.png');
        end
    %%%%%%%%%%%%%%%%%%%%%% PARTICLE LOCATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function get_particle_locations(obj, image_path, save_path)
            % get the particle locations from the image
            py.track.find_particle_locations(image_path=image_path, diam=int32(5), max_iterations=int32(10), minmass=int32(1), separation=int32(5), save_path=save_path);
            % load the saved csv from save_path
            particle_locations = readtable(save_path);
            obj.stack_info.particle_locations = particle_locations;
            obj.save_stack_callback();
            % assignin('base', 'particle_locations', particle_locations);
        end
        function toggle_pores(obj, ~,~)

        end
        function toggle_particle_locations(obj, ~,~)
            if obj.particle_locations_visible
                % hide the particle locations
                obj.particle_locations_visible = false;
                % get slider index
                slider_idx = round(get(obj.ui.controls.slider, 'Value'));
                obj.setFrame(slider_idx);
                set(gcf, 'WindowScrollWheelFcn', @obj.scrollWheelMoved);
            else
                % show the particle locations
                if ~isfield(obj.stack_info, 'particle_locations')
                    obj.display_warning('No particle locations found');
                    return;
                end
                obj.particle_locations_visible = true;
                obj.setFrame(1);
                % get the particle locations
                particle_locations = obj.stack_info.particle_locations;
                hold(obj.ui.controls.ax1, 'on');
                plot(obj.ui.controls.ax1, particle_locations.x, particle_locations.y, 'b*');
                hold(obj.ui.controls.ax1, 'off');
                set(gcf, 'WindowScrollWheelFcn', {});
            end
        end
    %%%%%%%%%%%%%%%%%%%%%% SHORTENING %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [start_index, end_index] = shorten_stack_callback(obj, ~, ~, start_, end_)
            obj.remove_shortening_callback();
            % check if stack_info has maxDiffIndex
            if ~isfield(obj.stack_info, 'maxDiffIndex')
                maxDiffIndex = matchImages(300);
                obj.display_warning(sprintf("Max difference index: %d", maxDiffIndex));
            end
            obj.setFrame(obj.stack_info.maxDiffIndex);
            if nargin > 2
                start_index = start_;
                end_index = end_;
            else
                obj.display_warning("go to start frame and press n");
                obj.wait_for_keypress("n");
                start_index = obj.stack_info.start_index + round(get(obj.ui.controls.slider, 'Value')) - 1;
                obj.display_warning("go to end frame and press n");
                obj.setFrame(obj.stack_info.end_index-obj.stack_info.start_index + 1);
                obj.wait_for_keypress("n");
                end_index = obj.stack_info.start_index + round(get(obj.ui.controls.slider, 'Value')) - 1;
            end
            obj.shorten_slider(start_index, end_index);
            obj.stack_info.shortened = true;
            obj.toggle_indicator(obj.ui.info.shortenedIndicator, true);
            % save the start and end indices in a obj.stack_info.mat file
            obj.save_stack_callback();
            % GPU accelerated function to compare first n images with previous image
            function maxDiffIndex = matchImages(n)
                differences = zeros(n, 1);
                % f = waitbar(0,'Please wait...','Name','Aligning stack...');
                WaitMessage = parfor_wait(n, 'Waitbar', true);
                img_data = obj.stack_info.img_data;
                roi = [200, 750, 400, 100]; % [x, y, width, height]
                for k = 2:n
                    % obj.display_warning(num2str(k));
                    baseFileName = img_data.img_files(k-1).name;
                    fullFileName = fullfile(img_data.img_files(k-1).folder, baseFileName);
                    template = imread(fullFileName);
                    template = imcrop(mat2gray(template), roi);
                    displaced_img1 = imtranslate(template, - obj.stack_info.displacements(k-1, :));
                    % imshow(displaced_img1, 'Parent', obj.ui.controls.ax1);
                    baseFileName = img_data.img_files(k).name;
                    fullFileName = fullfile(img_data.img_files(k).folder, baseFileName);
                    bwImage = imread(fullFileName);
                    bwImage = imcrop(mat2gray(bwImage), roi);
                    displaced_img2 = imtranslate(bwImage, - obj.stack_info.displacements(k, :));

                    difference = imabsdiff(displaced_img1, displaced_img2);
                    differences(k) = sum(difference(:));
                    WaitMessage.Send;
                end
                WaitMessage.Destroy
                maxDiffIndex = find(differences == max(differences));
                obj.stack_info.maxDiffIndex = maxDiffIndex;
                obj.save_stack_callback();
            end
        end
        function shorten_all_stack_callback(obj, ~,~)
            % iterate over all the stacks and align them
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.load_images_callback();
                if isfield(obj.stack_info, 'shortened')
                    if obj.stack_info.shortened == true
                        obj.logs{end+1} = sprintf("Trial %s already shortened",obj.stack_paths(i));
                        continue;
                    end
                end
                if contains(obj.path, 'time_control')
                    obj.stack_info.shortened = true;
                    obj.save_stack_callback();
                    continue;
                end
                obj.shorten_stack_callback();
            end
        end
        function remove_shortening_callback(obj, ~, ~)
            % remove the displacements file
            obj.stack_info.shortened = false;
            obj.stack_info.start_index = 1;
            obj.stack_info.end_index = numel(obj.stack_info.img_data.img_files);
            obj.shorten_slider(1, numel(obj.stack_info.img_data.img_files));
            obj.toggle_indicator(obj.ui.info.shortenedIndicator, false);
            % get current frame
            frame = get(obj.ui.controls.slider, 'Value');
            obj.setFrame(frame);
            % add star to save button
            obj.ui.controls.saveButton.String = 'Save *';
        end
    %%%%%%%%%%%%%%%%%%%%%% STACK ALIGNING %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function remove_alignment_callback(obj, ~, ~)
            % remove the displacements file
            obj.stack_info.aligned = false;
            obj.stack_info.displacements = [];
            obj.toggle_indicator(obj.ui.info.alignedIndicator, false);
            % get current frame
            frame = get(obj.ui.controls.slider, 'Value');
            obj.setFrame(frame);
            % add star to save button
            obj.ui.controls.saveButton.String = 'Save *';
        end
        function skip_alignment_callback(obj, ~, ~)
            % skip the current stack
            obj.skip_alignment = true;
        end
        function align_stack_callback(obj, varargin)
            obj.skip_alignment = false;
            obj.ui.controls.skipButton.Enable = 'on';
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
            if(obj.stack_info.img_data.num_imgs == 0)
                obj.display_warning("load some images first");
            else
                % Get the template
                frame_num = 1;
                [template,templatePosition] = getTemplate(mode, frame_num);
                % if the template is all zeros, get the next frame
                while all(template(:) == 0)
                    obj.display_warning("Template is all zeros, getting next frame");
                    frame_num = frame_num + 1;
                    obj.shorten_stack_callback(0, 0, frame_num, obj.stack_info.end_index);
                    [template, templatePosition] = getTemplate(mode, frame_num);
                end
                % match the template with each image in the stack
                displacements = matchTemplate(template, templatePosition);
                if ~obj.skip_alignment
                    % subtract the displacements of the first frame from all the displacements
                    displacements(:,1) = displacements(:,1)-displacements(obj.stack_info.start_index,1);
                    displacements(:,2) = displacements(:,2)-displacements(obj.stack_info.start_index,2);
                    obj.stack_info.displacements = displacements;
                    % Display the displacements
                    obj.plot_displacements();
                    % save displacements for later use
                    obj.stack_info.aligned = true;
                    obj.toggle_indicator(obj.ui.info.alignedIndicator, true);
                    obj.save_stack_callback();
                end
            end
            % disable skip button
            obj.ui.controls.skipButton.Enable = 'off';
            fprintf('Aligned stack %s\n', obj.path);
            function displacements = matchTemplate(template, templatePosition)
                displacements = zeros(obj.stack_info.img_data.num_imgs, 2);
                % f = waitbar(0,'Please wait...','Name','Aligning stack...');
                WaitMessage = parfor_wait(obj.stack_info.img_data.num_imgs, 'Waitbar', true);
                img_data = obj.stack_info.img_data;
                parfor k = 1:obj.stack_info.img_data.num_imgs
                    if obj.skip_alignment
                        error('Alignment skipped');
                    end
                    % obj.display_warning(num2str(k));
                    baseFileName = img_data.img_files(k).name;
                    fullFileName = fullfile(img_data.img_files(k).folder, baseFileName);
                    bwImageGPU = gpuArray(imread(fullFileName));
                    templateGPU = gpuArray(template); % Convert the template to a GPU array
                    % Extract the region of interest (ROI) with `n` pixel padding
                    x1 = max(1, round(templatePosition(1)) - obj.searchWindow); % X start
                    y1 = max(1, round(templatePosition(2)) - obj.searchWindow); % Y start
                    x2 = min(size(bwImageGPU, 2), round(templatePosition(1) + templatePosition(3)) + obj.searchWindow);
                    y2 = min(size(bwImageGPU, 1), round(templatePosition(2) + templatePosition(4)) + obj.searchWindow);
        
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
                image_idx = obj.stack_info.start_index + frame_num - 1;
                obj.display_warning(sprintf("Getting template from image %d frame %d", image_idx, frame_num));
                obj.setFrame(frame_num);
                windowSize = 100;
                x_offset = 5;
                y_offset = 10;
                h = drawrectangle('Parent', obj.ui.controls.ax1,'Position',[800-windowSize-x_offset,800-windowSize-y_offset,windowSize,windowSize]);
                if (mode == "manual")
                    % obj.display_warning("select a template by drawing a rectangle");
                    obj.display_warning("Press enter to confirm the template");
                    wait(h);
                end
                % Wait for the user to press the Enter key to confirm the template
                position = round(h.Position);
                template = imcrop(mat2gray(obj.stack_info.img_data.imgs{image_idx}), position);
            end
        end
        function plot_displacements(obj)
            displacements = obj.stack_info.displacements;
            % clear axis
            cla(obj.ui.controls.ax2);
            % plot the displacements on obj.ui.controls.ax2
            plot(obj.ui.controls.ax2, displacements(:,1), 'r');
            hold(obj.ui.controls.ax2, 'on');
            plot(obj.ui.controls.ax2, displacements(:,2), 'b');
            % draw vertical lines at start and end indices
            if obj.stack_info.shortened == true
                % disp(obj.stack_info.start_index);
                plot(obj.ui.controls.ax2, obj.stack_info.start_index, [-5, 5], 'g');
                plot(obj.ui.controls.ax2, obj.stack_info.end_index, [-5, 5], 'g');
            end
            title('Displacements');
            xlabel('Image number');
            ylabel('Displacement');
            legend('x', 'y');
            axis(obj.ui.controls.ax2, 'tight'); 
        end
        function align_all_stacks_callback(obj, ~,~)
            % iterate over all the stacks and align them
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.load_images_callback();
                if obj.stack_info.aligned == true
                    obj.logs{end+1} = sprintf("Trial %s already aligned",obj.stack_paths(i));
                    continue;
                end
                if contains(obj.path, 'time_control')
                    obj.stack_info.aligned = true;
                    obj.save_stack_callback();
                    continue;
                end
                obj.align_stack_callback('mode', 'auto');          
            end
        end
    %%%%%%%%%%%%%%%%%%%%%% MASKING %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function draw_masks_callback(obj, ~,~)
            % draw two right angle triangles on the image
            if isfield(obj.stack_info, 'masked')
                if obj.stack_info.masked == true
                    left_vertices = obj.stack_info.mask_vertices(1:3, :);
                    right_vertices = obj.stack_info.mask_vertices(4:6, :);
                end
            else
                left_vertices = [0 400;0 800;250 800];
                right_vertices = [800 400;800 800;550 800];
            end 

            left_triangle = drawpolygon(obj.ui.controls.ax1, 'Position', left_vertices);
            right_triangle = drawpolygon(obj.ui.controls.ax1, 'Position', right_vertices);
            obj.display_warning("Press n to save the mask");
            obj.wait_for_keypress("n");
            % get the mask from the drawn triangles
            mask = createMask(left_triangle) | createMask(right_triangle);
            % save the mask
            obj.stack_info.mask = mask;
            obj.stack_info.mask_vertices = [left_triangle.Position; right_triangle.Position];
            obj.stack_info.masked = true;
            obj.save_stack_callback();
            % clear the drawn triangles
            delete(left_triangle);
            delete(right_triangle);
            obj.setFrame(get(obj.ui.controls.slider, 'Value'));
            fprintf('Masked stack %s\n', obj.path);
        end
        function draw_all_masks_callback(obj, ~,~)
            % iterate over all the stacks and align them
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                obj.load_images_callback();
                if isfield(obj.stack_info, 'masked')
                    if obj.stack_info.masked == true
                        obj.logs{end+1} = sprintf("Trial %s already masked",obj.stack_paths(i));
                        continue;
                    end
                end
                if contains(obj.path, 'time_control')
                    obj.stack_info.masked = true;
                    obj.save_stack_callback();
                    continue;
                end
                obj.draw_masks_callback();
            end
        end
        function shorten_slider(obj, start_index, end_index)
            % obj.display_warning(sprintf("Shortening stack from %d to %d", start_index, end_index));
            % shorten stack
            obj.stack_info.start_index = start_index;
            obj.stack_info.end_index = end_index;
            set(obj.ui.controls.slider, 'Min', 1, 'Max', end_index - start_index + 1, 'Value', 1, ...
                'SliderStep', [1/(end_index-start_index) , 1/(end_index-start_index)], 'Visible', 'On');
            set(obj.ui.info.toFrame, 'String', num2str(end_index - start_index + 1));
            obj.setFrame(1);
        end
        function change_drive_callback(obj, ~, ~)
            % change the drive letter
            current_drive = 'E:';
            new_drive = 'F:';
            WaitMessage = parfor_wait(length(obj.stack_paths), 'Waitbar', true);
            for i = 1:length(obj.stack_paths)
                set(obj.ui.controls.stackDropdown, 'Value', i);
                current_idx = get(obj.ui.controls.stackDropdown, 'Value');
                obj.path = obj.stack_paths{current_idx};
                obj.update_info(obj.path);
                [iteration, parentDir] = obj.getIteration(obj.path);
                if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                    obj.stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
                    while isfield(obj.stack_info, 'stack_info')
                        obj.stack_info = obj.stack_info.stack_info;
                    end
                    % replace the parentDir with the new drive letter
                    obj.stack_info.parentDir = strrep(obj.stack_info.parentDir, current_drive, new_drive);
                    % replace the drive letter in each stack_info.img_data.img_files.folder
                    for j = 1:numel(obj.stack_info.img_data.img_files)
                        obj.stack_info.img_data.img_files(j).folder = strrep(obj.stack_info.img_data.img_files(j).folder, current_drive, new_drive);
                    end
                    obj.save_stack_callback();
                    % assignin('base', 'new_stack_info', stack_info);
                else
                    continue;
                end
                WaitMessage.Send;
            end
            WaitMessage.Destroy
        end
    %%%%%%%%%%%%%%%%%%%%%% JOIN STACKS %%%%%%%%%%%%%%%%%%%%%%
        % FUNCTION that takes two stack paths and rename the images in the second stack
        function combine_stacks_callback(obj)        
            % show popup to select the two stacks
            [primary_path, secondary_path] = obj.select_stacks();
            if isempty(primary_path) || isempty(secondary_path)
                return;
            end
            % combine the stacks
            obj.combine_stacks(primary_path, secondary_path);
            % update the stack paths
            obj.stack_paths = get_stack_paths();
        end
        
        function [primary_path, secondary_path] = select_stacks(~)
            % show a dialog to select the two stacks
            primary_path = uigetdir('F:\shake_table_data', 'Select the primary stack');
            if primary_path == 0
                primary_path = '';
            end
            secondary_path = uigetdir('F:\shake_table_data', 'Select the secondary stack');
            if secondary_path == 0
                secondary_path = '';
            end
        end
        function combine_stacks(obj, primary_path,secondary_path)
            fprintf('Combining %s and %s\n', primary_path, secondary_path);
            % get the last image name in the non cont path
            cont_img_files = dir(fullfile(secondary_path, '*.TIF'));
            last_img_name = obj.get_last_image_name(primary_path);
            parts = split(last_img_name, '_');
            last_number_string = strrep(parts{end}, '.TIF', '');
            prefix = parts{1};
            % rename the images in the secondary path        
            wait = waitbar(0, 'Renaming and moving images');
            for i = 1:numel(cont_img_files)
                img_name = cont_img_files(i).name;
                new_number_string = obj.incrementPaddedString(last_number_string, i);
                new_name = sprintf('%s_%s.TIF', prefix, new_number_string);
                fprintf('Copying images from %s to %s\n', fullfile(secondary_path, img_name), fullfile(primary_path, new_name));
                copyfile(fullfile(secondary_path, img_name), fullfile(primary_path, new_name));
                waitbar(i/numel(cont_img_files));
            end
            close(wait);
        end
        function last_img_name = get_last_image_name(path)
            % get the last image name in the path
            img_files = dir(fullfile(path, '*.TIF'));
            last_img_name = img_files(end).name;
        end
        function paddedResult = incrementPaddedString(obj, paddedString, increment)
            % incrementPaddedString - Increments a zero-padded string and preserves padding.
            %
            %   paddedResult = incrementPaddedString(paddedString, increment)
            %
            %   Inputs:
            %       paddedString - A string representing a number with leading zeros (e.g., '00097').
            %       increment    - The amount to increment the number by (e.g., 1, 5, -2).
            %
            %   Output:
            %       paddedResult - A string representing the incremented number, with the same
            %                      leading zero padding as the input string.  Returns empty string
            %                      if input is invalid.
            
                % Check if the input string is a valid number
                if ~ischar(paddedString) || any(~ismember(paddedString, '0123456789'))
                    paddedResult = '';
                    warning('Input paddedString is not a valid numeric string.');
                    return;
                end
            
                % Convert the padded string to a number
                number = str2double(paddedString);
            
                % Add the increment
                newNumber = number + increment;
            
                % Determine the number of leading zeros required
                numDigits = length(paddedString);
            
                % Format the new number back into a string with leading zeros
                paddedResult = sprintf(['%0' num2str(numDigits) 'd'], newNumber);
            
            end
    %%%%%%%%%%%%%%%%%%%%%% UTILITY CALLBACKS %%%%%%%%%%%%%%%%%%%%%%
        function load_stack_info(obj, ~,~)
            current_idx = get(obj.ui.controls.stackDropdown, 'Value');
            obj.path = obj.stack_paths{current_idx};
            [iteration, parentDir] = obj.getIteration(obj.path);

            if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                obj.stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
                while isfield(obj.stack_info, 'stack_info')
                    obj.stack_info = obj.stack_info.stack_info;
                end
            end
        end
        function setFrame(obj, k)
            persistent counter;
            % initialize the counter
            if isempty(counter)
                counter = 0;
            end
            % reset imgs after 100 frames
            if counter > 100
                % obj.display_warning("Resetting images");
                counter = 0;
                obj.stack_info.img_data.imgs = cell(1, obj.stack_info.img_data.num_imgs);
            end
            image_idx = obj.stack_info.start_index + k - 1;
            % fprintf('Frame number: %d\n image idx %d', k, image_idx);
            % Get the current slider value
            set(obj.ui.controls.frameNumber, 'String', num2str(k)); % Update the frame number display
            set(obj.ui.controls.slider, 'Value', k);  % update the slider value
            % load the image if it's not already loaded
            if isempty(obj.stack_info.img_data.imgs{image_idx})            
                % fprintf('Loading image %d as its empty', image_idx);
                obj.stack_info.img_data.imgs{image_idx} = mat2gray(imread(fullfile(obj.stack_info.img_data.img_files(image_idx).folder, obj.stack_info.img_data.img_files(image_idx).name)));
                counter = counter + 1;
            end
            % if aligned displacements exist, apply them to the image
            if obj.stack_info.aligned == true
                displaced_img = imtranslate(obj.stack_info.img_data.imgs{image_idx}, - obj.stack_info.displacements(image_idx, :));
                if isfield(obj.stack_info, 'masked') && obj.stack_info.masked == true && isfield(obj.stack_info, 'mask')
                    displaced_img = imoverlay(displaced_img, obj.stack_info.mask, 'r');
                end 
                imshow(displaced_img, 'Parent', obj.ui.controls.ax1);
            else
                imshow(obj.stack_info.img_data.imgs{image_idx}, 'Parent', obj.ui.controls.ax1);
            end
            if isfield(obj.stack_info, 'timestamps')
                timestamp = obj.stack_info.timestamps{image_idx};
                if ~isempty(timestamp)
                    % draw the timestamp on the image
                    % disp(image_idx);
                    text(obj.ui.controls.ax1, 'Units', 'normalized', 'Position', [0.99, 0.03], ...
                        'String', sprintf("%.2f Sec", ...
                        obj.time_2_sec(timestamp)-obj.time_2_sec(obj.stack_info.timestamps{obj.stack_info.start_index})), ...
                        'Color', 'white', 'FontSize', 18, 'HorizontalAlignment', 'right');
                end
            end
        end
        function toggle_indicator(~, indicator, status)
            if status
                set(indicator, 'BackgroundColor', 'green');
            else
                set(indicator, 'BackgroundColor', 'red');
            end
        end
        function zoom_callback(obj, ~, event,x,y)
            % Get the scroll direction
            if event.VerticalScrollCount > 0
                % Zoom in
                zoomFactor = 1.1;
            else
                % Zoom out
                zoomFactor = 1/1.1;
            end
            % Get the current limits
            xLim = xlim(obj.ui.controls.ax1);
            yLim = ylim(obj.ui.controls.ax1);
            % Get the current point
            % x = event.IntersectionPoint(1);
            % y = event.IntersectionPoint(2);
            % Calculate the new limits
            xLim = (xLim - x) * zoomFactor + x;
            yLim = (yLim - y) * zoomFactor + y;
            % Set the new limits
            xlim(obj.ui.controls.ax1, xLim);
            ylim(obj.ui.controls.ax1, yLim);
        end
        function execute_selected_function(obj, ~, ~)
            % Get the selected function from the dropdown
            selected_index = get(obj.ui.controls.functionDropdown, 'Value');
            selected_function = obj.ui.function_list{selected_index};
            
            % Execute the selected function based on its name
            switch selected_function
                case 'Change Drive'
                    obj.change_drive_callback();            
                % You can add more functions as needed
                case 'Plot all Gr'
                    obj.plot_all_gr();
                case 'Plot all TimeStamps'
                    obj.plot_all_timestamps();
                case 'Get Scales'
                    obj.get_scales();
                case 'Plot scales'
                    obj.plot_scales();
                case 'Local number density'
                    obj.plot_local_number_density();
                case 'Average timeStamps'
                    obj.average_timeStamps();
                case 'Combine stacks'
                    obj.combine_stacks_callback();
                case 'Plot durations'
                    obj.plot_stack_durations();
                case 'Set all Empty or Not'
                    obj.set_all_empty_or_not();
                case 'Find all nhood'
                    obj.find_all_nhood();
                case 'Find all Si6'
                    obj.find_all_psi();
                case 'Set all Jam or Not'
                    obj.set_all_jam_or_not();
            end
            
            obj.display_warning(['Executed: ' selected_function]);
        end
        function forced_callback(obj, ~, ~)
            obj.forced = ~obj.forced;
        end
        function next_stack_callback(obj, ~, ~)
            if get(obj.ui.controls.stackDropdown, 'Value') < length(get(obj.ui.controls.stackDropdown, 'String'))
                set(obj.ui.controls.stackDropdown, 'Value', get(obj.ui.controls.stackDropdown, 'Value') + 1);
                obj.load_images_callback();
            else
                obj.display_warning("You're on the last stack");
            end
        end
        function prev_stack_callback(obj, ~, ~)
            if get(obj.ui.controls.stackDropdown, 'Value') < length(get(obj.ui.controls.stackDropdown, 'String'))
                set(obj.ui.controls.stackDropdown, 'Value', get(obj.ui.controls.stackDropdown, 'Value') - 1);
                obj.load_images_callback();
            else
                obj.display_warning("You're on the last stack");
            end
        end
        function goto_callback(obj, ~, ~)
            % go to the specified stack
            N = str2double(get(obj.ui.info.N_info, 'String'));
            f = str2double(get(obj.ui.info.f_info, 'String'));
            i = str2double(get(obj.ui.info.i_info, 'String'));
            obj.path = sprintf('F:\\shake_table_data\\N%d\\%dhz_hopperflow\\60deg\\10cm\\%d', N, f, i);
            % get index of the matching path from stack_paths
            idx = find(contains(obj.stack_paths, obj.path));
            if isempty(idx)
                obj.display_warning("Invalid stack number");
            else
                set(obj.ui.controls.stackDropdown, 'Value', idx);
                obj.load_images_callback();
            end
        end
        function play_pause_callback(obj, ~, ~)
            if obj.is_playing
                stop(obj.ui.info.playTimer);
                set(obj.ui.controls.playPauseButton, 'CData', obj.ui.play_icon);
            else
                start(obj.ui.info.playTimer);
                set(obj.ui.controls.playPauseButton, 'CData', obj.ui.pause_icon);
            end
            obj.is_playing = ~obj.is_playing;
            % fprintf('speed is : %d', speed);
        end
        function play_timer_callback(obj, ~, ~)
            current_value = get(obj.ui.controls.slider, 'Value');
            if current_value < get(obj.ui.controls.slider, 'Max')
                set(obj.ui.controls.slider, 'Value', current_value + obj.speed);
                obj.slider_callback();
            else
                stop(obj.ui.info.playTimer);
                set(obj.ui.controls.playPauseButton, 'CData', obj.ui.play_icon);
                obj.is_playing = false;
            end
        end
        function open_directory_callback(obj, ~, ~)
            % open the directory in windows explorer
            obj.path = obj.stack_paths{get(obj.ui.controls.stackDropdown, 'Value')};
            [~, parentDir] = obj.getIteration(obj.path);
            winopen(parentDir);
        end
        function slider_callback(obj, ~, ~)
            if ~evalin('base', 'exist(''stack_info'', ''var'')')
                obj.display_warning("load some images first");
            else        
                % Get the current slider value
                slider_value = round(get(obj.ui.controls.slider, 'Value'));
                % if isfield(stack_info, 'dead_zone')
                %     if slider_value >= stack_info.dead_zone(1) && slider_value <= stack_info.dead_zone(2)
                %         obj.display_warning("Dead zone, skipping");
                %         obj.setFrame(stack_info.dead_zone(2) + 1);
                %         set(obj.ui.controls.slider, 'Value', stack_info.dead_zone(2) + 1);
                %     end
                % end
                obj.setFrame(slider_value)
            end
        end
        function set_speed_callback(obj, ~,~)
            speed_idx = get(obj.ui.controls.speedDropdown, 'Value');
            obj.speed = obj.speeds{speed_idx};
        end
        function toggle_get_time_ui(obj)
            if contains(obj.path, 'time_control')
                obj.ui.controls.get_time = uicontrol(obj.ui.buttonPanel, 'Style', 'pushbutton', 'String', 'Get time', ...
                    'Units', 'normalized','Position', [0.2 0.1 0.6 0.06], 'Callback', @get_times, 'Enable', 'on');
                % % add forced checkbox beside get_time
                % uicontrol(buttonPanel, 'Style', 'checkbox', 'String', 'Forced', ...
                %     'Units', 'normalized', 'Position', [0.8 0.1 0.2 0.06]);
            else
                if ~isempty(obj.ui.controls.get_time)
                    delete(obj.ui.controls.get_time);
                end
            end
        end
        function update_info(obj, path)
            % cancel any ongoing operation
            [iteration, ~] = obj.getIteration(path);
            [N, fs] = obj.get_info(path);
            set(obj.ui.info.N_info, 'String', N);
            set(obj.ui.info.f_info, 'String', fs);
            set(obj.ui.info.i_info, 'String', iteration);
        end
        function display_warning(obj, msg)
            duration = 3;
            % Add the message to the logs
            obj.logs{end+1} = msg;
            % Create a text UI control to display the warning message
            warning_text = uicontrol('Style', 'text', 'String', msg, ...
                                        'Units', 'normalized','Position', [0.3 0.91 0.6 0.08], 'BackgroundColor', 'red', ...
                                        'ForegroundColor', 'white', 'FontSize', 16, 'Visible', 'on');
        
            % Use a timer to hide the warning message after the specified duration
            t = timer('StartDelay', duration, 'TimerFcn', {@hide_warning, warning_text});
            start(t);
        end
        function hide_warning(obj, ~, warning_text)
            % Set the Visible property of the warning text to 'off'
            set(warning_text, 'Visible', 'off');
        
            % Delete the text UI control
            delete(warning_text);
        
            % Delete the timer
            delete(obj);
        end
        function scrollWheelMoved(obj, ~, event)
            persistent last_scroll_time;  % will retain its value between calls
            now = datetime('now');
            if isempty(last_scroll_time)
                last_scroll_time = now;  % initialize to the current time
            end
            time_between_scrolls = seconds(now - last_scroll_time);  % in seconds
            last_scroll_time = now;  % update for next time
        
            step_size = max(1, round(1 / time_between_scrolls * 0.5));  % larger step size for smaller time_between_scrolls
        
            current_value = get(obj.ui.controls.slider, 'Value');
            if event.VerticalScrollCount > 0  % if scrolling down
                new_value = current_value - step_size; % decrease the value
            else  % if scrolling up
                new_value = current_value + step_size; % increase the value
            end
            new_value = max(min(new_value, get(obj.ui.controls.slider, 'Max')), get(obj.ui.controls.slider, 'Min')); % ensure the new value is within the slider's range
            set(obj.ui.controls.slider, 'Value', new_value);  % update the slider value
            obj.slider_callback();  % call the slider's callback function to update the display
        end
        function show_logs_callback(obj, ~, ~)
            % Create a new figure for the logs overlay
            log_fig = figure('Name', 'Logs', 'NumberTitle', 'Off', 'Position', [1100 100 400 600]);
            % Create a listbox to display the logs
            uicontrol('Style', 'listbox', 'Parent', log_fig, 'Units', 'normalized', ...
                        'Position', [0 0 1 1], 'String', obj.logs, 'FontSize', 12);
        end
    %%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%
        % function to return a unique color from jet colormap for each N
        function color = get_color(~, N)
            switch N
                case 1
                    color = [0, 0, 0];
                case 4
                    color = [1, 0, 0];
                case 12
                    color = [0, 1, 0];
                case 24
                    color = [0, 0, 1];
                case 48
                    color = [1, 0, 1];
                otherwise
                    color = [0, 0, 0];
            end
        end
        function marker = get_marker(~, N)
            markers = ['+', 'p', 'h', '.', 'o', '*', 'd', '^', 'v', '>', 'x', '<', 's'];
            marker = markers(mod(N, numel(markers)) + 1);
            % fprintf('Marker for N = %d is %s\n',N, marker);
        end
        function wait_for_keypress(~, key_to_wait_for)
            keypressed = 0;
            fig = gcf; % Get current figure handle
            
            function myKeyPressFcn(~, event)
                if strcmp(event.Key, key_to_wait_for)
                    keypressed = 1;
                end
            end
            
            set(fig, 'KeyPressFcn', @myKeyPressFcn);
            
            while keypressed == 0
                drawnow;
                pause(0.05);
            end
        end
        function [N, fs] = get_info(obj, path)
            if isa(path, 'char')
                path = string(path);
            end
            % F:\shake_table_data\N12\10hz_hopperflow\60deg\10cm\1
            parts = path.split("\");
            % print parts
            N = sscanf(parts{3}, 'N%d');
            fs = sscanf(parts{4}, '%dhz_hopperflow');
        end
        function [trial_name, parentDir] = getIteration(~, path)
            if isa(path, 'char')
                path = string(path);
            end
            parts = path.split("\");
            trial_name = parts(end);
            parentDir = strjoin(parts(1:end-1), "\");
        end
        function save_stack_callback(obj, ~,~)
            % save the current stack
            if obj.stack_info.img_data.num_imgs == 0
                obj.display_warning("load some images first");
            else
                while isfield(obj.stack_info, 'stack_info')
                    obj.stack_info = obj.stack_info.stack_info;
                end
                temp_stack_info = obj.stack_info;
                % get the current stack path
                obj.path = obj.stack_paths{get(obj.ui.controls.stackDropdown, 'Value')};
                [iteration, parentDir] = obj.getIteration(obj.path);
                % remove the image data from the stack_info
                obj.stack_info.img_data.imgs = cell(1, obj.stack_info.img_data.num_imgs);
                % save the stack
                save(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'obj.stack_info');
                obj.display_warning(sprintf("Stack saved to %s//stack_info_%s.mat", parentDir, iteration));
                obj.ui.controls.saveButton.String = 'Save';
                obj.stack_info = temp_stack_info;
                assignin('base', 'stack_info', obj.stack_info);
            end
        end
        function stack_info = initialize_stack_info(obj, img_data)
            % Image data
            obj.path = obj.stack_paths{get(obj.ui.controls.stackDropdown, 'Value')};
            [iteration, parentDir] = obj.getIteration(obj.path);
            obj.display_warning("Initializing stack info");
            stack_info = struct('start_index', 1, 'end_index', numel(img_data.img_files),...
                'parentDir', parentDir, 'iteration', iteration, ...
                'aligned', false, 'shortened', false, 'masked', false,...
                'displacements', zeros(numel(img_data.img_files), 2), 'img_data', img_data);
        end

    end
end