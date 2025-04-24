classdef Utils < handle
    properties
        app
    end
    methods
        function obj = Utils(app)
            obj.app = app;
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
                    % get a color from the jet colormap
                    colormap = jet(20);
                    color = colormap(N, :);
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
        function [N, fs] = get_info(~, path)
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
        function save_all_stacks(obj)
            for i = 26:length(obj.app.stack_paths)
                % Set the current stack in the dropdown
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                path = obj.app.stack_paths{i};
                if contains(path, 'time_control') || contains(path, 'temp')
                    fprintf('Skipping %s\n', path);           
                    continue;
                end
                % Load the stack info
                obj.load_stack_info();
                
                % Save the stack info
                obj.save_stack_callback();
                fprintf('Saved stack %d/%d\n', i, length(obj.app.stack_paths));
            end
            obj.display_warning("All stacks have been saved.");
        end
        function save_stack_callback(obj, ~,~)
            % save the current stack
            if obj.app.stack_info.img_data.num_imgs == 0
                obj.display_warning("load some images first");
            else
                WaitMessage = parfor_wait(4, 'Waitbar', true);
                while isfield(obj.app.stack_info, 'stack_info')
                    obj.app.stack_info = obj.app.stack_info.stack_info;
                end
                WaitMessage.Send;
                temp_stack_info = obj.app.stack_info;
                % get the current stack path
                obj.app.path = obj.app.stack_paths{get(obj.app.ui.controls.stackDropdown, 'Value')};
                [iteration, parentDir] = obj.getIteration(obj.app.path);
                % remove the image data from the stack_info
                obj.app.stack_info.img_data.imgs = cell(1, obj.app.stack_info.img_data.num_imgs);
                WaitMessage.Send;
                % save the stack
                stack_info = obj.app.stack_info;
                save(sprintf('%s//stack_info_%s.mat', parentDir, iteration), '-struct', 'stack_info', '-v7.3');
                WaitMessage.Send;
                obj.display_warning(sprintf("Stack saved to %s//stack_info_%s.mat", parentDir, iteration));
                obj.app.ui.controls.saveButton.String = 'Save';
                obj.app.stack_info = temp_stack_info;
                fprintf('Saved stack info to %s//stack_info_%s.mat\n', parentDir, iteration);
                WaitMessage.Send;
                % assignin('base', 'stack_info', obj.app.stack_info);
                WaitMessage.Destroy;
            end
        end
        function stack_info = initialize_stack_info(obj, img_data)
            % Image data
            obj.app.path = obj.app.stack_paths{get(obj.app.ui.controls.stackDropdown, 'Value')};
            [iteration, parentDir] = obj.getIteration(obj.app.path);
            obj.display_warning("Initializing stack info");
            stack_info = struct('start_index', 1, 'end_index', numel(img_data.img_files),...
                'parentDir', parentDir, 'iteration', iteration, ...
                'aligned', false, 'shortened', false, 'masked', false,...
                'displacements', zeros(numel(img_data.img_files), 2), 'img_data', img_data);
        end
        %%%%%%%%%%%%%%%%%%%%%% UTILITY CALLBACKS %%%%%%%%%%%%%%%%%%%%%%
        function load_stack_info(obj, ~,~)
            current_idx = get(obj.app.ui.controls.stackDropdown, 'Value');
            obj.app.path = obj.app.stack_paths{current_idx};
            [iteration, parentDir] = obj.getIteration(obj.app.path);

            if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                obj.app.stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
                while isfield(obj.app.stack_info, 'stack_info')
                    obj.app.stack_info = obj.app.stack_info.stack_info;
                end
                fprintf('Loaded stack info from %s//stack_info_%s.mat\n', parentDir, iteration);
            end
            % assignin('base', 'stack_info', obj.app.stack_info);
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
                obj.app.stack_info.img_data.imgs = cell(1, obj.app.stack_info.img_data.num_imgs);
            end
            image_idx = obj.app.stack_info.start_index + k - 1;
            % fprintf('Frame number: %d\n image idx %d', k, image_idx);
            obj.app.current_image_idx = image_idx; % Update the current image index
            % Get the current slider value
            set(obj.app.ui.controls.frameNumber, 'String', num2str(k)); % Update the frame number display
            set(obj.app.ui.controls.slider, 'Value', k);  % update the slider value
            % load the image if it's not already loaded
            if isempty(obj.app.stack_info.img_data.imgs{image_idx})            
                % fprintf('Loading image %d as its empty', image_idx);
                obj.app.stack_info.img_data.imgs{image_idx} = mat2gray(imread(fullfile(obj.app.stack_info.img_data.img_files(image_idx).folder, obj.app.stack_info.img_data.img_files(image_idx).name)));
                counter = counter + 1;
            end
            % if aligned displacements exist, apply them to the image
            if obj.app.stack_info.aligned == true
                displaced_img = imtranslate(obj.app.stack_info.img_data.imgs{image_idx}, - obj.app.stack_info.displacements(image_idx, :));
                if isfield(obj.app.stack_info, 'masked') && obj.app.stack_info.masked == true && isfield(obj.app.stack_info, 'mask')
                    displaced_img = imoverlay(displaced_img, obj.app.stack_info.mask, 'r');
                end 
                imshow(displaced_img, 'Parent', obj.app.ui.controls.ax1);
            else
                imshow(obj.app.stack_info.img_data.imgs{image_idx}, 'Parent', obj.app.ui.controls.ax1);
            end
            if isfield(obj.app.stack_info, 'timestamps')
                timestamp = obj.app.stack_info.timestamps{image_idx};
                if ~isempty(timestamp)
                    % draw the timestamp on the image
                    % disp(image_idx);
                    text(obj.app.ui.controls.ax1, 'Units', 'normalized', 'Position', [0.99, 0.03], ...
                        'String', sprintf("%.2f Sec", ...
                        obj.app.timer.time_2_sec(timestamp)-obj.app.timer.time_2_sec(obj.app.stack_info.timestamps{obj.app.stack_info.start_index})), ...
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
            xLim = xlim(obj.app.ui.controls.ax1);
            yLim = ylim(obj.app.ui.controls.ax1);
            % Get the current point
            % x = event.IntersectionPoint(1);
            % y = event.IntersectionPoint(2);
            % Calculate the new limits
            xLim = (xLim - x) * zoomFactor + x;
            yLim = (yLim - y) * zoomFactor + y;
            % Set the new limits
            xlim(obj.app.ui.controls.ax1, xLim);
            ylim(obj.app.ui.controls.ax1, yLim);
        end
        function execute_selected_function(obj, ~, ~)
            % Get the selected function from the dropdown
            selected_index = get(obj.app.ui.controls.functionDropdown, 'Value');
            selected_function = obj.app.ui.function_list{selected_index};
            
            % Execute the selected function based on its name
            switch selected_function
                case 'Change Drive'
                    obj.app.masking.change_drive_callback();  
                case 'Save all stacks'
                    obj.save_all_stacks();
                % You can add more functions as needed
                case 'Plot all Gr'
                    obj.app.rdf.plot_all_gr();
                case 'Plot all TimeStamps'
                    obj.app.timer.plot_all_timestamps();
                case 'Get Scales'
                    obj.app.trial.get_scales();
                case 'Set Scales'
                    obj.app.trial.set_bd();
                case 'Plot scales'
                    obj.app.trial.plot_scales();
                case 'Local number density'
                    obj.app.rdf.plot_local_number_density();
                case 'Average timeStamps'
                    obj.app.timer.average_timeStamps();
                case 'Combine stacks'
                    obj.app.combine_stacks_callback();
                case 'Plot durations'
                    obj.app.timer.plot_stack_durations();
                case 'Set all Empty or Not'
                    obj.app.trial.set_all_empty_or_not();
                case 'Find all nhood'
                    obj.app.steinhardt.find_all_nhood();
                case 'Find all Si6'
                    obj.app.steinhardt.find_all_psi();
                case 'Set all Jam or Not'
                    obj.app.trial.set_all_jam_or_not();
                case 'Get LBOOP'
                    obj.app.steinhardt.get_LBOOP();
                case 'Get all LBOOP'
                    obj.app.steinhardt.get_all_LBOOP();
                case 'Plot LBOOP for all N'
                    obj.app.steinhardt.plot_LBOOP_per_chain();
                case 'Plot avg LBOOP for all N'
                    obj.app.steinhardt.plot_avg_LBOOP_per_chain();
                case 'Plot BOOP'
                    obj.app.steinhardt.plot_avg_BOOP();
                case 'Superpose all LBOOP'
                    obj.app.steinhardt.superpose_all_LBOOP();
                case 'Detect voids - stack'
                    obj.app.voids.detect_voids_stack();
                case 'Detect voids - all stacks'
                    obj.app.voids.detect_voids_all_stacks();
                case 'Analyze voids - stack'
                    obj.app.voids.analyze_stack_voids();
                case 'Analyze voids - all stacks'
                    obj.app.voids.analyze_all_stacks_voids();
                case 'Consolidate voids data'
                    obj.app.voids.consolidate_voids_data();
            end
            
            obj.display_warning(['Executed: ' selected_function]);
        end
        function forced_callback(obj, ~, ~)
            obj.app.forced = ~obj.app.forced;
        end
        function next_stack_callback(obj, ~, ~)
            if get(obj.app.ui.controls.stackDropdown, 'Value') < length(get(obj.app.ui.controls.stackDropdown, 'String'))
                set(obj.app.ui.controls.stackDropdown, 'Value', get(obj.app.ui.controls.stackDropdown, 'Value') + 1);
                obj.app.load_images_callback();
            else
                obj.display_warning("You're on the last stack");
            end
        end
        function prev_stack_callback(obj, ~, ~)
            if get(obj.app.ui.controls.stackDropdown, 'Value') < length(get(obj.app.ui.controls.stackDropdown, 'String'))
                set(obj.app.ui.controls.stackDropdown, 'Value', get(obj.app.ui.controls.stackDropdown, 'Value') - 1);
                obj.app.load_images_callback();
            else
                obj.display_warning("You're on the last stack");
            end
        end
        function goto_callback(obj, ~, ~)
            % go to the specified stack
            N = str2double(get(obj.app.ui.info.N_info, 'String'));
            f = str2double(get(obj.app.ui.info.f_info, 'String'));
            i = str2double(get(obj.app.ui.info.i_info, 'String'));
            obj.app.path = sprintf('F:\\shake_table_data\\N%d\\%dhz_hopperflow\\60deg\\10cm\\%d', N, f, i);
            % get index of the matching path from stack_paths
            idx = find(contains(obj.app.stack_paths, obj.app.path));
            if isempty(idx)
                obj.display_warning("Invalid stack number");
            else
                set(obj.app.ui.controls.stackDropdown, 'Value', idx);
                obj.app.load_images_callback();
            end
        end
        function play_pause_callback(obj, ~, ~)
            if obj.app.is_playing
                stop(obj.app.ui.info.playTimer);
                set(obj.app.ui.controls.playPauseButton, 'CData', obj.app.ui.play_icon);
            else
                start(obj.app.ui.info.playTimer);
                set(obj.app.ui.controls.playPauseButton, 'CData', obj.app.ui.pause_icon);
            end
            obj.app.is_playing = ~obj.app.is_playing;
            % fprintf('speed is : %d', speed);
        end
        function play_timer_callback(obj, ~, ~)
            current_value = get(obj.app.ui.controls.slider, 'Value');
            if current_value < get(obj.app.ui.controls.slider, 'Max')
                set(obj.app.ui.controls.slider, 'Value', current_value + obj.app.speed);
                obj.slider_callback();
            else
                stop(obj.app.ui.info.playTimer);
                set(obj.app.ui.controls.playPauseButton, 'CData', obj.app.ui.play_icon);
                obj.app.is_playing = false;
            end
        end
        function open_directory_callback(obj, ~, ~)
            % open the directory in windows explorer
            obj.app.path = obj.app.stack_paths{get(obj.app.ui.controls.stackDropdown, 'Value')};
            [~, parentDir] = obj.getIteration(obj.app.path);
            winopen(parentDir);
        end
        function slider_callback(obj, ~, ~)
            if obj.app.stack_info.img_data.num_imgs == 0
                obj.display_warning("load some images first");
            else        
                % Get the current slider value
                slider_value = round(get(obj.app.ui.controls.slider, 'Value'));
                % if isfield(stack_info, 'dead_zone')
                %     if slider_value >= stack_info.dead_zone(1) && slider_value <= stack_info.dead_zone(2)
                %         obj.display_warning("Dead zone, skipping");
                %         obj.setFrame(stack_info.dead_zone(2) + 1);
                %         set(obj.app.ui.controls.slider, 'Value', stack_info.dead_zone(2) + 1);
                %     end
                % end
                obj.setFrame(slider_value)
            end
        end
        function set_speed_callback(obj, ~,~)
            speed_idx = get(obj.app.ui.controls.speedDropdown, 'Value');
            obj.app.speed = obj.app.speeds{speed_idx};
        end
        function toggle_get_time_ui(obj)
            if contains(obj.app.path, 'time_control')
                obj.app.ui.controls.get_time = uicontrol(obj.app.ui.buttonPanel, 'Style', 'pushbutton', 'String', 'Get time', ...
                    'Units', 'normalized','Position', [0.2 0.1 0.6 0.06], 'Callback', @obj.app.timer.get_times, 'Enable', 'on');
                % % add forced checkbox beside get_time
                % uicontrol(buttonPanel, 'Style', 'checkbox', 'String', 'Forced', ...
                %     'Units', 'normalized', 'Position', [0.8 0.1 0.2 0.06]);
            else
                if ~isempty(obj.app.ui.controls.get_time)
                    delete(obj.app.ui.controls.get_time);
                end
            end
        end
        function update_info(obj, path)
            % cancel any ongoing operation
            [iteration, ~] = obj.getIteration(path);
            [N, fs] = obj.get_info(path);
            set(obj.app.ui.info.N_info, 'String', N);
            set(obj.app.ui.info.f_info, 'String', fs);
            set(obj.app.ui.info.i_info, 'String', iteration);
        end
        function display_warning(obj, msg)
            duration = 3;
            % Add the message to the logs
            obj.app.logs{end+1} = msg;
            % Create a text UI control to display the warning message
            warning_text = uicontrol('Style', 'text', 'String', msg, ...
                                        'Units', 'normalized','Position', [0.3 0.91 0.6 0.08], 'BackgroundColor', 'red', ...
                                        'ForegroundColor', 'white', 'FontSize', 16, 'Visible', 'on');

            % Use a timer to hide the warning message after the specified duration
            t = timer('StartDelay', duration, 'TimerFcn', {@obj.hide_warning, warning_text});
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

            current_value = get(obj.app.ui.controls.slider, 'Value');
            if event.VerticalScrollCount > 0  % if scrolling down
                new_value = current_value - step_size; % decrease the value
            else  % if scrolling up
                new_value = current_value + step_size; % increase the value
            end
            new_value = max(min(new_value, get(obj.app.ui.controls.slider, 'Max')), get(obj.app.ui.controls.slider, 'Min')); % ensure the new value is within the slider's range
            set(obj.app.ui.controls.slider, 'Value', new_value);  % update the slider value
            obj.slider_callback();  % call the slider's callback function to update the display
        end
        function show_logs_callback(obj, ~, ~)
            % Create a new figure for the logs overlay
            log_fig = figure('Name', 'Logs', 'NumberTitle', 'Off', 'Position', [1100 100 400 600]);
            % Create a listbox to display the logs
            uicontrol('Style', 'listbox', 'Parent', log_fig, 'Units', 'normalized', ...
                        'Position', [0 0 1 1], 'String', obj.app.logs, 'FontSize', 12);
        end
        % function all_distances = calculate_all_distances(locations)
        %     x_gpu = gpuArray(locations.x);
        %     y_gpu = gpuArray(locations.y);
        %     n_particles = size(locations, 1);
        %     % Pre-allocate distances array in GPU memory
        %     all_distances = zeros(n_particles*(n_particles-1)/2, 1, 'gpuArray');
        %     % calculate all distances if not already calculated
        %     disp('calculating distances');
        %     % Create indices for vectorized distance calculation
        %     idx = 1;
        %     for i = 1:n_particles-1
        %         % Calculate distances between particle i and all particles j>i
        %         dx = x_gpu(i) - x_gpu(i+1:end);
        %         dy = y_gpu(i) - y_gpu(i+1:end);
                
        %         % Calculate Euclidean distances
        %         d = sqrt(dx.^2 + dy.^2);
                
        %         % Store in the pre-allocated array
        %         n_dists = length(d);
        %         all_distances(idx:idx+n_dists-1) = d;
        %         idx = idx + n_dists;
        %     end
        % end
        % function uniform_locations = get_uniform_distribution(~, particle_locations, n)
        %     % # find box corners
        %     xmin = int32(min(particle_locations.x));
        %     xmax = int32(max(particle_locations.x));
        %     ymin = int32(min(particle_locations.y));
        %     ymax = int32(max(particle_locations.y));
        %     x = randi([xmin, xmax], n, 1);
        %     y = randi([ymin, ymax], n, 1);
        %     uniform_locations = table(x, y);
        % end
        % function mark_dead_zone_callback(~,~)
        %     % get the start and end frame of the dead timeline and save it in the stack_info
        %     display_warning("Select the start frame of the dead timeline");
        %     wait_for_keypress("n");
        %     start_frame = round(get(slider, 'Value'));
        %     display_warning("Select the end frame of the dead timeline");
        %     wait_for_keypress("n");
        %     end_frame = round(get(slider, 'Value'));   
        %     stack_info.dead_zone = [start_frame, end_frame];
        %     save_stack_callback();
        % end
        % function lastNonConstIndex = findLastNonConstIndex(arr)
        %     displacements = stack_info.displacements;
        %     displacements(:,1) = displacements(:,1) - displacements(end,1);
        %     displacements(:,2) = displacements(:,2) - displacements(end,2);
        %     % find the last non zero index
        %     lastNonConstIndex = find(displacements(:,2), 1, 'last');
        % end
    end
end