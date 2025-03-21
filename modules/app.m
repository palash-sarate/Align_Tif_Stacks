classdef app 
    properties
        stack_paths
        ui
        utils
        
        stack_info
        Gr_module
        forced
        speed
        speeds
    end
    methods
        % constructor
        function obj = app(stack_info, utils)
            obj.stack_info = stack_info;
            obj.utils = utils;
            obj.Gr_module = Gr_module(stack_info);
            obj.forced = false;
            obj.speed = 1;
            obj.speeds = {1,2,4,8};
            obj.stack_paths = get_stack_paths();
            obj.ui = ui_builder();
        end

        function obj = forced_callback(obj, ~, ~)
            obj.forced = ~obj.forced;
        end
        function obj = set_speed_callback(obj,~,~)
            speed_idx = get(ui.speed_dropdown, 'Value');
            obj.speed = obj.speeds{speed_idx};
        end
        function load_images_callback(~, ~)
            WaitMessage = parfor_wait(4, 'Waitbar', true);
            current_idx = get(ui.stack_dropdown, 'Value');
            path = stack_paths{current_idx};
            % if path has time_control in it load the get_times button
            toggle_get_time_ui();
    
            update_info(path);
            [iteration, parentDir] = getIteration(path);
            set(ui.stack_label, 'String', sprintf('Stack #%d of %d', current_idx, length(stack_paths)));
    
            if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
                while isfield(stack_info, 'stack_info')
                    stack_info = stack_info.stack_info;
                end
                WaitMessage.Send;
                assignin('base', 'stack_info', stack_info);
            else
                img_data.img_files = dir(fullfile(path, '*.tif'));
                stack_info = initialize_stack_info(img_data);
            end
            
            stack_info.img_data.num_imgs = numel(stack_info.img_data.img_files);
            stack_info.img_data.imgs = cell(1, stack_info.img_data.num_imgs);
            WaitMessage.Send;
            % if start and end indices are set, set the shortened indicator to green
            if stack_info.shortened == true
                toggle_indicator(ui.shortened_indicator, true);
            else
                toggle_indicator(ui.shortened_indicator, false);
            end
            WaitMessage.Send;
            % if displacement_n.mat exists, set the aligned indicator to green
            if stack_info.aligned == true && ~isempty(stack_info.displacements)
                toggle_indicator(ui.aligned_indicator, true);
                % plot the displacements on ui.ax2
                plot_displacements();
            else
                toggle_indicator(ui.aligned_indicator, false);
                % clear axis
                cla(ui.ax2);
            end
            shorten_slider(stack_info.start_index, stack_info.end_index);
            WaitMessage.Send;
            fprintf('Loaded stack %s\n', path);
            WaitMessage.Destroy;
        end
        function load_stack_info(~,~)
            current_idx = get(ui.stack_dropdown, 'Value');
            path = stack_paths{current_idx};
            [iteration, parentDir] = getIteration(path);
    
            if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
                while isfield(stack_info, 'stack_info')
                    stack_info = stack_info.stack_info;
                end
            end
        end
        function play_pause_callback(~, ~)
            if is_playing
                stop(ui.play_timer);
                set(ui.play_pause_button, 'CData', ui.play_icon);
            else
                start(ui.play_timer);
                set(ui.play_pause_button, 'CData', ui.pause_icon);
            end
            is_playing = ~is_playing;
            % fprintf('speed is : %d', speed);
        end
        function play_timer_callback(~, ~)
            current_value = get(ui.slider, 'Value');
            if current_value < get(ui.slider, 'Max')
                set(ui.slider, 'Value', current_value + speed);
                slider_callback();
            else
                stop(ui.play_timer);
                set(ui.play_pause_button, 'CData', ui.play_icon);
                is_playing = false;
            end
        end
        function next_stack_callback(~, ~)
            if get(ui.stack_dropdown, 'Value') < length(get(ui.stack_dropdown, 'String'))
                set(ui.stack_dropdown, 'Value', get(ui.stack_dropdown, 'Value') + 1);
                load_images_callback();
            else
                display_warning("You're on the last stack");
            end
        end
        function prev_stack_callback(~, ~)
            if get(ui.stack_dropdown, 'Value') < length(get(ui.stack_dropdown, 'String'))
                set(ui.stack_dropdown, 'Value', get(ui.stack_dropdown, 'Value') - 1);
                load_images_callback();
            else
                display_warning("You're on the last stack");
            end
        end
        function slider_callback(~, ~)
            if ~evalin('base', 'exist(''stack_info'', ''var'')')
                display_warning("load some images first");
            else        
                % Get the current slider value
                slider_value = round(get(ui.slider, 'Value'));
                if isfield(stack_info, 'dead_zone')
                    if slider_value >= stack_info.dead_zone(1) && slider_value <= stack_info.dead_zone(2)
                        display_warning("Dead zone, skipping");
                        setFrame(stack_info.dead_zone(2) + 1);
                        set(ui.slider, 'Value', stack_info.dead_zone(2) + 1);
                    end
                end
                setFrame(slider_value)
            end
        end
        function goto_callback(~, ~)
            % go to the specified stack
            N = str2double(get(ui.N_info, 'String'));
            f = str2double(get(ui.f_info, 'String'));
            i = str2double(get(ui.i_info, 'String'));
            path = sprintf('F:\\shake_table_data\\N%d\\%dhz_hopperflow\\60deg\\10cm\\%d', N, f, i);
            % get index of the matching path from stack_paths
            idx = find(contains(stack_paths, path));
            if isempty(idx)
                display_warning("Invalid stack number");
            else
                set(ui.stack_dropdown, 'Value', idx);
                load_images_callback();
            end
        end
        function [N, fs] = get_info(path)
            if isa(path, 'char')
                path = string(path);
            end
            % F:\shake_table_data\N12\10hz_hopperflow\60deg\10cm\1
            parts = path.split("\");
            % print parts
            N = sscanf(parts{3}, 'N%d');
            fs = sscanf(parts{4}, '%dhz_hopperflow');
        end
        function update_info(path)
            % cancel any ongoing operation
            [iteration, ~] = getIteration(path);
            [N, fs] = get_info(path);
            set(ui.N_info, 'String', N);
            set(ui.f_info, 'String', fs);
            set(ui.i_info, 'String', iteration);
        end
        function execute_selected_function(~, ~)
            % Get the selected function from the dropdown
            selected_index = get(ui.function_dropdown, 'Value');
            selected_function = function_list{selected_index};
            
            % Execute the selected function based on its name
            switch selected_function
                case 'Change Drive'
                    change_drive_callback();            
                % You can add more functions as needed
                case 'Plot all Gr'
                    plot_all_gr();
                case 'Plot all TimeStamps'
                    plot_all_timestamps();
            end
            
            display_warning(['Executed: ' selected_function]);
        end
        function show_logs_callback(~, ~)
            % Create a new figure for the logs overlay
            log_fig = figure('Name', 'Logs', 'NumberTitle', 'Off', 'Position', [1100 100 400 600]);
            % Create a listbox to display the logs
            uicontrol('Style', 'listbox', 'Parent', log_fig, 'Units', 'normalized', ...
                        'Position', [0 0 1 1], 'String', logs, 'FontSize', 12);
        end
        function display_warning(msg)
            duration = 3;
            % Add the message to the logs
            logs{end+1} = msg;
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
        function scrollWheelMoved(~, event)
            persistent last_scroll_time;  % will retain its value between calls
            now = datetime('now');
            if isempty(last_scroll_time)
                last_scroll_time = now;  % initialize to the current time
            end
            time_between_scrolls = seconds(now - last_scroll_time);  % in seconds
            last_scroll_time = now;  % update for next time
        
            step_size = max(1, round(1 / time_between_scrolls * 0.5));  % larger step size for smaller time_between_scrolls
        
            current_value = get(ui.slider, 'Value');
            if event.VerticalScrollCount > 0  % if scrolling down
                new_value = current_value - step_size; % decrease the value
            else  % if scrolling up
                new_value = current_value + step_size; % increase the value
            end
            new_value = max(min(new_value, get(ui.slider, 'Max')), get(ui.slider, 'Min')); % ensure the new value is within the ui.slider's range
            set(ui.slider, 'Value', new_value);  % update the slider value
            slider_callback();  % call the slider's callback function to update the display
        end
        function open_directory_callback(~, ~)
            % open the directory in windows explorer
            path = stack_paths{get(ui.stack_dropdown, 'Value')};
            [~, parentDir] = getIteration(path);
            winopen(parentDir);
        end
        function toggle_get_time_ui()
            if contains(path, 'time_control')
                get_time = uicontrol(ui.buttonPanel, 'Style', 'pushbutton', 'String', 'Get time', ...
                    'Units', 'normalized','Position', [0.2 0.1 0.6 0.06], 'Callback', @get_times, 'Enable', 'on');
                % % add forced checkbox beside get_time
                % uicontrol(ui.buttonPanel, 'Style', 'checkbox', 'String', 'Forced', ...
                %     'Units', 'normalized', 'Position', [0.8 0.1 0.2 0.06]);
            else
                if ~isempty(get_time)
                    delete(get_time);
                end
            end
        end
        
    end 
end