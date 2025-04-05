clc;
close all;
clear all;  % Clear cached classes

imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
if count(py.sys.path,pwd) == 0
    insert(py.sys.path,int32(0),pwd);
end
setenv('TCL_LIBRARY', 'C:/Python312/tcl/tcl8.6');
setenv('TK_LIBRARY', 'C:/Python312/tcl/tk8.6');
% workspace;  % Make sure the workspace panel is showing.
% start parallel pool
% parpool(4);

stack_paths = get_stack_paths();
start_image_viewer(stack_paths);
% TODO: Time stamp
% TODO: Provision to combine stacks
% TODO: Time stamp from OCR
% TODO: plot phase space from the total time taken for the trial
% TODO:  

function stack_paths = get_stack_paths()
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
                path = sprintf("F:\\shake_table_data\\N%d\\%dhz_hopperflow\\%ddeg\\%dcm\\",Ns(n),fs(freq),deg,wd);
                % find folders in the path directory
                subFolders = dir(path);
                subFolders = subFolders([subFolders.isdir]);  % Keep only directories
                subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % Remove '.' and '..' directories

                for k = 1:length(subFolders)
                    stack_paths = [stack_paths; fullfile(path, subFolders(k).name)];
                end
            end
        end
    end

    stack_paths = [stack_paths; "F:\shake_table_data\time_control\\1"];
    stack_paths = [stack_paths; "F:\shake_table_data\time_control\\2"];
    stack_paths = [stack_paths; "F:\shake_table_data\time_control\\3"];
end
% function to open the given list of tif images in a scrollable window
function start_image_viewer(stack_paths)
    close all;
    % Get information about all monitors
    monitorPositions = get(0, 'MonitorPositions'); 
    monitorChoice = 1;
    % Check if there's a second monitor
    if size(monitorPositions, 1) > 1 && monitorChoice == 2
        % Get second monitor position
        secondMonitor = monitorPositions(2, :);
        
        % Calculate figure position on left side of second monitor
        figWidth = min(1200, secondMonitor(3) * 0.9);
        figHeight = min(800, secondMonitor(4) * 0.9);
        left = secondMonitor(1); % Left edge of second monitor
        bottom = secondMonitor(2) + (secondMonitor(4) - figHeight)/2; % Centered vertically
        
        % Create figure with specific position
        f = figure('Name', 'Stack Viewer', 'NumberTitle', 'Off', ...
            'Position', [left, bottom, figWidth, figHeight], ...
            'WindowScrollWheelFcn', @scrollWheelMoved);
    else
        % Fallback if only one monitor
        f = figure('Name', 'Stack Viewer', 'NumberTitle', 'Off', ...
            'WindowState', 'maximized', 'WindowScrollWheelFcn', @scrollWheelMoved);
    end
    % Create panels on the figure
    buttonPanel = uipanel(f, 'Units', 'normalized', 'Position', [0 0 0.2 1]); % [left bottom width height]
    axesPanel = uipanel(f, 'Units', 'normalized', 'Position', [0.2 0 0.8 1]);
    searchWindow = 50;
    folder_ico = imread('./folder.png');
    folder_ico = imresize(folder_ico, [20, 20]);
    next_ico = imread('./next.png');
    next_ico = imresize(next_ico, [20, 20]);
    prev_ico = imread('./prev.png');
    prev_ico = imresize(prev_ico, [20, 20]);
    % Create axes on the axes panel
    ax1 = axes(axesPanel, 'Position', [0.01 0.2 0.49 0.7]);
    ax2 = axes(axesPanel, 'Position', [0.5 0.2 0.49 0.7]);
    % Create button on the panel
    stack_dropdown = uicontrol(buttonPanel, 'Style', 'popupmenu', 'String', stack_paths, ...
        'Units', 'normalized', 'Position', [0.2 0.9 0.6 0.06], 'Callback', @load_images_callback);
    % add open directory button beside the dropdown
    open_dir_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', '', ...
        'Units', 'normalized','Position', [0.8 0.93 0.08 0.03], 'CData', folder_ico,...
        'Callback', @open_directory_callback, 'Enable', 'on');
    drawMasks_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Draw mask', ...
        'Units', 'normalized','Position', [0.2 0.7 0.4 0.06], 'Callback', @draw_masks_callback, 'Enable', 'on');
    drawAllMasks_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Do all', ...
        'Units', 'normalized','Position', [0.6 0.7 0.2 0.06], 'Callback', @draw_all_masks_callback, 'Enable', 'on');
    % deadZone_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Mark dead zone', ...
    %     'Units', 'normalized','Position', [0.2 0.64 0.6 0.06], 'Callback', @mark_dead_zone_callback, 'Enable', 'on');
    Gr_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Get Gr', ...
        'Units', 'normalized','Position', [0.2 0.64 0.4 0.06], 'Callback', @get_Gr, 'Enable', 'on');
    Gr_AllStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Get all', ...
        'Units', 'normalized','Position', [0.6 0.64 0.2 0.06], 'Callback', @Gr_all_stacks_callback, 'Enable', 'on');
    shortenStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Shorten stack', ...
        'Units', 'normalized','Position', [0.2 0.58 0.4 0.06], 'Callback', @shorten_stack_callback, 'Enable', 'on');
    shortenAllStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Do all', ...
        'Units', 'normalized','Position', [0.6 0.58 0.2 0.06], 'Callback', @shorten_all_stack_callback, 'Enable', 'on');
    alignStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Align stack', ...
        'Units', 'normalized','Position', [0.2 0.52 0.4 0.06], 'Callback', @align_stack_callback, 'Enable', 'on');
    alignAllStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Align all', ...
        'Units', 'normalized','Position', [0.6 0.52 0.2 0.06], 'Callback', @align_all_stacks_callback, 'Enable', 'on');

    timestampStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Timestamp stack', ...
    'Units', 'normalized','Position', [0.2 0.46 0.4 0.06], 'Callback', @timestamp_stack, 'Enable', 'on');
    timestampAllStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Do all', ...
        'Units', 'normalized','Position', [0.6 0.46 0.2 0.06], 'Callback', @timestamp_all_stacks, 'Enable', 'on');

    logs_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Logs', ...
        'Units', 'normalized','Position', [0.51 0.32 0.3 0.06], 'Callback', @show_logs_callback, 'Enable', 'on');
    save_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Save', ...
        'Units', 'normalized','Position', [0.2 0.32 0.3 0.06], 'Callback', @save_stack_callback, 'Enable', 'on');
    uicontrol(buttonPanel, 'Style', 'text', 'String', 'From Frame:', ...
        'Units', 'normalized', 'Position', [0.1 0.22 0.4 0.08]);    %[left bottom width height]
    uicontrol(buttonPanel, 'Style', 'text', 'String', 'To Frame:', ...
        'Units', 'normalized', 'Position', [0.5 0.22 0.4 0.08]);
    from_frame = uicontrol(buttonPanel, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.2 0.22 0.2 0.04], 'String', '1');
    to_frame = uicontrol(buttonPanel, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.6 0.22 0.2 0.04], 'String', '1');
    slider = uicontrol(axesPanel, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.1 0.06 0.6 0.08], ...
        'Callback', @slider_callback, 'Visible', 'on');
    frame_number = uicontrol(axesPanel, 'Style', 'text','Units', 'normalized', 'Position', [0.7 0.04 0.1 0.08], 'String', '1','FontSize', 14);
    % Create a "Next stack" button
    nextStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', '', ...
        'Units', 'normalized', 'Position', [0.1 0.93 0.08 0.03], 'CData', next_ico,...
         'Callback', @next_stack_callback, 'Enable', 'on');
    prevStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', '', ...
        'Units', 'normalized', 'Position', [0.03 0.93 0.08 0.03], 'CData', prev_ico,...
            'Callback', @prev_stack_callback, 'Enable', 'on');
    skip_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Skip', ...
        'Units', 'normalized','Position', [0.2 0.15 0.6 0.06], 'Callback', @skip_alignment_callback, 'Enable', 'off');
    stack_label = uicontrol(buttonPanel, 'Style', 'text', 'String', 'Stack #1', ...
        'Units', 'normalized', 'Position', [0.2 0.96 0.6 0.03]);
    % create info UI
    uicontrol(buttonPanel, 'Style', 'text', 'String', 'N:', ...
        'Units', 'normalized', 'Position', [-0.1 0.82 0.4 0.08]);
    uicontrol(buttonPanel, 'Style', 'text', 'String', 'f:', ...
        'Units', 'normalized', 'Position', [0.15 0.82 0.4 0.08]);
    uicontrol(buttonPanel, 'Style', 'text', 'String', 'Iter:', ...
        'Units', 'normalized', 'Position', [0.4 0.82 0.4 0.08]);
    N_info = uicontrol(buttonPanel, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.15 0.87 0.1 0.04], 'String', '4');
    f_info = uicontrol(buttonPanel, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.4 0.87 0.1 0.04], 'String', '4');
    i_info = uicontrol(buttonPanel, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.65 0.87 0.1 0.04], 'String', '1');
    goto_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'load', ...
    'Units', 'normalized','Position', [0.77 0.87 0.15 0.04], 'Callback', @goto_callback);
    function_list = {'Change Drive', 'Plot all Gr', 'Plot all TimeStamps', 'Get Scales','Plot scales'...
        'Local number density','Average timeStamps','Combine stacks','Plot durations','Set all Empty or Not'};
    function_dropdown = uicontrol(buttonPanel, 'Style', 'popupmenu', ...
        'String', function_list, ...
        'Units', 'normalized', 'Position', [0.2 0.04 0.6 0.06]);

    execute_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', '', ...
        'Units', 'normalized', 'Position', [0.8 0.07 0.08 0.03], 'CData', next_ico,...
        'Callback', @execute_selected_function);
        % uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', '', ...
        %     'Units', 'normalized', 'Position', [0.1 0.93 0.08 0.03], 'CData', next_ico,...
        %      'Callback', @next_stack_callback, 'Enable', 'on');
    get_time = [];path = [];
    % Create play/pause button beside the scrollbar
    play_icon = imread('./play.png');
    play_icon = imresize(play_icon, [40, 40]);
    pause_icon = imread('./pause.png');
    pause_icon = imresize(pause_icon, [40, 40]);
    particle_icon = imread('./particles_icon.png');
    particle_icon = imresize(particle_icon, [40, 40]);
    play_pause_button = uicontrol(axesPanel, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.01 0.06 0.08 0.08], 'CData', play_icon, 'Callback', @play_pause_callback);
    forced_checkbox = uicontrol(buttonPanel, 'Style', 'checkbox', 'String', 'Forced', ...
        'Units', 'normalized', 'Position', [0.8 0.7 0.2 0.06], 'Callback', @forced_callback, 'Enable', 'on');
    speeds = {1,2,4,8};
    speed_dropdown = uicontrol(axesPanel, 'Style', 'popupmenu', 'String', speeds, ...
        'Units', 'normalized', 'Position', [0.01 0.001 0.08 0.08], 'Callback', @set_speed_callback);
    particle_button = uicontrol(axesPanel, 'Style', 'pushbutton', 'CData', particle_icon, ...
        'Units', 'normalized','Position', [0 0.95 0.04 0.05], 'Callback', @toggle_particle_locations);
    particle_locations_visible = false;
    is_playing = false;
    play_timer = timer('ExecutionMode', 'fixedRate', 'Period', 0.1, 'TimerFcn', @play_timer_callback);
    speed = 1;
    % Create indicators for aligned and shortened statuses
    aligned_indicator = create_indicator(buttonPanel, [0 0.8 0.3 0.05], "Aligned", @remove_alignment_callback);
    shortened_indicator = create_indicator(buttonPanel, [0.5 0.8 0.3 0.05], "Shortened", @remove_shortening_callback);
    logs = {}; % Initialize logs list
    stack_info = struct();
    skip_alignment = false;
    forced = false;
    % goto_callback();

    function load_images_callback(~, ~)
        WaitMessage = parfor_wait(4, 'Waitbar', true);
        current_idx = get(stack_dropdown, 'Value');
        path = stack_paths{current_idx};
        % if path has time_control in it load the get_times button
        toggle_get_time_ui();

        update_info(path);
        [iteration, parentDir] = getIteration(path);
        set(stack_label, 'String', sprintf('Stack #%d of %d', current_idx, length(stack_paths)));

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
            toggle_indicator(shortened_indicator, true);
        else
            toggle_indicator(shortened_indicator, false);
        end
        WaitMessage.Send;
        % if displacement_n.mat exists, set the aligned indicator to green
        if stack_info.aligned == true && ~isempty(stack_info.displacements)
            toggle_indicator(aligned_indicator, true);
            % plot the displacements on ax2
            plot_displacements();
        else
            toggle_indicator(aligned_indicator, false);
            % clear axis
            cla(ax2);
        end
        shorten_slider(stack_info.start_index, stack_info.end_index);
        WaitMessage.Send;
        fprintf('Loaded stack %s\n', path);
        WaitMessage.Destroy;
    end
%%%%%%%%%%%%%%%%%%%%%%%% SCALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function get_scales(~,~)
        num_particles_to_check = 5;
        for i = 1:length(stack_paths)
            if contains(stack_paths{i}, 'time_control') || contains(stack_paths{i}, 'temp') || contains(stack_paths{i}, 'cont')
                continue;
            end
            set(stack_dropdown, 'Value', i);
            load_images_callback();
            if isfield(stack_info, 'scale')
                fprintf('Scales already exist for stack %s\n', path);
                if ~forced
                    fprintf('Skipping stack %s\n', path);
                    continue;
                end
            end
            get_scale(num_particles_to_check);
            % average the stack_info.radii
            radii_values = cellfun(@(x) x{2}, stack_info.radii);
            avg_radius = mean(radii_values);
            std_radius = std(radii_values);
            stack_info.scale = 2 * avg_radius;
            stack_info.scale_std = std_radius;
            save_stack_callback();
        end
        function get_scale(n)
            % set the viewer to first image
            setFrame(1);
            fprintf('Getting scales for stack %s\n', path);
            % get the particle locations
            particle_locations = stack_info.particle_locations;
            % zoom onto 10 particles selected randomly from the image
            for i = 1:n
                random_index = randi(size(particle_locations, 1));
                % get the particle location
                x = particle_locations.x(random_index);
                y = particle_locations.y(random_index);
                radius = particle_locations.size(random_index) * 3;
                set(gcf, 'WindowScrollWheelFcn', {@zoom_callback,x,y});
                % zoom onto the particle
                xlim(ax1, [x - 50, x + 50]);
                ylim(ax1, [y - 50, y + 50]);
    
                % highlight the particle
                hold(ax1, 'on');
                plot(ax1, x, y, 'r*');
                hold(ax1, 'off');
                % get the scale from the user by asking to draw a circle around the particle
                h = drawcircle(ax1, 'Center', [x, y], 'Radius', radius);
                wait(h);
                % save the circle to the stack_info
                % add the scale to the stack_info
                if ~isfield(stack_info, 'radii')
                    stack_info.radii = {};
                end
                stack_info.radii{end+1} = {h.Center, h.Radius};
            end
            % delete the circle
            delete(h);
            % reset the zoom
            % zoom(ax1,'reset');
            % enable the scroll callback
            set(gcf, 'WindowScrollWheelFcn', @scrollWheelMoved);
            % assignin('base', 'stack_info', stack_info);
        end
    end
    function plot_scales()
        scales = zeros(1, length(stack_paths));
        sizes = zeros(1, length(stack_paths));
        scales_std = zeros(1, length(stack_paths));
        sizes_std = zeros(1, length(stack_paths));
        % check if csv already exists
        if isfile('F:\shake_table_data\Results\scales.csv')
            data = readmatrix('F:\shake_table_data\Results\scales.csv');
            scales = data(:,1);
            scales_std = data(:,2);
            sizes = data(:,3);
            sizes_std = data(:,4);
        else
            for i = 1:length(stack_paths)
                if contains(stack_paths{i}, 'time_control') || contains(stack_paths{i}, 'temp') || contains(stack_paths{i}, 'cont')
                    continue;
                end
                set(stack_dropdown, 'Value', i);
                load_images_callback();
                if isfield(stack_info, 'scale')                
                    scales(i) = stack_info.scale;
                    scales_std(i) = stack_info.scale_std;
                end
                sizes(i) = mean(stack_info.particle_locations.size);
                sizes_std(i) = std(stack_info.particle_locations.size);
            end
        end
        cla(ax2);

        errorbar(ax2, scales,scales_std, 'b', 'LineWidth', 2);
        hold(ax2, 'on');
        errorbar(ax2, sizes * 3 * 2,sizes_std, 'g', 'LineWidth', 2);
        hold(ax2, 'off');
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
    function validity = isValidTimeStamp(timestamp)
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
    function timestamp_stack(~,~)
        if ~isfield(stack_info, 'timestamps')
            stack_info.timestamps = cell(1, stack_info.img_data.num_imgs);
        end
        binary_timestamp_available = true;
        WaitMessage = parfor_wait(stack_info.img_data.num_imgs, 'Waitbar', true);
        for i = 1:stack_info.img_data.num_imgs
            if binary_timestamp_available
                timestamp = get_binary_timestamp(i);
                if ~isempty(timestamp)
                    stack_info.timestamps{i} = timestamp;
                else
                    binary_timestamp_available = false;
                    fprintf('Binary timestamp not available for stack %s\n', path);
                end 
            end
            if ~binary_timestamp_available
                % if one time stamp is invalid, then binary stamps don't exist
                % predict timestamp from the averaged time graph
                timestamp = predict_timeStamp(i);
                stack_info.timestamps{i} = timestamp;
            end
            WaitMessage.Send;
        end
        save_stack_callback();
        WaitMessage.Destroy;
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
    end
    function timestamp_all_stacks = average_timeStamps()
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
    function timestamp_all_stacks(~,~)
        WaitMessage = parfor_wait(length(stack_paths), 'Waitbar', true);
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
            path = stack_paths{i};
            if contains(path, 'temp')
                fprintf('Skipping %s\n', path);           
                continue;
            end
            load_images_callback();
            timestamp_stack();
            WaitMessage.Send;
        end
        WaitMessage.Destroy;
    end
    % function to start getting the times for every nth image and also has the option to stop the current operation
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
        cla(ax2);hold on;
        WaitMessage = parfor_wait(length(stack_paths), 'Waitbar', true);
        % iterate over all the stacks
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
            current_idx = get(stack_dropdown, 'Value');
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
                plot(ax2, cellfun(@(t) t.time_us, timestamps) - timestamps{1}.time_us, 'Color', get_color(N), 'DisplayName', 'None');
            else
                fprintf('Timestamps not found in timestamp_all_stacks for %s\n now checking stack_info\n', path);
                if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                    load_stack_info();
                    if isfield(stack_info, 'timestamps') && isValidTimeStamp(stack_info.timestamps{1})
                        timestamps = stack_info.timestamps;
                        timestamp_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('timestamps') = timestamps;
                        % fprintf('Timestamps found in stack_info for %s\n plotting in colour %s', path,mat2str(get_color(N)));
                        plot(ax2, cellfun(@(t) t.time_us, timestamps) - timestamps{1}.time_us, 'Color', get_color(N), 'DisplayName', 'None');
                    end
                else
                    fprintf('Timestamps not found in stack_info for %s\n', path);
                end
            end
            WaitMessage.Send;
        end
        hold off;
        title('Timestamps');
        xlabel('Frame number');
        % average the timestamps
        timestamp_all_stacks = average_timeStamps();
        % save the gr_all_stacks
        save('F:\shake_table_data\Results\timestamp_all_stacks.mat', 'timestamp_all_stacks');
        % save the plot
        saveas(ax2, 'F:\shake_table_data\Results\timestamps.png');
        WaitMessage.Destroy;
    end
    function timestamp = predict_timeStamp(image_idx)
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
    function get_stack_durations()
        % loop over all stacks and get duration by 
        % subtracting the start_index time stamp from the end_index timestamp
        stack_durations = struct();
        WaitMessage = parfor_wait(length(stack_paths), 'Waitbar', true);
        
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
            path = stack_paths{i};
            [N, fs] = get_info(path);
            [iteration, ~] = getIteration(path);
            if contains(path, 'temp') || contains(path, 'time_control')
                fprintf('Skipping %s\n', path);           
                continue;
            end
            load_stack_info();
            stack_info.duration = time_2_sec(stack_info.timestamps{stack_info.end_index}) - ...
                time_2_sec(stack_info.timestamps{stack_info.start_index});
            stack_durations.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)) = stack_info.duration;
            fprintf('Duration for stack %s is %d',path,stack_info.duration);
            % save_stack_callback();
            % assignin('base', 'stack_info', stack_info);
            WaitMessage.Send;
        end
        save('F:\shake_table_data\Results\stack_durations.mat', 'stack_durations');
        WaitMessage.Destroy;
    end
    function plot_stack_durations()
        % check if struct with gr of all stacks exists at F:\shake_table_data\Results
        if exist('F:\shake_table_data\Results\stack_durations.mat', 'file')
            stack_durations = load('F:\shake_table_data\Results\stack_durations.mat');
            while isfield(stack_durations, 'stack_durations')
                stack_durations = stack_durations.stack_durations;
            end
        else
            fprintf('No stack_durations found. Please run get_stack_durations\n');
            get_stack_durations();
            plot_stack_durations();
            return;
        end
        % clear axis
        cla(ax2);hold on;
        % iterate over all the stacks
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
            current_idx = get(stack_dropdown, 'Value');
            path = stack_paths{current_idx};
            [iteration, ~] = getIteration(path);
            [N, fs] = get_info(path);
            if contains(path, 'time_control') || contains(path, 'temp')
                fprintf('Skipping %s\n', path);           
                continue;
            end
            if isfield(stack_durations, sprintf('N%d', N)) && ...
                    isfield(stack_durations.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                    isfield(stack_durations.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration))
                fprintf('Duration found in stack_durations for %s\n', path);
                duration = stack_durations.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration));
                plot(ax2, fs, duration, get_marker(N), 'Color', get_color(N),'HandleVisibility','off');
            else
                fprintf('Duration not found in stack_durations for %s\n', path);
            end
        end
        unique_Ns = [4,12,24,48];
        legend_handles = zeros(1, numel(unique_Ns));
        legend_entries = cell(1, numel(unique_Ns));
        
        for j = 1:numel(unique_Ns)
            N_val = unique_Ns(j);
            % Create a "dummy" line just for the legend with the right color
            legend_handles(j) = plot(ax2, NaN, NaN, get_marker(N_val), 'Color', get_color(N_val));
            legend_entries{j} = sprintf('N = %d', N_val);
        end
        
        % Create legend using only our dummy lines
        legend(ax2, legend_handles, legend_entries, 'Location', 'best');
        hold off;
        xlim([0, 22]);
        ylim([0, 1200]);
        title('Stack durations');
        xlabel('Frequency (Hz)');
        ylabel('Duration (s)');
        exportgraphics(ax2, 'F:\shake_table_data\Results\stack_durations.png');
    end
    function secs = time_2_sec(timestamp)
        secs = timestamp.h * 3600 + timestamp.min * 60 + timestamp.s + timestamp.us / 1e6;
    end
%%%%%%%%%%%%%%%%%%%%%% TRIAL CHARACTERIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
 function is_late = is_rec_start_late(n,fs,iter)
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
%  function to loop over all stack and set if the stack emptied fully or not
% should load the stack and go to the last frame and provide option to user to 
% state whether the hopper has emptied or not
function set_all_empty_or_not()
    % loop over all stacks and set the emptyornot field in stack_info
    % WaitMessage = parfor_wait(length(stack_paths), 'Waitbar', true);
    for i = 1:length(stack_paths)
        set(stack_dropdown, 'Value', i);
        path = stack_paths{i};
        % [N, fs] = get_info(path);
        % [iteration, ~] = getIteration(path);
        if contains(path, 'time_control') || contains(path, 'temp')
            fprintf('Skipping %s\n', path);           
            continue;
        end
        load_images_callback();
        setFrame(stack_info.end_index - stack_info.start_index + 1);
        % show dialog box to user to set the emptyornot field
        options = {'Empty', 'Not Empty'};
        answer = questdlg('Is the hopper empty?', 'Empty or Not', options{1}, options{2}, options{1});
        if strcmp(answer, options{1})
            stack_info.empty = true;
        else
            stack_info.empty = false;
        end
        save_stack_callback();
        % WaitMessage.Send;
    end
    % WaitMessage.Destroy;
end

%%%%%%%%%%%%%%%%%%%%%% GR %%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Gr_all_stacks_callback(~,~)
        % iterate over all the stacks
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
            path = stack_paths{i};
            if contains(path, 'time_control') || contains(path, 'temp') || contains(path, 'cont')
                continue;
            end
            load_images_callback();
            if isfield(stack_info, 'gr')
                fprintf('Gr already exists for stack %s\n', path);
                if ~forced
                    fprintf('Skipping stack %s\n', path);
                    continue;
                end
            end
            get_Gr('mode', 'auto');
        end
    end
    function get_Gr(~,~)
        radius = 15 * 5;
        start = stack_info.start_index;
        image_path = fullfile(stack_info.img_data.img_files(start).folder, stack_info.img_data.img_files(start).name);
        [iter, parentDir] = getIteration(path);
        save_path = fullfile(parentDir, sprintf('particle_locations_%s.csv', iter));

        % if save_path doesn't exist, get the particle locations from the first image
        if ~isfile(save_path)
            disp('getting particle locations');
            get_particle_locations(image_path, save_path);
        end

        % check if gr exists in stack_info
        if ~isfield(stack_info, 'gr') || forced
            disp('calculating gr');
            % calculate gr for the first image
            bin_width = 3;
            calculate_gr(radius, bin_width);
            % assignin('base', 'gr', gr);
        end
        % plot the gr
        plot_gr(iter, parentDir);
        fprintf('Gr calculated for stack %s\n', path);
    end
    function gr = calculate_gr(r_max, dr)
        % calculate the radial distribution function
        % load the particle locations
        particle_locations = stack_info.particle_locations;
        % bins
        bin_centers = dr:0.1:r_max-dr;
        % calculate the histogram
        disp('calculating radial distribution function');

        % Option 2: For precise window counting (if dr != bin width/2)
        counts = get_counts(bin_centers, dr, particle_locations);
        fprintf('Counts before normalizing: %d\n', counts);
        gr = normalize_count(counts, particle_locations);

        stack_info.gr = gr;
        stack_info.gr_bins = bin_centers;
        % assignin('base', 'stack_info', stack_info);
        save_stack_callback();
        function counts = get_counts(bins_centers, dr, particle_locations)
            counts = zeros(1, numel(bins_centers));
            wait = waitbar(0, 'Calculating Count for each Particle');
            for i = 1:numel(particle_locations.x)            
                % get the particle location
                x = particle_locations.x(i);
                y = particle_locations.y(i);
                if stack_info.mask(floor(y), floor(x)) == 1
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
            percent = calculate_unmasked_ring_percentage([y, x], bin_center, dr, stack_info.mask);
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
            if ~isfield(stack_info, 'number_density')
                calculate_number_density();
            end
            norm_count = norm_count ./ stack_info.number_density;
        end
        function calculate_number_density()
            % calculate the number density of the particles
            % get the particle locations
            particle_locations = stack_info.particle_locations;
            % get the box dimensions
            xmin = min(particle_locations.x);
            xmax = max(particle_locations.x);
            ymin = min(particle_locations.y);
            ymax = max(particle_locations.y);
            % get the mask area
            mask_area = sum(stack_info.mask(:));
            % calculate the area of the box
            area = (xmax - xmin) * (ymax - ymin) - mask_area;
            % calculate the number density
            number_density = numel(particle_locations.x) / area;
            stack_info.number_density = number_density;
        end
    end
    function plot_local_number_density()
        % plot the local number density as a colormap only for 1st image
        if ~isfield(stack_info, 'local_number_density')
            calculate_local_number_density();
        end
        local_number_density = stack_info.local_number_density;
        % save_stack_callback();
        % plot the local number density
        cla(ax2);
        imagesc(ax2, local_number_density);
        set(ax2, 'YDir', 'normal'); % Ensure the y-axis is not flipped
        title('Local number density');
        colorbar;
        caxis([-1, 1]); % Set colorbar limits to 0 and 1
    end
    function calculate_local_number_density()
        % get the particle locations
        particle_locations = stack_info.particle_locations;
        bd = stack_info.bd;
        % number density at each pixel is the...
        % number of particles in a box of size 20bd x 20bd
        % divided by the area of the box 20bd x 20bd - mask area
        setFrame(1);
        [image_height,image_width] = size(stack_info.img_data.imgs{stack_info.start_index});
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
                mask_area = sum(stack_info.mask(xmin:xmax, ymin:ymax), 'all');
                area = (xmax - xmin) * (ymax - ymin) - mask_area;
                local_number_density(i,j) = sum(particles_in_box) / area;
                waitbar(i/image_height);
            end
        end
        wait.close();
        stack_info.local_number_density = local_number_density;
        assignin('base', 'stack_info', stack_info);
        % save_stack_callback();
    end
    function plot_gr(iter, parentDir)
        % clear axis
        cla(ax2);
        % plot the gr on ax2        
        plot(ax2, stack_info.gr_bins(1:end)/stack_info.bd, stack_info.gr);
        title('Radial distribution function');
        xlabel('r/bd');
        ylabel('g(r)');
        axis(ax2, 'tight');
        % save the plot
        temp = figure(visible='off');
        plot(stack_info.gr_bins(1:end)/stack_info.bd, stack_info.gr);
        title('Radial distribution function');
        xlabel('r/bd');
        ylabel('g(r)');
        saveas(temp, fullfile(parentDir, sprintf('gr_%s.png', iter)));
        close(temp);
    end
    function plot_all_gr(~,~)
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
        cla(ax2);hold on;
        % iterate over all the stacks
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
            current_idx = get(stack_dropdown, 'Value');
            path = stack_paths{current_idx};
            [iteration, parentDir] = getIteration(path);
            [N, fs] = get_info(path);
            Ns = [Ns, N];
            if contains(path, 'time_control') || contains(path, 'temp') || contains(path, 'cont')
                fprintf('Skipping %s\n', path);           
                continue;
            end
            if is_rec_start_late(N,fs,iteration)
                fprintf('Ignoring %d,%d,%s\n as rec started late',N,fs,iteration);
                continue;
            end
            % check if gr_all_stacks has gr, gr_bins data for N, fs, iteration
            if isfield(gr_all_stacks, sprintf('N%d', N)) && ...
                    isfield(gr_all_stacks.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                    isfield(gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration)) && ~forced
                fprintf('gr data found in gr_all_stacks for %s\n', path);
                gr = gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr');
                gr_bins = gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr_bins');
            else
                if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                    load_stack_info();
                    if isfield(stack_info, 'gr')
                        fprintf('gr data found in stack_info for %s\n', path);
                        % add the data to gr_all_stacks
                        gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr') = stack_info.gr;
                        gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr_bins') = stack_info.gr_bins;
                        gr = stack_info.gr;
                        gr_bins = stack_info.gr_bins;
                    else
                        fprintf('gr data not found in stack_info for %s calculating now\n', path);
                        get_Gr('mode', 'auto');
                    end
                end
            end
            % plot the gr on ax2
            plot(ax2, gr_bins(1:end)/stack_info.bd, gr, "DisplayName", "None", "Color", get_color(N));
        end
        % add legend with different colors for each N
        unique_Ns = unique(Ns);
        legend_handles = zeros(1, numel(unique_Ns));
        legend_entries = cell(1, numel(unique_Ns));
        
        for j = 1:numel(unique_Ns)
            N_val = unique_Ns(j);
            % Create a "dummy" line just for the legend with the right color
            legend_handles(j) = plot(ax2, NaN, NaN, '-', 'Color', get_color(N_val));
            legend_entries{j} = sprintf('N = %d', N_val);
        end
        
        % Create legend using only our dummy lines
        legend(ax2, legend_handles, legend_entries, 'Location', 'best');
        hold off;
        axis(ax2, 'tight');
        title('Radial distribution function');
        xlabel('r/bd');
        ylabel('g(r)');
        % save the gr_all_stacks
        save('F:\shake_table_data\Results\gr_all_stacks.mat', 'gr_all_stacks');
        % save the plot
        exportgraphics(ax2, 'F:\shake_table_data\Results\gr_all_stacks.png');
    end
%%%%%%%%%%%%%%%%%%%%%% PARTICLE LOCATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
    function get_particle_locations(image_path, save_path)
        % get the particle locations from the image
        py.track.find_particle_locations(image_path=image_path, diam=int32(5), max_iterations=int32(10), minmass=int32(1), separation=int32(5), save_path=save_path);
        % load the saved csv from save_path
        particle_locations = readtable(save_path);
        stack_info.particle_locations = particle_locations;
        save_stack_callback();
        % assignin('base', 'particle_locations', particle_locations);
    end
    function toggle_particle_locations(~,~)
        if particle_locations_visible
            % hide the particle locations
            particle_locations_visible = false;
            % get slider index
            slider_idx = round(get(slider, 'Value'));
            setFrame(slider_idx);
            set(gcf, 'WindowScrollWheelFcn', @scrollWheelMoved);
        else
            % show the particle locations
            if ~isfield(stack_info, 'particle_locations')
                display_warning('No particle locations found');
                return;
            end
            particle_locations_visible = true;
            setFrame(1);
            % get the particle locations
            particle_locations = stack_info.particle_locations;
            hold(ax1, 'on');
            plot(ax1, particle_locations.x, particle_locations.y, 'b*');
            hold(ax1, 'off');
            set(gcf, 'WindowScrollWheelFcn', {});
        end
    end
%%%%%%%%%%%%%%%%%%%%%% SHORTENING %%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            start_index = stack_info.start_index + round(get(slider, 'Value')) - 1;
            display_warning("go to end frame and press n");
            setFrame(stack_info.end_index-stack_info.start_index + 1);
            wait_for_keypress("n");
            end_index = stack_info.start_index + round(get(slider, 'Value')) - 1;
        end
        shorten_slider(start_index, end_index);
        stack_info.shortened = true;
        toggle_indicator(shortened_indicator, true);
        % save the start and end indices in a stack_info.mat file
        save_stack_callback();
        % GPU accelerated function to compare first n images with previous image
        function maxDiffIndex = matchImages(n)
            differences = zeros(n, 1);
            % f = waitbar(0,'Please wait...','Name','Aligning stack...');
            WaitMessage = parfor_wait(n, 'Waitbar', true);
            img_data = stack_info.img_data;
            roi = [200, 750, 400, 100]; % [x, y, width, height]
            for k = 2:n
                % display_warning(num2str(k));
                baseFileName = img_data.img_files(k-1).name;
                fullFileName = fullfile(img_data.img_files(k-1).folder, baseFileName);
                template = imread(fullFileName);
                template = imcrop(mat2gray(template), roi);
                displaced_img1 = imtranslate(template, - stack_info.displacements(k-1, :));
                % imshow(displaced_img1, 'Parent', ax1);
                baseFileName = img_data.img_files(k).name;
                fullFileName = fullfile(img_data.img_files(k).folder, baseFileName);
                bwImage = imread(fullFileName);
                bwImage = imcrop(mat2gray(bwImage), roi);
                displaced_img2 = imtranslate(bwImage, - stack_info.displacements(k, :));

                difference = imabsdiff(displaced_img1, displaced_img2);
                differences(k) = sum(difference(:));
                WaitMessage.Send;
            end
            WaitMessage.Destroy
            maxDiffIndex = find(differences == max(differences));
            stack_info.maxDiffIndex = maxDiffIndex;
            save_stack_callback();
        end
    end
    function shorten_all_stack_callback(~,~)
        % iterate over all the stacks and align them
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
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
    function remove_shortening_callback(~, ~)
        % remove the displacements file
        stack_info.shortened = false;
        stack_info.start_index = 1;
        stack_info.end_index = numel(stack_info.img_data.img_files);
        shorten_slider(1, numel(stack_info.img_data.img_files));
        toggle_indicator(shortened_indicator, false);
        % get current frame
        frame = get(slider, 'Value');
        setFrame(frame);
        % add star to save button
        save_button.String = 'Save *';
    end
%%%%%%%%%%%%%%%%%%%%%% PARTICLE ALIGNING %%%%%%%%%%%%%%%%%%%%%%%%%%%
    function remove_alignment_callback(~, ~)
        % remove the displacements file
        stack_info.aligned = false;
        stack_info.displacements = [];
        toggle_indicator(aligned_indicator, false);
        % get current frame
        frame = get(slider, 'Value');
        setFrame(frame);
        % add star to save button
        save_button.String = 'Save *';
    end
    function skip_alignment_callback(~, ~)
        % skip the current stack
        skip_alignment = true;
    end
    function align_stack_callback(varargin)
        skip_alignment = false;
        skip_button.Enable = 'on';
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
        if(stack_info.img_data.num_imgs == 0)
            display_warning("load some images first");
        else
            % Get the template
            frame_num = 1;
            [template,templatePosition] = getTemplate(mode, frame_num);
            % if the template is all zeros, get the next frame
            while all(template(:) == 0)
                display_warning("Template is all zeros, getting next frame");
                frame_num = frame_num + 1;
                shorten_stack_callback(0, 0, frame_num, stack_info.end_index);
                [template, templatePosition] = getTemplate(mode, frame_num);
            end
            % match the template with each image in the stack
            displacements = matchTemplate(template, templatePosition);
            if ~skip_alignment
                % subtract the displacements of the first frame from all the displacements
                displacements(:,1) = displacements(:,1)-displacements(stack_info.start_index,1);
                displacements(:,2) = displacements(:,2)-displacements(stack_info.start_index,2);
                stack_info.displacements = displacements;
                % Display the displacements
                plot_displacements();
                % save displacements for later use
                stack_info.aligned = true;
                toggle_indicator(aligned_indicator, true);
                save_stack_callback();
            end
        end
        % disable skip button
        skip_button.Enable = 'off';
        fprintf('Aligned stack %s\n', path);
        function displacements = matchTemplate(template, templatePosition)
            displacements = zeros(stack_info.img_data.num_imgs, 2);
            % f = waitbar(0,'Please wait...','Name','Aligning stack...');
            WaitMessage = parfor_wait(stack_info.img_data.num_imgs, 'Waitbar', true);
            img_data = stack_info.img_data;
            parfor k = 1:stack_info.img_data.num_imgs
                if skip_alignment
                    error('Alignment skipped');
                end
                % display_warning(num2str(k));
                baseFileName = img_data.img_files(k).name;
                fullFileName = fullfile(img_data.img_files(k).folder, baseFileName);
                bwImageGPU = gpuArray(imread(fullFileName));
                templateGPU = gpuArray(template); % Convert the template to a GPU array
                % Extract the region of interest (ROI) with `n` pixel padding
                x1 = max(1, round(templatePosition(1)) - searchWindow); % X start
                y1 = max(1, round(templatePosition(2)) - searchWindow); % Y start
                x2 = min(size(bwImageGPU, 2), round(templatePosition(1) + templatePosition(3)) + searchWindow);
                y2 = min(size(bwImageGPU, 1), round(templatePosition(2) + templatePosition(4)) + searchWindow);
    
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
            image_idx = stack_info.start_index + frame_num - 1;
            display_warning(sprintf("Getting template from image %d frame %d", image_idx, frame_num));
            setFrame(frame_num);
            windowSize = 100;
            x_offset = 5;
            y_offset = 10;
            h = drawrectangle('Parent', ax1,'Position',[800-windowSize-x_offset,800-windowSize-y_offset,windowSize,windowSize]);
            if (mode == "manual")
                % display_warning("select a template by drawing a rectangle");
                display_warning("Press enter to confirm the template");
                wait(h);
            end
            % Wait for the user to press the Enter key to confirm the template
            position = round(h.Position);
            template = imcrop(mat2gray(stack_info.img_data.imgs{image_idx}), position);
        end
    end
    function align_all_stacks_callback(~,~)
        % iterate over all the stacks and align them
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
            load_images_callback();
            if stack_info.aligned == true
                logs{end+1} = sprintf("Trial %s already aligned",stack_paths(i));
                continue;
            end
            if contains(path, 'time_control')
                stack_info.aligned = true;
                save_stack_callback();
                continue;
            end
            align_stack_callback('mode', 'auto');          
        end
    end
    function plot_displacements()
        displacements = stack_info.displacements;
        % clear axis
        cla(ax2);
        % plot the displacements on ax2
        plot(ax2, displacements(:,1), 'r');
        hold(ax2, 'on');
        plot(ax2, displacements(:,2), 'b');
        % draw vertical lines at start and end indices
        if stack_info.shortened == true
            % disp(stack_info.start_index);
            plot(ax2, stack_info.start_index, [-5, 5], 'g');
            plot(ax2, stack_info.end_index, [-5, 5], 'g');
        end
        title('Displacements');
        xlabel('Image number');
        ylabel('Displacement');
        legend('x', 'y');
        axis(ax2, 'tight'); 
    end
%%%%%%%%%%%%%%%%%%%%%% MASKING %%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        left_triangle = drawpolygon(ax1, 'Position', left_vertices);
        right_triangle = drawpolygon(ax1, 'Position', right_vertices);
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
        setFrame(get(slider, 'Value'));
        fprintf('Masked stack %s\n', path);
    end
    function draw_all_masks_callback(~,~)
        % iterate over all the stacks and align them
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
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
        set(slider, 'Min', 1, 'Max', end_index - start_index + 1, 'Value', 1, ...
            'SliderStep', [1/(end_index-start_index) , 1/(end_index-start_index)], 'Visible', 'On');
        set(to_frame, 'String', num2str(end_index - start_index + 1));
        setFrame(1);
    end
    function change_drive_callback(~, ~)
        % change the drive letter
        current_drive = 'E:';
        new_drive = 'F:';
        WaitMessage = parfor_wait(length(stack_paths), 'Waitbar', true);
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
            current_idx = get(stack_dropdown, 'Value');
            path = stack_paths{current_idx};
            update_info(path);
            [iteration, parentDir] = getIteration(path);
            if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
                while isfield(stack_info, 'stack_info')
                    stack_info = stack_info.stack_info;
                end
                % replace the parentDir with the new drive letter
                stack_info.parentDir = strrep(stack_info.parentDir, current_drive, new_drive);
                % replace the drive letter in each stack_info.img_data.img_files.folder
                for j = 1:numel(stack_info.img_data.img_files)
                    stack_info.img_data.img_files(j).folder = strrep(stack_info.img_data.img_files(j).folder, current_drive, new_drive);
                end
                save_stack_callback();
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
    function combine_stacks_callback()        
        % show popup to select the two stacks
        [primary_path, secondary_path] = select_stacks();
        if isempty(primary_path) || isempty(secondary_path)
            return;
        end
        % combine the stacks
        combine_stacks(primary_path, secondary_path);
        % update the stack paths
        stack_paths = get_stack_paths();
    end
    function [primary_path, secondary_path] = select_stacks()
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
    function combine_stacks(primary_path,secondary_path)
        fprintf('Combining %s and %s\n', primary_path, secondary_path);
        % get the last image name in the non cont path
        cont_img_files = dir(fullfile(secondary_path, '*.TIF'));
        last_img_name = get_last_image_name(primary_path);
        parts = split(last_img_name, '_');
        last_number_string = strrep(parts{end}, '.TIF', '');
        prefix = parts{1};
        % rename the images in the secondary path        
        wait = waitbar(0, 'Renaming and moving images');
        for i = 1:numel(cont_img_files)
            img_name = cont_img_files(i).name;
            new_number_string = incrementPaddedString(last_number_string, i);
            new_name = sprintf('%s_%s.TIF', prefix, new_number_string);
            fprintf('Copying images from %s to %s\n', fullfile(secondary_path, img_name), fullfile(primary_path, new_name));
            copyfile(fullfile(secondary_path, img_name), fullfile(primary_path, new_name));
            waitbar(i/numel(cont_img_files));
        end
        close(wait);
        function last_img_name = get_last_image_name(path)
            % get the last image name in the path
            img_files = dir(fullfile(path, '*.TIF'));
            last_img_name = img_files(end).name;
        end
    end
    function paddedResult = incrementPaddedString(paddedString, increment)
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
    function load_stack_info(~,~)
        current_idx = get(stack_dropdown, 'Value');
        path = stack_paths{current_idx};
        [iteration, parentDir] = getIteration(path);

        if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
            stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
            while isfield(stack_info, 'stack_info')
                stack_info = stack_info.stack_info;
            end
        end
    end
    function setFrame(k)
        persistent counter;
        % initialize the counter
        if isempty(counter)
            counter = 0;
        end
        % reset imgs after 100 frames
        if counter > 100
            % display_warning("Resetting images");
            counter = 0;
            stack_info.img_data.imgs = cell(1, stack_info.img_data.num_imgs);
        end
        image_idx = stack_info.start_index + k - 1;
        % fprintf('Frame number: %d\n image idx %d', k, image_idx);
        % Get the current slider value
        set(frame_number, 'String', num2str(k)); % Update the frame number display
        set(slider, 'Value', k);  % update the slider value
        % load the image if it's not already loaded
        if isempty(stack_info.img_data.imgs{image_idx})            
            % fprintf('Loading image %d as its empty', image_idx);
            stack_info.img_data.imgs{image_idx} = mat2gray(imread(fullfile(stack_info.img_data.img_files(image_idx).folder, stack_info.img_data.img_files(image_idx).name)));
            counter = counter + 1;
        end
        % if aligned displacements exist, apply them to the image
        if stack_info.aligned == true
            displaced_img = imtranslate(stack_info.img_data.imgs{image_idx}, - stack_info.displacements(image_idx, :));
            if isfield(stack_info, 'masked') && stack_info.masked == true && isfield(stack_info, 'mask')
                displaced_img = imoverlay(displaced_img, stack_info.mask, 'r');
            end 
            imshow(displaced_img, 'Parent', ax1);
        else
            imshow(stack_info.img_data.imgs{image_idx}, 'Parent', ax1);
        end
        if isfield(stack_info, 'timestamps')
            timestamp = stack_info.timestamps{image_idx};
            if ~isempty(timestamp)
                % draw the timestamp on the image
                % disp(image_idx);
                text(ax1, 'Units', 'normalized', 'Position', [0.99, 0.03], ...
                    'String', sprintf("%.2f Sec", ...
                    time_2_sec(timestamp)-time_2_sec(stack_info.timestamps{stack_info.start_index})), ...
                    'Color', 'white', 'FontSize', 18, 'HorizontalAlignment', 'right');
            end
        end
    end
    function indicator = create_indicator(parent, position, label, delete_callback)
        indicator = uicontrol(parent, 'Style', 'text', 'String', label, ...
        'Units', 'normalized', 'Position', [position(1)+0.1 position(2) position(3) position(4)], 'BackgroundColor', 'red');
        % if delete_callback is provided, set the callback for the indicator
        if exist('delete_callback', 'var')
            % add a dustbin icon to delete the indicator
            dustbin = imread('./bin_icon.png');
            dustbin = imresize(dustbin, [20, 20]);
            uicontrol(parent, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [position(1) + position(3) + 0.05 position(2) 0.08 0.05], 'CData', dustbin, ...
                'Callback', delete_callback);
        end
    end
    function toggle_indicator(indicator, status)
        if status
            set(indicator, 'BackgroundColor', 'green');
        else
            set(indicator, 'BackgroundColor', 'red');
        end
    end
    % function to zoom in and out at x,y with scroll
    function zoom_callback(~, event,x,y)
        % Get the scroll direction
        if event.VerticalScrollCount > 0
            % Zoom in
            zoomFactor = 1.1;
        else
            % Zoom out
            zoomFactor = 1/1.1;
        end
        % Get the current limits
        xLim = xlim(ax1);
        yLim = ylim(ax1);
        % Get the current point
        % x = event.IntersectionPoint(1);
        % y = event.IntersectionPoint(2);
        % Calculate the new limits
        xLim = (xLim - x) * zoomFactor + x;
        yLim = (yLim - y) * zoomFactor + y;
        % Set the new limits
        xlim(ax1, xLim);
        ylim(ax1, yLim);
    end
    function execute_selected_function(~, ~)
        % Get the selected function from the dropdown
        selected_index = get(function_dropdown, 'Value');
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
            case 'Get Scales'
                get_scales();
            case 'Plot scales'
                plot_scales();
            case 'Local number density'
                plot_local_number_density();
            case 'Average timeStamps'
                average_timeStamps();
            case 'Combine stacks'
                combine_stacks_callback();
            case 'Plot durations'
                plot_stack_durations();
            case 'Set all Empty or Not'
                set_all_empty_or_not();
        end
        
        display_warning(['Executed: ' selected_function]);
    end
    function forced_callback(~, ~)
        forced = ~forced;
    end
    function next_stack_callback(~, ~)
        if get(stack_dropdown, 'Value') < length(get(stack_dropdown, 'String'))
            set(stack_dropdown, 'Value', get(stack_dropdown, 'Value') + 1);
            load_images_callback();
        else
            display_warning("You're on the last stack");
        end
    end
    function prev_stack_callback(~, ~)
        if get(stack_dropdown, 'Value') < length(get(stack_dropdown, 'String'))
            set(stack_dropdown, 'Value', get(stack_dropdown, 'Value') - 1);
            load_images_callback();
        else
            display_warning("You're on the last stack");
        end
    end
    function goto_callback(~, ~)
        % go to the specified stack
        N = str2double(get(N_info, 'String'));
        f = str2double(get(f_info, 'String'));
        i = str2double(get(i_info, 'String'));
        path = sprintf('F:\\shake_table_data\\N%d\\%dhz_hopperflow\\60deg\\10cm\\%d', N, f, i);
        % get index of the matching path from stack_paths
        idx = find(contains(stack_paths, path));
        if isempty(idx)
            display_warning("Invalid stack number");
        else
            set(stack_dropdown, 'Value', idx);
            load_images_callback();
        end
    end
    function play_pause_callback(~, ~)
        if is_playing
            stop(play_timer);
            set(play_pause_button, 'CData', play_icon);
        else
            start(play_timer);
            set(play_pause_button, 'CData', pause_icon);
        end
        is_playing = ~is_playing;
        % fprintf('speed is : %d', speed);
    end
    function play_timer_callback(~, ~)
        current_value = get(slider, 'Value');
        if current_value < get(slider, 'Max')
            set(slider, 'Value', current_value + speed);
            slider_callback(slider);
        else
            stop(play_timer);
            set(play_pause_button, 'CData', play_icon);
            is_playing = false;
        end
    end
    function open_directory_callback(~, ~)
        % open the directory in windows explorer
        path = stack_paths{get(stack_dropdown, 'Value')};
        [~, parentDir] = getIteration(path);
        winopen(parentDir);
    end
    function slider_callback(~, ~)
        if ~evalin('base', 'exist(''stack_info'', ''var'')')
            display_warning("load some images first");
        else        
            % Get the current slider value
            slider_value = round(get(slider, 'Value'));
            % if isfield(stack_info, 'dead_zone')
            %     if slider_value >= stack_info.dead_zone(1) && slider_value <= stack_info.dead_zone(2)
            %         display_warning("Dead zone, skipping");
            %         setFrame(stack_info.dead_zone(2) + 1);
            %         set(slider, 'Value', stack_info.dead_zone(2) + 1);
            %     end
            % end
            setFrame(slider_value)
        end
    end
    function set_speed_callback(~,~)
        speed_idx = get(speed_dropdown, 'Value');
        speed = speeds{speed_idx};
    end
    function toggle_get_time_ui()
        if contains(path, 'time_control')
            get_time = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Get time', ...
                'Units', 'normalized','Position', [0.2 0.1 0.6 0.06], 'Callback', @get_times, 'Enable', 'on');
            % % add forced checkbox beside get_time
            % uicontrol(buttonPanel, 'Style', 'checkbox', 'String', 'Forced', ...
            %     'Units', 'normalized', 'Position', [0.8 0.1 0.2 0.06]);
        else
            if ~isempty(get_time)
                delete(get_time);
            end
        end
    end
    function update_info(path)
        % cancel any ongoing operation
        [iteration, ~] = getIteration(path);
        [N, fs] = get_info(path);
        set(N_info, 'String', N);
        set(f_info, 'String', fs);
        set(i_info, 'String', iteration);
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
    
        current_value = get(slider, 'Value');
        if event.VerticalScrollCount > 0  % if scrolling down
            new_value = current_value - step_size; % decrease the value
        else  % if scrolling up
            new_value = current_value + step_size; % increase the value
        end
        new_value = max(min(new_value, get(slider, 'Max')), get(slider, 'Min')); % ensure the new value is within the slider's range
        set(slider, 'Value', new_value);  % update the slider value
        slider_callback(slider);  % call the slider's callback function to update the display
    end
    function show_logs_callback(~, ~)
        % Create a new figure for the logs overlay
        log_fig = figure('Name', 'Logs', 'NumberTitle', 'Off', 'Position', [1100 100 400 600]);
        % Create a listbox to display the logs
        uicontrol('Style', 'listbox', 'Parent', log_fig, 'Units', 'normalized', ...
                    'Position', [0 0 1 1], 'String', logs, 'FontSize', 12);
    end
%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%
    % function to return a unique color from jet colormap for each N
    function color = get_color(N)
        % get the color from the jet colormap
        color = jet(50);
        color = color(N, :);
        % fprintf('Color for N = %d is %s\n', N, mat2str(color));
    end
    function marker = get_marker(N)
        markers = ['x', '+', 'h', '.', 'o', 's', 'd', '^', 'v', '>', '<', 'p', '*'];
        marker = markers(mod(N, numel(markers)) + 1);
        % fprintf('Marker for N = %d is %s\n',N, marker);
    end
    function wait_for_keypress(key_to_wait_for)
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
    function [trial_name, parentDir] = getIteration(path)
        if isa(path, 'char')
            path = string(path);
        end
        parts = path.split("\");
        trial_name = parts(end);
        parentDir = strjoin(parts(1:end-1), "\");
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
            path = stack_paths{get(stack_dropdown, 'Value')};
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
        path = stack_paths{get(stack_dropdown, 'Value')};
        [iteration, parentDir] = getIteration(path);
        display_warning("Initializing stack info");
        stack_info = struct('start_index', 1, 'end_index', numel(img_data.img_files),...
            'parentDir', parentDir, 'iteration', iteration, ...
            'aligned', false, 'shortened', false, 'masked', false,...
            'displacements', zeros(numel(img_data.img_files), 2), 'img_data', img_data);
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