clc;
close all;
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
% workspace;  % Make sure the workspace panel is showing.
% start parallel pool
% parpool(4);
% directory of the tif stacks
folder = 'E:\shake_table_data\';

% populate the list of paths to the tiff stacks
Ns = [4,12,24,48];
fs = [10,12,14,16,18,20];
deg = 60;
wd = 10;
stack_paths = [];

for n = 1:length(Ns)
    for freq = 1:length(fs)
        for w = 1:length(wd)
            % Construct the path with escaped backslashes
            path = sprintf("E:\\shake_table_data\\N%d\\%dhz_hopperflow\\%ddeg\\%dcm\\",Ns(n),fs(freq),deg,wd);
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
start_image_viewer(stack_paths);
% function to open the given list of tif images in a scrollable window
function start_image_viewer(stack_paths)
    close all;
    % Create figure and panel on it
    f = figure('Name', 'Simple Image Viewer', 'NumberTitle', 'Off','WindowState','maximized', 'WindowScrollWheelFcn', @scrollWheelMoved);
    % Create panels on the figure
    buttonPanel = uipanel(f, 'Units', 'normalized', 'Position', [0 0 0.2 1]);%[left bottom width height]
    axesPanel = uipanel(f, 'Units', 'normalized', 'Position', [0.2 0 0.8 1]);
    searchWindow = 50;
    folder_ico = imread('./folder.png');
    folder_ico = imresize(folder_ico, [20, 20]);
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
    shortenStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Shorten stack', ...
        'Units', 'normalized','Position', [0.2 0.58 0.6 0.06], 'Callback', @shorten_stack_callback, 'Enable', 'on');
    alignStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Align stack', ...
        'Units', 'normalized','Position', [0.2 0.5 0.4 0.06], 'Callback', @align_stack_callback, 'Enable', 'on');
    alignAllStack_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Align all', ...
        'Units', 'normalized','Position', [0.6 0.5 0.2 0.06], 'Callback', @align_all_stacks_callback, 'Enable', 'on');
    logs_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Logs', ...
        'Units', 'normalized','Position', [0.2 0.42 0.6 0.06], 'Callback', @show_logs_callback, 'Enable', 'on');
    save_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Save', ...
        'Units', 'normalized','Position', [0.2 0.32 0.6 0.06], 'Callback', @save_stack_callback, 'Enable', 'on');
    uicontrol(buttonPanel, 'Style', 'text', 'String', 'From Frame:', ...
        'Units', 'normalized', 'Position', [0.1 0.22 0.4 0.08]);    %[left bottom width height]
    uicontrol(buttonPanel, 'Style', 'text', 'String', 'To Frame:', ...
        'Units', 'normalized', 'Position', [0.5 0.22 0.4 0.08]);
    from_frame = uicontrol(buttonPanel, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.2 0.22 0.2 0.04], 'String', '1');
    to_frame = uicontrol(buttonPanel, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.6 0.22 0.2 0.04], 'String', '1');
    frame_number = uicontrol('Parent', f, 'Style', 'text', 'Position', [530 60 50 20], 'String', '1','FontSize', 14);
    slider = uicontrol(axesPanel, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.1 0.06 0.6 0.08], ...
        'Callback', @slider_callback, 'Visible', 'on');
    
    % Create indicators for aligned and shortened statuses
    aligned_indicator = create_indicator(buttonPanel, [0.2 0.8 0.1 0.05], "Aligned", @remove_alignment_callback);
    shortened_indicator = create_indicator(buttonPanel, [0.2 0.7 0.1 0.05], "Shortened");
    logs = {}; % Initialize logs list
    stack_info = struct();
    % load_images_callback();
    function indicator = create_indicator(parent, position, label, delete_callback)
        indicator = uicontrol(parent, 'Style', 'text', 'String', label, ...
        'Units', 'normalized', 'Position', [position(1)+0.1 position(2) 0.4 0.05], 'BackgroundColor', 'red');
        % if delete_callback is provided, set the callback for the indicator
        if exist('delete_callback', 'var')
            % add a dustbin icon to delete the indicator
            dustbin = imread('./bin_icon.png');
            dustbin = imresize(dustbin, [20, 20]);
            uicontrol(parent, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [position(1)+0.5 position(2) 0.08 0.05], 'CData', dustbin, ...
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
    function load_images_callback(~, ~)
        path = stack_paths{get(stack_dropdown, 'Value')};
        [iteration, parentDir] = getIteration(path);
        % if path is null load the first one
        if isempty(path)
            path = stack_paths{1};
        end

        if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
            stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
            assignin('base', 'stack_info', stack_info);
            stack_info = stack_info.stack_info;
        else
            img_data.img_files = dir(fullfile(path, '*.tif'));
            stack_info = initialize_stack_info(img_data);
        end

        % if start and end indices are set, set the shortened indicator to green
        if stack_info.shortened == true
            toggle_indicator(shortened_indicator, true);
            % Load the images
            stack_info.img_data.num_imgs = numel(stack_info.img_data.img_files(stack_info.start_index:stack_info.end_index));        
        else
            toggle_indicator(shortened_indicator, false);
            stack_info.img_data.num_imgs = numel(stack_info.img_data.img_files);
        end
        stack_info.img_data.imgs = cell(1, stack_info.img_data.num_imgs);

        % if displacement_n.mat exists, set the aligned indicator to green
        if stack_info.aligned == true && ~isempty(stack_info.displacements)
            toggle_indicator(aligned_indicator, true);
            % plot the displacements on ax2
            plot_displacements();
        else
            toggle_indicator(aligned_indicator, false);
        end

        set(to_frame, 'String', num2str(stack_info.img_data.num_imgs));
        % Initialize the slider
        set(slider, 'Min', 1, 'Max', stack_info.img_data.num_imgs, 'Value', 1, ...
            'SliderStep', [1/(stack_info.img_data.num_imgs-1) , 1/(stack_info.img_data.num_imgs-1)], 'Visible', 'On');
            
        % Display the first image
        setFrame(1);
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
            disp(stack_info.start_index);
            plot(ax2, stack_info.start_index, [-5, 5], 'g');
            plot(ax2, stack_info.end_index, [-5, 5], 'g');
        end
        title('Displacements');
        xlabel('Frame number');
        ylabel('Displacement');
        legend('x', 'y');
    end
    function align_stack_callback(varargin)
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
            if stack_info.shortened == true
                end_ = stack_info.end_index;
            else
                end_ = stack_info.img_data.num_imgs;
            end
            while all(template(:) == 0)
                frame_num = frame_num + 1;
                shorten_stack_callback(0, 0, frame_num, end_);
                [template,templatePosition] = getTemplate(mode, frame_num);
            end
            % match the template with each image in the stack
            displacements = matchTemplate(template, templatePosition);
            displacements(:,1) = displacements(:,1)-displacements(1,1);
            displacements(:,2) = displacements(:,2)-displacements(1,2);
            stack_info.displacements = displacements;
            % Display the displacements
            plot_displacements();
            % save displacements for later use
            stack_info.aligned = true;
            toggle_indicator(aligned_indicator, true);
            save_stack_callback();
        end
    end
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
    function displacements = matchTemplate(template, templatePosition)
        displacements = zeros(stack_info.img_data.num_imgs, 2);
        % f = waitbar(0,'Please wait...','Name','Aligning stack...');
        WaitMessage = parfor_wait(stack_info.img_data.num_imgs, 'Waitbar', true);
        img_data = stack_info.img_data;
        parfor k = 1:stack_info.img_data.num_imgs
            % setFrame(k);
            % display_warning(num2str(k));
            baseFileName = img_data.img_files(k).name;
            fullFileName = fullfile(img_data.img_files(k).folder, baseFileName);
            bwImage = imread(fullFileName);
            % Extract the region of interest (ROI) with `n` pixel padding
            x1 = max(1, round(templatePosition(1)) - searchWindow); % X start
            y1 = max(1, round(templatePosition(2)) - searchWindow); % Y start
            x2 = min(size(bwImage, 2), round(templatePosition(1) + templatePosition(3)) + searchWindow);
            y2 = min(size(bwImage, 1), round(templatePosition(2) + templatePosition(4)) + searchWindow);

            % Crop the padded region
            croppedImage = bwImage(y1:y2, x1:x2);
            % find the template in the image
            c = normxcorr2(template, croppedImage);
            [ypeak, xpeak] = find(c==max(c(:)));
            yoffSet = ypeak-size(template,1);
            xoffSet = xpeak-size(template,2);
            displacements(k, :) = [xoffSet, yoffSet];
            % Update the progress bar
            % progress = k / img_data.num_imgs;
            % waitbar(progress,f,sprintf("%d %% done",int32(progress*100)));
            WaitMessage.Send;
        end
        WaitMessage.Destroy
    end
    function [template,position] = getTemplate(mode, frame_num)
        % set slider to first image
        setFrame(frame_num);
        windowSize = 100;
        h = drawrectangle('Parent', ax1,'Position',[800-windowSize,800-windowSize,windowSize,windowSize]);
        if (mode == "manual")
            % display_warning("select a template by drawing a rectangle");
            display_warning("Press enter to confirm the template");
            wait(h);
        end
        % Wait for the user to press the Enter key to confirm the template
        position = round(h.Position);
        template = imcrop(mat2gray(stack_info.img_data.imgs{frame_num}), position);
    end
    function open_directory_callback(~, ~)
        % open the directory in windows explorer
        winopen(parentDir);
    end

    function setFrame(k)
        % Get the current slider value
        set(frame_number, 'String', num2str(k)); % Update the frame number display
        set(slider, 'Value', k);  % update the slider value
        % Display the corresponding image
        if isempty(stack_info.img_data.imgs{k})
            stack_info.img_data.imgs{k} = mat2gray(imread(fullfile(stack_info.img_data.img_files(k).folder, stack_info.img_data.img_files(k).name)));
        end
        % if aligned displacements exist, apply them to the image
        if exist('displacements', 'var')
            displaced_img = imtranslate(stack_info.img_data.imgs{k}, -stack_info.displacements(k, :));
            imshow(displaced_img, 'Parent', ax1);
        else
            imshow(stack_info.img_data.imgs{k}, 'Parent', ax1);
        end
    end
    function show_logs_callback(~, ~)
        % Create a new figure for the logs overlay
        log_fig = figure('Name', 'Logs', 'NumberTitle', 'Off', 'Position', [1100 100 400 600]);
        % Create a listbox to display the logs
        uicontrol('Style', 'listbox', 'Parent', log_fig, 'Units', 'normalized', ...
                    'Position', [0 0 1 1], 'String', logs, 'FontSize', 12);
    end
    function [start_index, end_index] = shorten_stack_callback(~, ~, start_, end_)
        setFrame(1);
        if nargin > 2
            start_index = start_;
            end_index = end_;
        else
            display_warning("go to start frame and press n");
            wait_for_keypress("n");
            start_index = round(get(slider, 'Value'));
            display_warning("go to end frame and press n");
            setFrame(stack_info.img_data.num_imgs);
            wait_for_keypress("n");
            end_index = round(get(slider, 'Value'));
        end
        shorten_stack(start_index, end_index);
        stack_info.shortened = true;
        toggle_indicator(shortened_indicator, true);
        % save the start and end indices in a stack_info.mat file
        save_stack_callback();

    end
    function shorten_stack(start_index, end_index)
        % shorten stack
        stack_info.start_index = start_index;
        stack_info.end_index = end_index;
        set(slider, 'Min', start_index, 'Max', end_index, 'Value', 1, ...
            'SliderStep', [1/(end_index-start_index) , 1/(end_index-start_index)], 'Visible', 'On');
        set(to_frame, 'String', num2str(end_index));
        setFrame(1);
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
        if isempty(last_scroll_time)
            last_scroll_time = now;  % initialize to the current time
        end
        time_between_scrolls = (now - last_scroll_time) * 24 * 60 * 60;  % in seconds
        last_scroll_time = now;  % update for next time
    
        step_size = max(1, round(1 / time_between_scrolls*0.1));  % larger step size for smaller time_between_scrolls
    
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
    function slider_callback(~, ~)
        if(stack_info.img_data.num_imgs == 0)
            display_warning("load some images first");
        else        
            % Get the current slider value
            slider_value = round(get(slider, 'Value'));
            setFrame(slider_value)
        end
    end
    function wait_for_keypress(key_to_wait_for)
        keypressed = 0;  % UserData is 0 before key press
        set(f, 'KeyPressFcn', @myKeyPressFcn)

        function myKeyPressFcn(src, event)
            if strcmp(event.Key, key_to_wait_for)
                keypressed = 1;  % Set UserData to 1 after key press
            end
        end

        % Wait until UserData is 1
        while keypressed == 0
            pause(0.1);  % Short pause to avoid overloading the CPU
        end
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
            temp_stack_info = stack_info;
            % get the current stack path
            path = stack_paths{get(stack_dropdown, 'Value')};
            [iteration, parentDir] = getIteration(path);
            % remove the image data from the stack_info
            stack_info.img_data.imgs = [];
            % save the stack
            save(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'stack_info');
            display_warning(sprintf("Stack saved to %s//stack_info_%s.mat", parentDir, iteration));
            save_button.String = 'Save';
            stack_info = temp_stack_info;
        end
    end
    function stack_info = initialize_stack_info(img_data)
        % Image data
        display_warning("Initializing stack info");
        stack_info = struct('start_index', 1, 'end_index', numel(img_data.img_files),...
            'parentDir', '', 'iteration', '', ...
            'aligned', false, 'shortened', false,...
            'displacements', zeros(1, 2), 'img_data', img_data);
    end
    function align_all_stacks_callback(~,~)
        % iterate over all the stacks and align them
        for i = 1:length(stack_paths)
            set(stack_dropdown, 'Value', i);
            load_images_callback();
            if exist(sprintf('%s//displacements_%s.mat', parentDir, iteration), 'file')
                logs{end+1} = sprintf("%s Trial %s already aligned",parentDir, iteration);
                continue;
            end
            align_stack_callback('mode', 'auto');
        end
    end
end