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
fs = [4,6,8,10,12,14,16,18,20];
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
        'Units', 'normalized','Position', [0.2 0.7 0.6 0.06], 'Callback', @draw_masks_callback, 'Enable', 'on');
    deadZone_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'Mark dead zone', ...
        'Units', 'normalized','Position', [0.2 0.64 0.6 0.06], 'Callback', @mark_dead_zone_callback, 'Enable', 'on');
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
        'Units', 'normalized', 'Position', [0.15 0.87 0.1 0.04], 'String', '1');
    f_info = uicontrol(buttonPanel, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.4 0.87 0.1 0.04], 'String', '1');
    i_info = uicontrol(buttonPanel, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.65 0.87 0.1 0.04], 'String', '1');
    goto_button = uicontrol(buttonPanel, 'Style', 'pushbutton', 'String', 'go to', ...
    'Units', 'normalized','Position', [0.77 0.87 0.15 0.04], 'Callback', @goto_callback);
    % Create play/pause button beside the scrollbar
    play_icon = imread('./play.png');
    play_icon = imresize(play_icon, [40, 40]);
    pause_icon = imread('./pause.png');
    pause_icon = imresize(pause_icon, [40, 40]);
    play_pause_button = uicontrol(axesPanel, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.01 0.06 0.08 0.08], 'CData', play_icon, 'Callback', @play_pause_callback);
    is_playing = false;
    play_timer = timer('ExecutionMode', 'fixedRate', 'Period', 0.1, 'TimerFcn', @play_timer_callback);
   
    % Create indicators for aligned and shortened statuses
    aligned_indicator = create_indicator(buttonPanel, [0 0.8 0.3 0.05], "Aligned", @remove_alignment_callback);
    shortened_indicator = create_indicator(buttonPanel, [0.5 0.8 0.3 0.05], "Shortened", @remove_shortening_callback);
    logs = {}; % Initialize logs list
    stack_info = struct();
    skip_alignment = false;

    function skip_alignment_callback(~, ~)
        % skip the current stack
        skip_alignment = true;
    end
    % load_images_callback();
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
    function load_images_callback(~, ~)
        current_idx = get(stack_dropdown, 'Value');
        path = stack_paths{current_idx};
        update_info(path);
        [iteration, parentDir] = getIteration(path);
        set(stack_label, 'String', sprintf('Stack #%d of %d', current_idx, length(stack_paths)));

        if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
            stack_info = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration));
            assignin('base', 'stack_info', stack_info);
            stack_info = stack_info.stack_info;
        else
            img_data.img_files = dir(fullfile(path, '*.tif'));
            stack_info = initialize_stack_info(img_data);
        end
        
        stack_info.img_data.num_imgs = numel(stack_info.img_data.img_files);
        stack_info.img_data.imgs = cell(1, stack_info.img_data.num_imgs);
        
        % if start and end indices are set, set the shortened indicator to green
        if stack_info.shortened == true
            toggle_indicator(shortened_indicator, true);
        else
            toggle_indicator(shortened_indicator, false);
        end

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
        xlabel('Image number');
        ylabel('Displacement');
        legend('x', 'y');
        axis(ax2, 'tight'); 
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
    function lastNonConstIndex = findLastNonConstIndex(arr)
        displacements = stack_info.displacements;
        displacements(:,1) = displacements(:,1) - displacements(end,1);
        displacements(:,2) = displacements(:,2) - displacements(end,2);
        % find the last non zero index
        lastNonConstIndex = find(displacements(:,2), 1, 'last');
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
            if isfield(stack_info, 'masked')
                if stack_info.masked == true
                    displaced_img = imoverlay(displaced_img, stack_info.mask, 'r');
                end
            end 
            imshow(displaced_img, 'Parent', ax1);
        else
            imshow(stack_info.img_data.imgs{image_idx}, 'Parent', ax1);
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
    function slider_callback(~, ~)
        if(stack_info.img_data.num_imgs == 0)
            display_warning("load some images first");
        else        
            % Get the current slider value
            slider_value = round(get(slider, 'Value'));
            if isfield(stack_info, 'dead_zone')
                if slider_value >= stack_info.dead_zone(1) && slider_value <= stack_info.dead_zone(2)
                    display_warning("Dead zone, skipping");
                    setFrame(stack_info.dead_zone(2) + 1);
                    set(slider, 'Value', stack_info.dead_zone(2) + 1);
                end
            end
            setFrame(slider_value)
        end
    end
    function goto_callback(~, ~)
        % go to the specified stack
        N = str2double(get(N_info, 'String'));
        f = str2double(get(f_info, 'String'));
        i = str2double(get(i_info, 'String'));
        path = sprintf('E:\\shake_table_data\\N%d\\%dhz_hopperflow\\60deg\\10cm\\%d', N, f, i);
        % get index of the matching path from stack_paths
        idx = find(contains(stack_paths, path));
        if isempty(idx)
            display_warning("Invalid stack number");
        else
            set(stack_dropdown, 'Value', idx);
            load_images_callback();
        end
    end
    function wait_for_keypress(key_to_wait_for)
        keypressed = 0;  % UserData is 0 before key press
        set(f, 'KeyPressFcn', @myKeyPressFcn)

        function myKeyPressFcn(~, event)
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
        end
    end
    function stack_info = initialize_stack_info(img_data)
        % Image data
        path = stack_paths{get(stack_dropdown, 'Value')};
        [iteration, parentDir] = getIteration(path);
        display_warning("Initializing stack info");
        stack_info = struct('start_index', 1, 'end_index', numel(img_data.img_files),...
            'parentDir', parentDir, 'iteration', iteration, ...
            'aligned', false, 'shortened', false,...
            'displacements', zeros(numel(img_data.img_files), 2), 'img_data', img_data);
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
            align_stack_callback('mode', 'auto');
        end
    end
    function open_directory_callback(~, ~)
        % open the directory in windows explorer
        path = stack_paths{get(stack_dropdown, 'Value')};
        [~, parentDir] = getIteration(path);
        winopen(parentDir);
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
    end
    function play_timer_callback(~, ~)
        current_value = get(slider, 'Value');
        if current_value < get(slider, 'Max')
            set(slider, 'Value', current_value + 1);
            slider_callback(slider);
        else
            stop(play_timer);
            set(play_pause_button, 'CData', play_icon);
            is_playing = false;
        end
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
    function mark_dead_zone_callback(~,~)
        % get the start and end frame of the dead timeline and save it in the stack_info
        display_warning("Select the start frame of the dead timeline");
        wait_for_keypress("n");
        start_frame = round(get(slider, 'Value'));
        display_warning("Select the end frame of the dead timeline");
        wait_for_keypress("n");
        end_frame = round(get(slider, 'Value'));   
        stack_info.dead_zone = [start_frame, end_frame];
        save_stack_callback();
    end
    function [N, fs] = get_info(path)
        if isa(path, 'char')
            path = string(path);
        end
        % E:\shake_table_data\N12\10hz_hopperflow\60deg\10cm\1
        parts = path.split("\");
        N = sscanf(parts{3}, 'N%d');
        fs = sscanf(parts{4}, '%dhz_hopperflow');
    end
    function update_info(path)
        % cancel any ongoing operation
        [iteration, ~] = getIteration(path);
        [N, fs] = get_info(path);
        set(N_info, 'String', N);
        set(f_info, 'String', fs);
        set(i_info, 'String', iteration);
    end
end