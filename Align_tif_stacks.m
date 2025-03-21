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
% TODO: Plot all timestamps
% TODO: 


% function to open the given list of tif images in a scrollable window
function start_image_viewer(stack_paths)
    close all;

    searchWindow = 50;
    function_list = {'Change Drive', 'Plot all Gr', 'Plot all TimeStamps'};
    get_time = [];path = [];
    is_playing = false;
    speed = 1;
    logs = {}; % Initialize logs list
    stack_info = struct();
    skip_alignment = false;
    forced = false;
    ui = ui_builder();
    goto_callback();
    



    % function to start getting the times for every nth image and also has the option to stop the current operation

    function remove_shortening_callback(~, ~)
        % remove the displacements file
        stack_info.shortened = false;
        stack_info.start_index = 1;
        stack_info.end_index = numel(stack_info.img_data.img_files);
        shorten_slider(1, numel(stack_info.img_data.img_files));
        toggle_indicator(ui.shortened_indicator, false);
        % get current frame
        frame = get(ui.slider, 'Value');
        setFrame(frame);
        % add star to save button
        save_button.String = 'Save *';
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
        set(ui.frame_number, 'String', num2str(k)); % Update the frame number display
        set(ui.slider, 'Value', k);  % update the slider value
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
            imshow(displaced_img, 'Parent', ui.ax1);
        else
            imshow(stack_info.img_data.imgs{image_idx}, 'Parent', ui.ax1);
        end
        if isfield(stack_info, 'timestamps')
            timestamp = stack_info.timestamps{image_idx};
            if ~isempty(timestamp)
                % draw the timestamp on the image
                % disp(image_idx);
                text(ui.ax1, 'Units', 'normalized', 'Position', [0.99, 0.03], ...
                    'String', sprintf("%.2f Sec", ...
                    time_2_sec(timestamp)-time_2_sec(stack_info.timestamps{1})), ...
                    'Color', 'white', 'FontSize', 18, 'HorizontalAlignment', 'right');
            end
        end
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
            start_index = stack_info.start_index + round(get(ui.slider, 'Value')) - 1;
            display_warning("go to end frame and press n");
            setFrame(stack_info.end_index-stack_info.start_index + 1);
            wait_for_keypress("n");
            end_index = stack_info.start_index + round(get(ui.slider, 'Value')) - 1;
        end
        shorten_slider(start_index, end_index);
        stack_info.shortened = true;
        toggle_indicator(ui.shortened_indicator, true);
        % save the start and end indices in a stack_info.mat file
        save_stack_callback();
    end
    function shorten_all_stack_callback(~,~)
        % iterate over all the stacks and align them
        for i = 1:length(stack_paths)
            set(ui.stack_dropdown, 'Value', i);
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

        left_triangle = drawpolygon(ui.ax1, 'Position', left_vertices);
        right_triangle = drawpolygon(ui.ax1, 'Position', right_vertices);
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
        setFrame(get(ui.slider, 'Value'));
        fprintf('Masked stack %s\n', path);
    end
    function draw_all_masks_callback(~,~)
        % iterate over all the stacks and align them
        for i = 1:length(stack_paths)
            set(ui.stack_dropdown, 'Value', i);
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
        set(ui.slider, 'Min', 1, 'Max', end_index - start_index + 1, 'Value', 1, ...
            'SliderStep', [1/(end_index-start_index) , 1/(end_index-start_index)], 'Visible', 'On');
        set(ui.to_frame, 'String', num2str(end_index - start_index + 1));
        setFrame(1);
    end

    % function to return a unique color from jet colormap for each N

    % function mark_dead_zone_callback(~,~)
    %     % get the start and end frame of the dead timeline and save it in the stack_info
    %     display_warning("Select the start frame of the dead timeline");
    %     wait_for_keypress("n");
    %     start_frame = round(get(ui.slider, 'Value'));
    %     display_warning("Select the end frame of the dead timeline");
    %     wait_for_keypress("n");
    %     end_frame = round(get(ui.slider, 'Value'));   
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