classdef Ui < handle
    properties
        app % Reference to the App instance
        fig
        buttonPanel
        axesPanel
        controls % Struct to hold UI controls
        info
        function_list
        folder_ico
        next_ico
        prev_ico
        play_icon
        pause_icon
        particle_icon
        dustbin_icon
    end

    methods
        function obj = Ui(app)
            obj.app = app;
            obj.function_list = {'Change Drive','Save all stacks', 'Plot all Gr', 'Plot all TimeStamps', 'Get Scale', 'Set Scales','Plot scales',...
                'Local number density','Average timeStamps','Combine stacks','Plot durations','Set all Empty or Not','Set all Jam or Not',...
                'Find all nhood', 'Find all Si6', 'Get LBOOP', 'Get all LBOOP','Plot LBOOP for all N','Plot avg LBOOP for all N','Plot BOOP',...
                'Superpose all LBOOP', 'Detect voids - stack', 'Detect voids - all stacks', 'Analyze voids - stack', 'Analyze voids - all stacks',...
                'Consolidate voids data','Visualize voids data','Create voids images'};
            obj.createUi();
        end

        function createUi(obj)
            close all;
            % Get information about all monitors
            monitorPositions = get(0, 'MonitorPositions'); 
            % Check if there's a second monitor
            if size(monitorPositions, 1) > 1 && obj.app.monitorChoice == 2
                % Get second monitor position
                secondMonitor = monitorPositions(2, :);
                
                % Calculate figure position on left side of second monitor
                figWidth = min(1200, secondMonitor(3) * 0.9);
                figHeight = min(800, secondMonitor(4) * 0.9);
                left = secondMonitor(1); % Left edge of second monitor
                bottom = secondMonitor(2) + (secondMonitor(4) - figHeight)/2; % Centered vertically
                
                % Create figure with specific position
                obj.fig = figure('Name', 'Stack Viewer', 'NumberTitle', 'Off', ...
                    'Position', [left, bottom, figWidth, figHeight], ...
                    'WindowScrollWheelFcn', @(src, event) obj.app.utils.scrollWheelMoved(src, event));
            else
                % Fallback if only one monitor
                obj.fig = figure('Name', 'Stack Viewer', 'NumberTitle', 'Off', ...
                    'WindowState', 'maximized', 'WindowScrollWheelFcn', @(src, event) obj.app.utils.scrollWheelMoved(src, event));
            end

            % Load icons
            obj.folder_ico = imresize(imread('./images/folder.png'), [20, 20]);
            obj.next_ico = imresize(imread('./images/next.png'), [20, 20]);
            obj.prev_ico = imresize(imread('./images/prev.png'), [20, 20]);
            obj.play_icon = imresize(imread('./images/play.png'), [40, 40]);
            obj.pause_icon = imresize(imread('./images/pause.png'), [40, 40]);
            obj.particle_icon = imresize(imread('./images/particles_icon.png'), [40, 40]);
            obj.dustbin_icon = imresize(imread('./images/bin_icon.png'), [20, 20]);

            % Create panels
            obj.buttonPanel = uipanel(obj.fig, 'Units', 'normalized', 'Position', [0 0 0.2 1]);
            obj.axesPanel = uipanel(obj.fig, 'Units', 'normalized', 'Position', [0.2 0 0.8 1]);
            % Axes and sliders
            obj.controls.ax1 = axes(obj.axesPanel, 'Position', [0.01 0.2 0.49 0.7]);
            obj.controls.ax2 = axes(obj.axesPanel, 'Position', [0.5 0.2 0.49 0.7]);

            % Create dropdown and buttons
            obj.controls.stackDropdown = uicontrol(obj.buttonPanel, 'Style', 'popupmenu', ...
                'String', obj.app.stack_paths, 'Units', 'normalized', ...
                'Position', [0.2 0.9 0.6 0.06], 'Callback', @(src, event) obj.app.load_images_callback());

            obj.controls.openDirButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', '', 'Units', 'normalized', 'Position', [0.8 0.93 0.08 0.03], ...
                'CData', obj.folder_ico, 'Callback', @(src, event) obj.app.utils.open_directory_callback());

            obj.controls.drawMasksButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Draw mask', 'Units', 'normalized', 'Position', [0.2 0.7 0.4 0.06], ...
                'Callback', @(src, event) obj.app.masking.draw_masks_callback());

            obj.controls.drawAllMasksButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Do all', 'Units', 'normalized', 'Position', [0.6 0.7 0.2 0.06], ...
                'Callback', @(src, event) obj.app.masking.draw_all_masks_callback());

            obj.controls.GrButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Get Gr', 'Units', 'normalized', 'Position', [0.2 0.64 0.4 0.06], ...
                'Callback', @(src, event) obj.app.rdf.get_Gr());

            obj.controls.GrAllStackButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Get all', 'Units', 'normalized', 'Position', [0.6 0.64 0.2 0.06], ...
                'Callback', @(src, event) obj.app.rdf.Gr_all_stacks_callback());

            obj.controls.shortenStackButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Shorten stack', 'Units', 'normalized', 'Position', [0.2 0.58 0.4 0.06], ...
                'Callback', @(src, event) obj.app.shortener.shorten_stack_callback());

            obj.controls.shortenAllStackButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Do all', 'Units', 'normalized', 'Position', [0.6 0.58 0.2 0.06], ...
                'Callback', @(src, event) obj.app.shortener.shorten_all_stack_callback());

            obj.controls.alignStackButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Align stack', 'Units', 'normalized', 'Position', [0.2 0.52 0.4 0.06], ...
                'Callback', @(src, event) obj.app.aligner.align_stack_callback());

            obj.controls.alignAllStackButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Align all', 'Units', 'normalized', 'Position', [0.6 0.52 0.2 0.06], ...
                'Callback', @(src, event) obj.app.aligner.align_all_stacks_callback());

            obj.controls.timestampStackButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Timestamp stack', 'Units', 'normalized', 'Position', [0.2 0.46 0.4 0.06], ...
                'Callback', @(src, event) obj.app.timer.timestamp_stack());

            obj.controls.timestampAllStackButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Do all', 'Units', 'normalized', 'Position', [0.6 0.46 0.2 0.06], ...
                'Callback', @(src, event) obj.app.timer.timestamp_all_stacks());

            obj.controls.logsButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Logs', 'Units', 'normalized', 'Position', [0.51 0.32 0.3 0.06], ...
                'Callback', @(src, event) obj.app.utils.show_logs_callback());

            obj.controls.saveButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', ...
                'String', 'Save', 'Units', 'normalized', 'Position', [0.2 0.32 0.3 0.06], ...
                'Callback', @(src, event) obj.app.utils.save_stack_callback());

            obj.info.fromFrameLabel = uicontrol(obj.buttonPanel, 'Style', 'text', 'String', 'From Frame:', ...
                'Units', 'normalized', 'Position', [0.1 0.22 0.4 0.08]);    %[left bottom width height]
            obj.info.toFrameLabel = uicontrol(obj.buttonPanel, 'Style', 'text', 'String', 'To Frame:', ...
                'Units', 'normalized', 'Position', [0.5 0.22 0.4 0.08]);
            obj.info.fromFrame = uicontrol(obj.buttonPanel, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [0.2 0.22 0.2 0.04], 'String', '1');
            obj.info.toFrame = uicontrol(obj.buttonPanel, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [0.6 0.22 0.2 0.04], 'String', '1');

            obj.controls.slider = uicontrol(obj.axesPanel, 'Style', 'slider', ...
                'Units', 'normalized', 'Position', [0.1 0.06 0.6 0.08], ...
                'Callback', @(src, event) obj.app.utils.slider_callback(), 'Visible', 'on');

            obj.controls.frameNumber = uicontrol(obj.axesPanel, 'Style', 'text', ...
                'Units', 'normalized', 'Position', [0.7 0.04 0.1 0.08], ...
                'String', '1', 'FontSize', 14);

            obj.controls.nextStackButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', 'String', '', ...
                'Units', 'normalized', 'Position', [0.1 0.93 0.08 0.03], 'CData', obj.next_ico,...
                'Callback', @(src, event) obj.app.utils.next_stack_callback, 'Enable', 'on');
            obj.controls.prevStackButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', 'String', '', ...
                'Units', 'normalized', 'Position', [0.03 0.93 0.08 0.03], 'CData', obj.prev_ico,...
                'Callback', @(src, event) obj.app.utils.prev_stack_callback, 'Enable', 'on');
            obj.controls.skipButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', 'String', 'Skip', ...
                'Units', 'normalized','Position', [0.2 0.15 0.6 0.06], 'Callback', @(src, event) obj.app.aligner.skip_alignment_callback, 'Enable', 'off');
            obj.info.stackLabel = uicontrol(obj.buttonPanel, 'Style', 'text', 'String', 'Stack #1', ...
                'Units', 'normalized', 'Position', [0.2 0.96 0.6 0.03]);

            % create info UI
            uicontrol(obj.buttonPanel, 'Style', 'text', 'String', 'N:', ...
                'Units', 'normalized', 'Position', [-0.1 0.82 0.4 0.08]);
            uicontrol(obj.buttonPanel, 'Style', 'text', 'String', 'f:', ...
                'Units', 'normalized', 'Position', [0.15 0.82 0.4 0.08]);
            uicontrol(obj.buttonPanel, 'Style', 'text', 'String', 'Iter:', ...
                'Units', 'normalized', 'Position', [0.4 0.82 0.4 0.08]);
            obj.info.N_info = uicontrol(obj.buttonPanel, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [0.15 0.87 0.1 0.04], 'String', '4');
            obj.info.f_info = uicontrol(obj.buttonPanel, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [0.4 0.87 0.1 0.04], 'String', '4');
            obj.info.i_info = uicontrol(obj.buttonPanel, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [0.65 0.87 0.1 0.04], 'String', '1');
            obj.controls.goto_button = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', 'String', 'load', ...
                'Units', 'normalized','Position', [0.77 0.87 0.15 0.04], 'Callback', @(src, event) obj.app.utils.goto_callback);

            obj.controls.functionDropdown = uicontrol(obj.buttonPanel, 'Style', 'popupmenu', ...
                'String', obj.function_list, ...
                'Units', 'normalized', 'Position', [0.2 0.04 0.6 0.06]);

            obj.controls.executeButton = uicontrol(obj.buttonPanel, 'Style', 'pushbutton', 'String', '', ...
                'Units', 'normalized', 'Position', [0.8 0.07 0.08 0.03], 'CData', obj.next_ico,...
                'Callback', @(src, event) obj.app.utils.execute_selected_function);

            obj.controls.playPauseButton = uicontrol(obj.axesPanel, 'Style', 'pushbutton', ...
                'Units', 'normalized', 'Position', [0.01 0.06 0.08 0.08], ...
                'CData', obj.play_icon, 'Callback', @(src, event) obj.app.utils.play_pause_callback());

            obj.controls.forcedCheckbox = uicontrol(obj.buttonPanel, 'Style', 'checkbox', 'String', 'Forced', ...
                'Units', 'normalized', 'Position', [0.8 0.7 0.2 0.06], 'Callback', @(src, event) obj.app.utils.forced_callback, 'Enable', 'on');
            obj.controls.speedDropdown = uicontrol(obj.axesPanel, 'Style', 'popupmenu', 'String', obj.app.speeds, ...
                'Units', 'normalized', 'Position', [0.01 0.001 0.08 0.08], 'Callback', @(src, event) obj.app.utils.set_speed_callback);

            obj.controls.particleButton = uicontrol(obj.axesPanel, 'Style', 'pushbutton', ...
                'CData', obj.particle_icon, 'Units', 'normalized', ...
                'Position', [0 0.95 0.04 0.05], 'Callback', @(src, event) obj.app.particle_locator.toggle_particle_locations());
            % Add UI elements for voids detection
            obj.app.ui.controls.particleButton = uicontrol(obj.axesPanel, 'Style', 'pushbutton', ...
                'CData', obj.particle_icon, 'Units', 'normalized', ...
                'Position', [0.05 0.95 0.04 0.05], 'Callback', @(src, event) obj.app.voids.detect_voids_callback());

            obj.info.playTimer = timer('ExecutionMode', 'fixedRate', 'Period', 0.1, 'TimerFcn', @(src, event) obj.app.utils.play_timer_callback);
            % Create indicators for aligned and shortened statuses
            obj.info.alignedIndicator = obj.create_indicator(obj.buttonPanel, [0 0.8 0.3 0.05], "Aligned");
            obj.add_delete_button(obj.buttonPanel, [0 0.8 0.3 0.05], @(src, event) obj.app.aligner.remove_alignment_callback);

            obj.info.shortenedIndicator = obj.create_indicator(obj.buttonPanel, [0.5 0.8 0.3 0.05], "Shortened");
            obj.add_delete_button(obj.buttonPanel, [0.5 0.8 0.3 0.05], @(src, event) obj.app.shortener.remove_shortening_callback);

            obj.controls.get_time = [];
        end
        function indicator = create_indicator(~, parent, position, label)
            % Create the indicator text
            indicator = uicontrol(parent, 'Style', 'text', 'String', label, ...
                'Units', 'normalized', 'Position', [position(1)+0.1 position(2) position(3) position(4)], ...
                'BackgroundColor', 'red');
        end

        function add_delete_button(obj, parent, position, callback)
            % Add a delete button next to the indicator
            uicontrol(parent, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [position(1) + position(3) + 0.05 position(2) 0.08 0.05], ...
                'CData', obj.dustbin_icon, 'Callback', callback);
        end
    end
end
