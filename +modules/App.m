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
        aligner
        shortener
        rdf
        timer
        masking
        steinhardt
        trial
        voids
        particle_locator

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
        current_image_idx
        monitorChoice = 2;
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
            app.stack_paths = app.get_stack_paths();

            app.utils = modules.Utils(app);
            app.voids = modules.Voids(app);
            app.timer = modules.Timer(app);
            app.trial = modules.Trial(app);
            app.steinhardt = modules.Steinhardt(app);
            app.masking = modules.Masking(app);
            app.aligner = modules.Aligner(app);
            app.shortener = modules.Shortener(app);
            app.rdf = modules.RDF(app);
            app.particle_locator = modules.Particle_locator(app);
            app.ui = modules.Ui(app);

            % load the first stack
            app.utils.goto_callback()
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
            obj.utils.toggle_get_time_ui();

            obj.utils.update_info(obj.path);
            [iteration, parentDir] = obj.utils.getIteration(obj.path);
            set(obj.ui.info.stackLabel, 'String', sprintf('Stack #%d of %d', current_idx, length(obj.stack_paths)));

            stack_info_path = sprintf('%s//stack_info_%s.mat', parentDir, iteration);

            if exist(stack_info_path, 'file')
                obj.stack_info = load(stack_info_path);
                while isfield(obj.stack_info, 'stack_info')
                    obj.stack_info = obj.stack_info.stack_info;
                end
                % obj.stack_info = modules.Stackinfo(stack_info_path);
                WaitMessage.Send;
                % assignin('base', 'stack_info', obj.stack_info);
            else
                img_data.img_files = dir(fullfile(obj.path, '*.tif'));
                obj.stack_info = obj.utils.initialize_stack_info(img_data);
            end
            
            obj.stack_info.img_data.num_imgs = numel(obj.stack_info.img_data.img_files);
            obj.stack_info.img_data.imgs = cell(1, obj.stack_info.img_data.num_imgs);
            WaitMessage.Send;
            % if start and end indices are set, set the shortened indicator to green
            if obj.stack_info.shortened == true
                obj.utils.toggle_indicator(obj.ui.info.shortenedIndicator, true);
            else
                obj.utils.toggle_indicator(obj.ui.info.shortenedIndicator, false);
            end
            WaitMessage.Send;
            % if displacement_n.mat exists, set the aligned indicator to green
            if obj.stack_info.aligned == true && ~isempty(obj.stack_info.displacements)
                obj.utils.toggle_indicator(obj.ui.info.alignedIndicator, true);
                % plot the displacements on obj.ui.controls.ax2
                obj.aligner.plot_displacements();
            else
                obj.utils.toggle_indicator(obj.ui.info.alignedIndicator, false);
                % clear axis
                cla(obj.ui.controls.ax2);
            end
            obj.masking.shorten_slider(obj.stack_info.start_index, obj.stack_info.end_index);
            WaitMessage.Send;
            fprintf('Loaded stack %s\n', obj.path);
            WaitMessage.Destroy;
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
        function paddedResult = incrementPaddedString(~, paddedString, increment)
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
    end
end