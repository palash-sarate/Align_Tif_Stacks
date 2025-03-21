classdef utils < handle
    methods
        % constructor
        function obj = utils()
        end
        function color = get_color(~, N, Ns)
            % get the number of unique Ns
            unique_Ns = unique(Ns);
            % get the index of the current N
            idx = find(unique_Ns == N);
            % get the color from the jet colormap
            color = jet(numel(unique_Ns));
            color = color(idx, :);
            % fprintf('Color for N = %d is %s\n', N, mat2str(color));
        end
        function out = pyList2cell(~,pyobj)
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
        function secs = time_2_sec(timestamp)
            secs = timestamp.min * 60 + timestamp.s + timestamp.us / 1e6;
        end
        function change_drive_callback(~, ~)
            % change the drive letter
            current_drive = 'E:';
            new_drive = 'F:';
            WaitMessage = parfor_wait(length(stack_paths), 'Waitbar', true);
            for i = 1:length(stack_paths)
                set(ui.stack_dropdown, 'Value', i);
                current_idx = get(ui.stack_dropdown, 'Value');
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
        function [trial_name, parentDir] = getIteration(path)
            if isa(path, 'char')
                path = string(path);
            end
            parts = path.split("\");
            trial_name = parts(end);
            parentDir = strjoin(parts(1:end-1), "\");
        end
    end
end