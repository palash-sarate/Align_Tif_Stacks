classdef Trial < handle
    properties
        app
    end
    methods
        function obj = Trial(app)
            obj.app = app;
        end
        %%%%%%%%%%%%%%%%%%%%%% SCALE %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function get_scales(obj, ~,~)
            num_particles_to_check = 5;
            for i = 1:length(obj.app.stack_paths)
                if contains(obj.app.stack_paths{i}, 'time_control') || contains(obj.app.stack_paths{i}, 'temp') || contains(obj.app.stack_paths{i}, 'cont')
                    continue;
                end
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.load_images_callback();
                if isfield(obj.app.stack_info, 'scale')
                    fprintf('Scales already exist for stack %s\n', obj.app.path);
                    if ~obj.app.forced
                        fprintf('Skipping stack %s\n', obj.app.path);
                        continue;
                    end
                end
                get_scale(num_particles_to_check);
                % average the obj.app.stack_info.radii
                radii_values = cellfun(@(x) x{2}, obj.app.stack_info.radii);
                avg_radius = mean(radii_values);
                std_radius = std(radii_values);
                obj.app.stack_info.scale = 2 * avg_radius;
                obj.app.stack_info.scale_std = std_radius;
                obj.app.utils.save_stack_callback();
            end
            function get_scale(n)
                % set the viewer to first image
                obj.app.utils.setFrame(1);
                fprintf('Getting scales for stack %s\n', obj.app.path);
                % get the particle locations
                particle_locations = obj.app.stack_info.particle_locations;
                % zoom onto 10 particles selected randomly from the image
                for i = 1:n
                    random_index = randi(size(particle_locations, 1));
                    % get the particle location
                    x = particle_locations.x(random_index);
                    y = particle_locations.y(random_index);
                    radius = particle_locations.size(random_index) * 3;
                    set(gcf, 'WindowScrollWheelFcn', {@obj.app.utils.zoom_callback,x,y});
                    % zoom onto the particle
                    xlim(obj.app.ui.controls.ax1, [x - 50, x + 50]);
                    ylim(obj.app.ui.controls.ax1, [y - 50, y + 50]);
        
                    % highlight the particle
                    hold(obj.app.ui.controls.ax1, 'on');
                    plot(obj.app.ui.controls.ax1, x, y, 'r*');
                    hold(obj.app.ui.controls.ax1, 'off');
                    % get the scale from the user by asking to draw a circle around the particle
                    h = drawcircle(obj.app.ui.controls.ax1, 'Center', [x, y], 'Radius', radius);
                    wait(h);
                    % save the circle to the stack_info
                    % add the scale to the stack_info
                    if ~isfield(obj.app.stack_info, 'radii')
                        obj.app.stack_info.radii = {};
                    end
                    obj.app.stack_info.radii{end+1} = {h.Center, h.Radius};
                end
                % delete the circle
                delete(h);
                % reset the zoom
                % zoom(obj.app.ui.controls.ax1,'reset');
                % enable the scroll callback
                set(gcf, 'WindowScrollWheelFcn', @obj.app.utils.scrollWheelMoved);
                % assignin('base', 'stack_info', obj.app.stack_info);
            end
        end
        function plot_scales(obj)
            scales = zeros(1, length(obj.app.stack_paths));
            sizes = zeros(1, length(obj.app.stack_paths));
            scales_std = zeros(1, length(obj.app.stack_paths));
            sizes_std = zeros(1, length(obj.app.stack_paths));
            % check if csv already exists
            if isfile('F:\shake_table_data\Results\scales.csv')
                data = readmatrix('F:\shake_table_data\Results\scales.csv');
                scales = data(:,1);
                scales_std = data(:,2);
                sizes = data(:,3);
                sizes_std = data(:,4);
            else
                for i = 1:length(obj.app.stack_paths)
                    if contains(obj.app.stack_paths{i}, 'time_control') || contains(obj.app.stack_paths{i}, 'temp') || contains(obj.app.stack_paths{i}, 'cont')
                        continue;
                    end
                    set(obj.app.ui.controls.stackDropdown, 'Value', i);
                    obj.app.load_images_callback();
                    if isfield(obj.app.stack_info, 'scale')                
                        scales(i) = obj.app.stack_info.scale;
                        scales_std(i) = obj.app.stack_info.scale_std;
                    end
                    sizes(i) = mean(obj.app.stack_info.particle_locations.size);
                    sizes_std(i) = std(obj.app.stack_info.particle_locations.size);
                end
            end
            cla(obj.app.ui.controls.ax2);

            errorbar(obj.app.ui.controls.ax2, scales,scales_std, 'b', 'LineWidth', 2);
            hold(obj.app.ui.controls.ax2, 'on');
            errorbar(obj.app.ui.controls.ax2, sizes * 3 * 2,sizes_std, 'g', 'LineWidth', 2);
            hold(obj.app.ui.controls.ax2, 'off');
            title('Scales and particle sizes');
            xlabel('Stack number');
            ylabel('Scale');
            % save data to csv
            data = [scales; scales_std; sizes; sizes_std];
            data = data';
            writematrix(data, 'F:\shake_table_data\Results\scales.csv');
            saveas(gcf, 'F:\shake_table_data\Results\scales.png');
        end
        %%%%%%%%%%%%%%%%%%%%%% TRIAL CHARACTERIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
        function is_late = is_rec_start_late(~, n,fs,iter)
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
            % WaitMessage = parfor_wait(length(obj.app.stack_paths), 'Waitbar', true);
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                % [N, fs] = obj.app.utils.get_info(path);
                % [iteration, ~] = obj.app.utils.getIteration(path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                obj.app.load_images_callback();
                obj.app.utils.setFrame(obj.app.stack_info.end_index - obj.app.stack_info.start_index + 1);
                % show dialog box to user to set the emptyornot field
                options = {'Empty', 'Not Empty'};
                answer = questdlg('Is the hopper empty?', 'Empty or Not', options{1}, options{2}, options{1});
                if strcmp(answer, options{1})
                    obj.app.stack_info.empty = true;
                else
                    obj.app.stack_info.empty = false;
                end
                obj.app.utils.save_stack_callback();
                % WaitMessage.Send;
            end
            % WaitMessage.Destroy;
        end
        function set_all_jam_or_not(obj)
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                [N, ~] = obj.app.utils.get_info(obj.app.path);
                % [iteration, ~] = obj.app.utils.getIteration(path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                obj.app.load_images_callback();
                obj.app.utils.setFrame(obj.app.stack_info.end_index - obj.app.stack_info.start_index + 1);
                if N == 4
                    obj.app.stack_info.jammed = true;
                else
                    % show dialog box to user to set the emptyornot field
                    options = {'Jammed', 'Not Jammed'};
                    answer = questdlg('Is the Jammed?', 'Jammed or Not', options{1}, options{2}, options{1});
                    if strcmp(answer, options{1})
                        obj.app.stack_info.jammed = true;
                    else
                        obj.app.stack_info.jammed = false;
                    end
                end
                obj.app.utils.save_stack_callback();
                % WaitMessage.Send;
            end
            % WaitMessage.Destroy;
        end
    end
end