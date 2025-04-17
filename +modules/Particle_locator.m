classdef Particle_locator < handle
    properties
        app
    end
    methods
        function obj = Particle_locator(app)
            obj.app = app;
            obj.addUi();
        end
        function addUi(obj)
            % Add UI elements for particle detection
            % obj.app.ui.controls.particleButton = uicontrol('Style', 'pushbutton', 'String', 'Detect Particles', ...
            %     'Position', [20, 20, 100, 30], 'Callback', @(src, event) obj.detect_particles_callback(src, event));
        end
        %%%%%%%%%%%%%%%%%%%%%% PARTICLE LOCATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function particle_locations = get_particle_locations(obj, image_idx)
            % check if obj.app.stack_info.particle_locations is a table
            if istable(obj.app.stack_info.particle_locations)
                fprintf('Particle locations are in old format\n');
                % means that its in old format
                % convert to cell array
                particle_locations = obj.app.stack_info.particle_locations;
                obj.app.stack_info.particle_locations = cell(1, numel(obj.app.stack_info.img_data.img_files));
                obj.app.stack_info.particle_locations{obj.app.stack_info.start_index} = particle_locations;

                obj.app.utils.save_stack_callback();
            end

            % check if particle_locations exist for index i
            if ~isfield(obj.app.stack_info, 'particle_locations') || isempty(obj.app.stack_info.particle_locations{image_idx})
                image_path = fullfile(obj.app.stack_info.img_data.img_files(image_idx).folder,...
                obj.app.stack_info.img_data.img_files(image_idx).name);
                % make a temp save path
                [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
                save_path = fullfile(parentDir, sprintf('particle_locations_%s.csv', iter));
                % get the particle locations from the image
                py.track.find_particle_locations(image_path=image_path, diam=int32(5),...
                max_iterations=int32(10), minmass=int32(1), separation=int32(5),...
                save_path=save_path);
                % load the saved csv from save_path
                particle_locations = readtable(save_path);
                obj.app.stack_info.particle_locations{image_idx} = particle_locations;
            else
                particle_locations = obj.app.stack_info.particle_locations{image_idx};
            end
            % clean up the particle locations
            % remove the locations that are in the mask
            if obj.app.stack_info.masked
                fprintf('Removing particle locations in the mask\n');
                mask = obj.app.stack_info.mask;
                % get the x and y coordinates of the particle locations
                x = particle_locations.x;
                y = particle_locations.y;
                % remove the locations that are in the mask
                rows_to_remove = [];
                for i = 1:length(x)
                    if mask(round(y(i)), round(x(i))) == 1
                        % fprintf('Removing particle location (%d, %d)\n',round(y(i)), round(x(i)));
                        rows_to_remove = [rows_to_remove, i];
                    end
                end
                % remove the rows from the particle locations
                particle_locations(rows_to_remove, :) = [];
            end

        end
        function draw_particle_locations(obj)
            idx = obj.app.current_image_idx;
            % get the particle locations from the image
            particle_locations = obj.get_particle_locations(idx);
            % fprintf('Drawing particle locations for image %d\n', idx);
            % imshow(image, 'Parent', obj.app.ui.controls.ax2);
            hold(obj.app.ui.controls.ax1, 'on');
            plot(obj.app.ui.controls.ax1, particle_locations.x, particle_locations.y, 'b*');
            hold(obj.app.ui.controls.ax1, 'off');
        end
        function toggle_particle_locations(obj, ~,~)
            if obj.app.particle_locations_visible
                % hide the particle locations
                obj.app.particle_locations_visible = false;
                % get slider index
                slider_idx = round(get(obj.app.ui.controls.slider, 'Value'));
                obj.app.utils.setFrame(slider_idx);
                set(gcf, 'WindowScrollWheelFcn', @obj.app.utils.scrollWheelMoved);
            else
                % show the particle locations
                if ~isfield(obj.app.stack_info, 'particle_locations')
                    obj.app.utils.display_warning('No particle locations found');
                    return;
                end
                obj.app.particle_locations_visible = true;
                obj.draw_particle_locations();
                set(gcf, 'WindowScrollWheelFcn', {});
            end
        end
    end
end