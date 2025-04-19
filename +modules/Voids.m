classdef Voids < handle
    properties
        app
    end
    methods
        function obj = Voids(app)
            obj.app = app;
        end
        function detect_voids_all_stacks(obj)
            % Create a progress bar
            h = waitbar(0, 'Processing stacks');
            % loop over all the stacks
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                % [N, fs] = obj.app.utils.get_info(path);
                % [iteration, ~] = obj.app.utils.getIteration(path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                obj.app.utils.load_stack_info();

                % detect voids in the current stack
                obj.detect_voids_stack();
                % Update progress bar
                progress = i / length(obj.app.stack_paths);
                waitbar(progress, h);
            end
            close(h);
        end
        function detect_voids_stack(obj)
            WaitMessage = parfor_wait(obj.app.stack_info.end_index - obj.app.stack_info.start_index + 1,...
             'Waitbar', true);
            % loop over all the frames
            for i = obj.app.stack_info.start_index:obj.app.stack_info.end_index
                % find the voids in the image
                obj.detect_voids(i);
                % Update progress bar
                WaitMessage.Send;
            end
            % save the current stack info
            obj.app.utils.save_stack_callback();
            % Close the progress bar
            WaitMessage.Destroy;
        end
        function detect_voids(obj, img_idx)
            % Find voids in the image using morphological operations
            particle_locations = obj.app.particle_locator.get_particle_locations(img_idx);
            % find the white holes in the image
            bw_img = obj.get_locations_image(particle_locations);
            % find the voids in the image
            [B,L,N,A] = bwboundaries(bw_img, 'noholes');
            % add the B L N A to the app.stack_info.voids at the current index
            obj.append_voids(B, L, N, A, img_idx);
        end
        function detect_voids_callback(obj)
            % Find voids in the image using morphological operations
            particle_locations = obj.app.particle_locator.get_particle_locations(obj.app.current_image_idx);
            % find the white holes in the image
            bw_img = obj.get_locations_image(particle_locations);
            % find the voids in the image
            [B,L,N,A] = bwboundaries(bw_img, 'noholes');
            % add the B L N A to the app.stack_info.voids at the current index
            obj.append_voids(B, L, N, A, obj.app.current_image_idx);

            % only keep boundaries that don't touch the edge of the image
            B = obj.remove_holes_on_edge(B, bw_img);
            % remove the holes that are too small
            B = B(cellfun(@(x) length(x) > 30, B));

            % plot the bw image on ax2
            cla(obj.app.ui.controls.ax2);
            % plot(obj.app.ui.controls.ax2, particle_locations.x, 800 - particle_locations.y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
            imshow(bw_img, 'Parent', obj.app.ui.controls.ax2);
            obj.overlay_boundaries(B);
        end
        function append_voids(obj, B, L, N, A, idx)
            % Append the voids to the app.stack_info.voids at the current index
            obj.app.stack_info.voids{idx}.B = B;
            obj.app.stack_info.voids{idx}.L = L;
            obj.app.stack_info.voids{idx}.N = N;
            obj.app.stack_info.voids{idx}.A = A;
        end
        function draw_boundaries(obj, B)
            colors=['b' 'g' 'r' 'c' 'm' 'y'];
            for k=1:length(B)
              boundary = B{k};
              cidx = mod(k,length(colors))+1;
              plot(obj.app.ui.controls.ax2, boundary(:,2), boundary(:,1),...
                   colors(cidx),'LineWidth', 2, 'HandleVisibility', 'off');
            
              %randomize text position for better visibility
            %   rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
            %   col = boundary(rndRow,2); row = boundary(rndRow,1);
            %   h = text(col+1, row-1, num2str(L(row,col)));
            %   set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
            end
        end
        function overlay_boundaries(obj, B)
            hold(obj.app.ui.controls.ax2, 'on');
            obj.draw_boundaries(B);
            hold(obj.app.ui.controls.ax2, 'off');
        end
        function B = remove_holes_on_edge(~, B, bw_img)
            img_height = size(bw_img, 1);
            img_width = size(bw_img, 2);
            edge_boundaries = false(length(B), 1);
            for k = 1:length(B)
                boundary = B{k};
                if any(boundary(:,1) == 1) || any(boundary(:,1) == img_height) || ...
                   any(boundary(:,2) == 1) || any(boundary(:,2) == img_width)
                    edge_boundaries(k) = true;
                end
            end
            B(edge_boundaries) = [];
        end
        function bw_img = get_locations_image(~, particle_locations)    
            % Define the size of the binary image (adjust as needed)
            img_height = 800; % Replace with the actual height of your image
            img_width = 800;  % Replace with the actual width of your image
            bw_img = ones(img_height, img_width); % Initialize a blank binary image
            % Define the radius of the circles
            radius = 5; % Adjust the radius as needed
            
            % Create a grid for the image
            [X, Y] = meshgrid(1:img_width, 1:img_height);
            
            % Iterate through particle locations and draw circles
            for i = 1:length(particle_locations.x)
                centerX = particle_locations.x(i);
                centerY = particle_locations.y(i);
                
                % Create a binary mask for the circle
                circleMask = (X - centerX).^2 + (Y - centerY).^2 <= radius^2;
                bw_img(circleMask) = 0; % Set circle pixels to 0 (black)
            end
        end
    end
end