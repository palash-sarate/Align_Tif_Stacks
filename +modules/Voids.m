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
            if ~obj.app.forced && isfield(obj.app.stack_info, 'voids')
                % Check if voids already exist in the stack
                % msg = 'Voids already exist in this stack. Do you want to overwrite them?';
                % choice = questdlg(msg, 'Overwrite Voids', 'Yes', 'No', 'No');
                % if strcmp(choice, 'No')
                %     return;
                % end
                fprintf('Voids already exist in this stack. Skipping %s\n', obj.app.path);
                return;
            end            

            h2 = waitbar(0, 'Processing stack');
            % select equal spaced 100 frames to process from stack
            images_to_process = round(linspace(obj.app.stack_info.start_index, obj.app.stack_info.end_index, 100));
            % loop over all the frames to process
            for i = 1:length(images_to_process)
                % get the current frame index
                image_idx = images_to_process(i);                 
                % find the voids in the image
                obj.detect_voids(image_idx);
                % Update progress bar
                progress = i / length(images_to_process);
                waitbar(progress, h2);
            end
            % save the current stack info
            obj.app.utils.save_stack_callback();
            % Close the progress bar
            close(h2);
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
            % Create a progress bar
            h = waitbar(0, 'Detecting voids in image');
            % Find voids in the image using morphological operations
            particle_locations = obj.app.particle_locator.get_particle_locations(obj.app.current_image_idx);
            % find the white holes in the image
            bw_img = obj.get_locations_image(particle_locations);
            % find the voids in the image
            [B,L,N,A] = bwboundaries(bw_img, 'noholes');
            % add the B L N A to the app.stack_info.voids at the current index
            obj.append_voids(B, L, N, A, obj.app.current_image_idx);

            % only keep boundaries that don't touch the edge of the image
            B = obj.remove_holes_on_edge(B, size(bw_img, 1), size(bw_img, 2));
            % remove the holes that are too small
            B = B(cellfun(@(x) length(x) > 30, B));

            % plot the bw image on ax2
            cla(obj.app.ui.controls.ax2);
            % plot(obj.app.ui.controls.ax2, particle_locations.x, 800 - particle_locations.y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
            imshow(bw_img, 'Parent', obj.app.ui.controls.ax2);
            obj.overlay_boundaries(B);
            close(h);
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
        function B = remove_holes_on_edge(~, B, height, width)
            edge_boundaries = false(length(B), 1);
            for k = 1:length(B)
                boundary = B{k};
                if any(boundary(:,1) == 1) || any(boundary(:,1) == height) || ...
                   any(boundary(:,2) == 1) || any(boundary(:,2) == width)
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
        % area fraction of voids over time
        function analyze_all_stacks_voids(obj)
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
                % calculate the area fraction of voids in the current stack
                obj.analyze_stack_voids();
                % Update progress bar
                progress = i / length(obj.app.stack_paths);
                waitbar(progress, h);
            end
            close(h);
        end
        function analyze_stack_voids(obj)
            if ~isfield(obj.app.stack_info, 'voids')
                fprintf('No voids found in %s\n', obj.app.path);
                return;
            end
            h2 = waitbar(0, 'Processing stack');
            % load first image to get the size
            image_path = fullfile(obj.app.stack_info.img_data.img_files(1).folder, obj.app.stack_info.img_data.img_files(1).name);
            image = imread(image_path);
            img_height = size(image, 1);
            img_width = size(image, 2);
            total_area = sum(sum(~obj.app.stack_info.mask));

            % loop over all the frames to process
            for image_idx = 1:length(obj.app.stack_info.voids)
                % fprintf('Processing frame %d of %s\n', image_idx, obj.app.path);
                void_data = obj.app.stack_info.voids{image_idx};
                if isempty(void_data)
                    % fprintf('No voids found in frame %d of %s\n', i, obj.app.path);
                    continue;
                end
                % clean up the void data
                void_data = obj.clean_void_data(void_data, img_height, img_width);
                % calculate the area fraction of voids in the image                
                void_area = obj.calculate_void_area(void_data);
                % Radial distribution function of voids
                % save the void data to the app.stack_info.voids at the current index
                obj.app.stack_info.voids{image_idx}.void_area = void_area;
                obj.app.stack_info.voids{image_idx}.void_area_frac = void_area / total_area;

                % rdf = obj.get_rdf(void_data, img_height, img_width);
                % obj.app.stack_info.voids{image_idx}.rdf = rdf;

                % Update progress bar
                progress = image_idx / length(obj.app.stack_info.voids);
                waitbar(progress, h2);
            end
            % save the current stack info
            obj.app.utils.save_stack_callback();
            close(h2);
        end
        function total_void_area = calculate_void_area(~, void_data)
            % loop over the voids and calculate the total area of voids
            total_void_area = 0;
            for k = 1:length(void_data.B)
                boundary = void_data.B{k};
                % calculate the area of the void
                total_void_area = total_void_area + polyarea(boundary(:,2), boundary(:,1));
            end            
        end
        function rdf_data =  get_rdf(obj, void_data, img_height, img_width)
            % number of voids
            num_voids = length(void_data.B);
            if num_voids < 2
                rdf_data = [];
                return; % Not enough voids to calculate RDF
            end
        
            % Calculate centroids of each void
            centroids = zeros(num_voids, 2); % [x, y] coordinates
            for k = 1:num_voids
                boundary = void_data.B{k};
                centroids(k, :) = [mean(boundary(:, 2)), mean(boundary(:, 1))]; % [x, y]
            end
            centroids = array2table(centroids, 'VariableNames', {'x', 'y'});
            max_distance = sqrt(img_width^2 + img_height^2); % Maximum possible distance
            dr = 0.1; % Bin width for RDF calculation
            rdf_data = obj.app.rdf.calculate_gr(max_distance, dr, centroids);
            % Compute pairwise distances between void centroids
            % distances = pdist(centroids); % Pairwise distances as a vector

            % Define the bins for the RDF
            % num_bins = 50; % Number of bins
            % bin_edges = linspace(0, max_distance, num_bins + 1);
            % bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
            % Calculate the histogram of distances
            % counts = histcounts(distances, bin_edges);

            % Normalize the RDF
            % Calculate the area of each annular bin
            % bin_areas = pi * (bin_edges(2:end).^2 - bin_edges(1:end-1).^2);
            % Normalize by the total number of void pairs and the bin area
            % total_pairs = num_voids * (num_voids - 1) / 2; % Total number of pairs
            % rdf_values = counts ./ (bin_areas * total_pairs);

            % Return the RDF as a struct
            % rdf.bin_centers = bin_centers;
            % rdf.values = rdf_values;
        end
        % rdf of voids
        % density of voids over time
        % total area of voids over time
        function consolidate_voids_data(obj)
            all_data = struct();
            h = waitbar(0, 'Processing stacks');
            % loop over the stacks
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                [N, fs] = obj.app.utils.get_info(obj.app.path);
                [iteration, ~] = obj.app.utils.getIteration(obj.app.path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                obj.app.utils.load_stack_info();

                % get number of non empty voids in the stack
                num_voids = sum(~cellfun(@isempty, obj.app.stack_info.voids));

                void_area = zeros(num_voids, 1);
                void_area_frac = zeros(num_voids, 1);
                image_indexes = zeros(num_voids, 1);
                timeStamps = cell(num_voids, 1);

                % loop over voids
                temp_idx = 1;
                for j = 1:length(obj.app.stack_info.voids)
                    void_data = obj.app.stack_info.voids{j};
                    if isempty(void_data)
                        continue;
                    end
                    % clean up the void data
                    void_area(temp_idx) = void_data.void_area;
                    void_area_frac(temp_idx) = void_data.void_area_frac;
                    image_indexes(temp_idx) = j;
                    timeStamps{temp_idx} = obj.app.timer.predict_timeStamp(j);
                    temp_idx = temp_idx + 1;
                end
                voids_data = table(void_area, void_area_frac, image_indexes, timeStamps);
                all_data.(sprintf('N%d',N)).(sprintf('f%d',fs)).(sprintf('iter%s',iteration)) = voids_data;
                % Update progress bar
                progress = i / length(obj.app.stack_paths);
                waitbar(progress, h);
            end
            % save the voids data to results folder);
            save('F:\shake_table_data\Results\voids_data.mat', 'all_data');
            close(h);
        end
        function visualize_voids_data(obj)
            % load the voids data from results folder
            voids_data = load('F:\shake_table_data\Results\voids_data.mat');
            f = figure('Name', 'Area frac Visualization');
            ax = axes(f);
            hold(ax, 'on');
            f2 = figure('Name', 'Total area Visualization');
            ax2 = axes(f2);
            hold(ax2, 'on');
            f3 = figure('Name', 'Initial Area frac');
            ax3 = axes(f3);
            hold(ax3, 'on');
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                [N, fs] = obj.app.utils.get_info(obj.app.path);
                [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')                
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                empty = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration), '-mat', 'empty');
                if ~isfield(empty, 'empty')
                    obj.app.trial.set_stack_empty_or_not();
                    empty = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration), '-mat', 'empty');
                end
                if ~empty.empty                  
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                data = voids_data.all_data.(sprintf('N%d',N)).(sprintf('f%d',fs)).(sprintf('iter%s',iteration));
                plot(ax3, N, data.void_area_frac(1), 'o', 'Color',...
                    obj.app.utils.get_color(N), 'DisplayName', 'None');
                plot(ax, data.image_indexes, data.void_area_frac, 'Color',...
                    obj.app.utils.get_color(N), 'DisplayName', 'None');
                plot(ax2, data.image_indexes, data.void_area, 'Color',...
                    obj.app.utils.get_color(N), 'DisplayName', 'None');
            end

            unique_Ns = [4, 12, 24, 48];
            legend_handles = zeros(1, numel(unique_Ns));
            legend_handles2 = zeros(1, numel(unique_Ns));
            legend_handles3 = zeros(1, numel(unique_Ns));
            legend_entries = cell(1, numel(unique_Ns));
            
            for j = 1:numel(unique_Ns)
                N_val = unique_Ns(j);
                legend_handles(j) = plot(ax, NaN, NaN, 'LineWidth', 2, 'Color', obj.app.utils.get_color(N_val));
                legend_handles2(j) = plot(ax2, NaN, NaN, 'LineWidth', 2, 'Color', obj.app.utils.get_color(N_val));
                legend_handles3(j) = plot(ax3, NaN, NaN, 'LineWidth', 2, 'Color', obj.app.utils.get_color(N_val));
                legend_entries{j} = sprintf('N = %d', N_val);
            end
            
            legend(ax, legend_handles, legend_entries, 'Location', 'best');
            legend(ax2, legend_handles2, legend_entries, 'Location', 'best');
            legend(ax3, legend_handles3, legend_entries, 'Location', 'SouthEast');
            hold(ax, 'off');
            hold(ax2, 'off');
            hold(ax3, 'off');
            ax3.YLim = [0, 0.2];
            % set uninque N to be x axis ticks of ax3
            ax3.XTick = unique_Ns;
            ax3.XTickLabel = unique_Ns;
        end
        function create_images(obj)
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                [N, fs] = obj.app.utils.get_info(obj.app.path);
                [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')                
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                empty = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration), '-mat', 'empty');
                if ~isfield(empty, 'empty')
                    obj.app.trial.set_stack_empty_or_not();
                    empty = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration), '-mat', 'empty');
                end
                if ~empty.empty
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                data = voids_data.all_data.(sprintf('N%d',N)).(sprintf('f%d',fs)).(sprintf('iter%s',iteration));
            end
        end
        function data = clean_void_data(obj, data, height, width)
            % only keep boundaries that don't touch the edge of the image
            data.B = obj.remove_holes_on_edge(data.B, height, width);
            % remove the holes that are too small
            data.B = data.B(cellfun(@(x) length(x) > 30, data.B));
        end
    end
end