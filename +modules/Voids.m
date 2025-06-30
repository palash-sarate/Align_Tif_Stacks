classdef Voids < handle
    properties
        app
        voids_data
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
            % clear the voids field in the stack_info
            obj.app.stack_info.voids = cell(length(obj.app.stack_info.img_data.img_files), 1);
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
                % Fill the hole with black
                % fill(obj.app.ui.controls.ax2, 800 - boundary(:,2), 800 - boundary(:,1), ...
                %     'k', 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(obj.app.ui.controls.ax2, boundary(:,2), boundary(:,1),...
                    colors(cidx), 'LineWidth', 2, 'HandleVisibility', 'off');
            
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
            bw_img = double(bw_img);
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
                particle_locations = obj.app.particle_locator.get_particle_locations(image_idx);
                bw_image = obj.get_locations_image(particle_locations);
                % clean up the void data
                void_data = obj.clean_void_data(void_data, img_height, img_width);
                % calculate the area fraction of voids in the image                
                void_area = obj.calculate_void_area(void_data);
                chain_area = sum(bw_image(:) == 0);
                % Radial distribution function of voids
                % save the void data to the app.stack_info.voids at the current index
                obj.app.stack_info.voids{image_idx}.void_area = void_area;
                obj.app.stack_info.voids{image_idx}.void_area_frac = void_area / total_area;
                obj.app.stack_info.voids{image_idx}.chain_area = chain_area;
                obj.app.stack_info.voids{image_idx}.chain_area_frac = chain_area / total_area;
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
                voids = obj.app.stack_info.voids;
                % get number of non empty voids in the stack
                num_voids = sum(~cellfun(@isempty, voids));

                void_area = zeros(num_voids, 1);
                void_area_frac = zeros(num_voids, 1);
                image_indexes = zeros(num_voids, 1);
                timeStamps = cell(num_voids, 1);

                % loop over voids
                temp_idx = 1;
                for j = 1:length(voids)
                    void_data = voids{j};
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
                voids_table = table(void_area, void_area_frac, image_indexes, timeStamps);
                all_data.(sprintf('N%d',N)).(sprintf('f%d',fs)).(sprintf('iter%s',iteration)) = voids_table;
                % Update progress bar
                progress = i / length(obj.app.stack_paths);
                waitbar(progress, h);
            end
            % save the voids data to results folder);
            save('F:\shake_table_data\Results\voids_data.mat', 'all_data');
            close(h);
        end
        function visualize_voids_characteristics(obj)
            % load the voids data from results folder
            obj.voids_data = load('F:\shake_table_data\Results\voids_data.mat');
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
                data = obj.voids_data.all_data.(sprintf('N%d',N)).(sprintf('f%d',fs)).(sprintf('iter%s',iteration));
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
            % save the figures
            save_dir = "F://shake_table_data//Results//voids_images//";
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            exportgraphics(ax, sprintf('%s//voids_area_frac.png', save_dir));
            exportgraphics(ax2, sprintf('%s//voids_area.png', save_dir));
            exportgraphics(ax3, sprintf('%s//voids_initial_area_frac.png', save_dir));
        end
        function plot_void_area_fraction_over_time_stacks(obj)
            % Load the voids data from results folder
            obj.voids_data = load('F:\shake_table_data\Results\voids_data.mat');
            unique_Ns = [4, 12, 24, 48];
            
            % Loop through each unique N to create separate plots
            for j = 1:numel(unique_Ns)
                N_val = unique_Ns(j);
                
                % Create a new figure for each N
                f = figure('Name', sprintf('Area frac Visualization for N = %d', N_val));
                ax = axes(f);
                hold(ax, 'on');
                
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
                        fprintf('Skipping as empty %s\n', obj.app.path);
                        continue;
                    end
                    
                    if N == N_val
                        data = obj.voids_data.all_data.(sprintf('N%d', N)).(sprintf('f%d', fs)).(sprintf('iter%s', iteration));
                        normalized_x = ((data.image_indexes-min(data.image_indexes)) / max(data.image_indexes)) * 100;
                        plot(ax, normalized_x, data.void_area_frac, 'Color', ...
                            obj.app.utils.get_color(fs));
                    end
                end
                
                % Add legend and labels for unique fs
                % Get all unique fs values for this N
                fs_list = [4,6,8,10,12,14,16,18,20];
                unique_fs = unique(fs_list);

                % Create legend handles with correct colors for each fs
                legend_handles = gobjects(1, numel(unique_fs));
                legend_entries = cell(1, numel(unique_fs));
                for k = 1:numel(unique_fs)
                    legend_handles(k) = plot(ax, NaN, NaN, 'Color', obj.app.utils.get_color(0,unique_fs(k),max(unique_fs)), 'LineWidth', 2);
                    legend_entries{k} = sprintf('fs = %d', unique_fs(k));
                end
                legend(ax, legend_handles, legend_entries, 'Location', 'best');
                xlabel(ax, 'Percent Completion (%)');
                ylabel(ax, 'Area Fraction');
                title(ax, sprintf('Area Fraction Visualization for N = %d', N_val));
                hold(ax, 'off');
                
                % Save the figure
                save_dir = "F://shake_table_data//Results//voids_images//";
                if ~exist(save_dir, 'dir')
                    mkdir(save_dir);
                end
                exportgraphics(ax, sprintf('%s//voids_area_frac_N%d.png', save_dir, N_val));
            end
        end
        function create_void_images(obj)
            h = waitbar(0, 'Processing stacks');
            for i = 72%30:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                % [N, fs] = obj.app.utils.get_info(obj.app.path);
                [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')                
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                obj.app.utils.load_stack_info();
                start_frame_timestamp = obj.app.timer.time_2_sec(obj.app.stack_info.timestamps{obj.app.stack_info.start_index});

                for j = 24821%1:10:length(obj.app.stack_info.voids)
                    voids = obj.app.stack_info.voids{j};
                    if isempty(voids)
                        continue;
                    end
                    fprintf("processing voids for image %d\n", j);
                    particle_locations = obj.app.particle_locator.get_particle_locations(j);
                    % find the white holes in the image
                    bw_img = obj.get_locations_image(particle_locations);

                    %%%%%%%%%%%%%%%%% find the voids in the image
                    % [B,~,~,~] = bwboundaries(bw_img, 'noholes');
                    % % add the B L N A to the app.stack_info.voids at the current index
                    % % obj.append_voids(B, L, N, A, obj.app.current_image_idx);

                    % % only keep boundaries that don't touch the edge of the image
                    % B = obj.remove_holes_on_edge(B, size(bw_img, 1), size(bw_img, 2));
                    % % remove the holes that are too small
                    % B = B(cellfun(@(x) length(x) > 30, B));
                    %%%%%%%%%%%%%%%%%
                    
                    % find the voids in the image
                    clean_voids = obj.clean_void_data(voids, size(bw_img, 1), size(bw_img, 2));

                    % plot the bw image on ax2
                    cla(obj.app.ui.controls.ax2);
                    imshow(bw_img, 'Parent', obj.app.ui.controls.ax2);
                    plot(obj.app.ui.controls.ax2, particle_locations.x, 800 - particle_locations.y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
                    obj.overlay_boundaries(clean_voids.B);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % % add the original image to the ax2 to the bottom right corner
                    % image = imread(fullfile(obj.app.stack_info.img_data.img_files(j).folder, obj.app.stack_info.img_data.img_files(j).name));
                    % resizedGray = imresize(image, [overlay_height, overlay_width]);
                    % startRow = height - overlay_height + 1;
                    % startCol = width - overlay_width + 1;
                    % endRow = height;
                    % endCol = width;
                    % bw_img(startRow:endRow, startCol:endCol) = resizedGray;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % add timestamp to the image
                    timeStamp = obj.app.timer.predict_timeStamp(j);
                    secs = obj.app.timer.time_2_sec(timeStamp);
                    text(obj.app.ui.controls.ax2, 'Units', 'normalized', 'Position', [0.99, 0.03], ...
                        'String', sprintf("%.2f Sec", secs-start_frame_timestamp), ...
                        'Color', 'black', 'FontSize', 18, 'HorizontalAlignment', 'right');
                    % save the image
                    % remove axes
                    obj.app.ui.controls.ax2.XColor = 'none';
                    obj.app.ui.controls.ax2.YColor = 'none';
                    % Set the axis to 800x800 pixels before saving
                    set(obj.app.ui.controls.ax2, 'Units', 'pixels', 'Position', [1 1 800 800]);
                    set(gcf, 'Units', 'pixels', 'Position', [100 100 800 800]);

                    save_dir = sprintf("%s//voids_images_%s//", parentDir, iteration);
                    if ~exist(save_dir, 'dir')
                        mkdir(save_dir);
                    end                    
                    exportgraphics(obj.app.ui.controls.ax2, sprintf('%s//voids_%d.png', save_dir, j));
                end
                progress = i / length(obj.app.stack_paths);
                waitbar(progress, h);
            end
            close(h);
        end
        function data = clean_void_data(obj, data, height, width)
            % only keep boundaries that don't touch the edge of the image
            data.B = obj.remove_holes_on_edge(data.B, height, width);
            % remove the holes that are too small
            data.B = data.B(cellfun(@(x) length(x) > 30, data.B));
        end
        function plot_size_distribution(obj)
            f = figure('Name', 'Voids Size Distribution');
            ax = axes(f);
            hold(ax, 'on');
            % if hist_data already exists, load it
            if evalin('base', 'exist(''hist_data'', ''var'')')
                hist_data = evalin('base', 'hist_data');
            else
                hist_data = struct();
                for i = 1:length(obj.app.stack_paths)
                    set(obj.app.ui.controls.stackDropdown, 'Value', i);
                    obj.app.path = obj.app.stack_paths{i};
                    [N, fs] = obj.app.utils.get_info(obj.app.path);
                    [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
                    if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')                
                        fprintf('Skipping %s\n', obj.app.path);           
                        continue;
                    end
                    % obj.app.utils.load_stack_info();
                    voids_struct = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration), '-mat', 'voids');
                    % assignin('base', 'voids_struct', voids_struct);
                    if ~isfield(voids_struct, 'voids')
                        fprintf('No voids found in %s\n', obj.app.path);
                        continue;
                    end
                    voids = voids_struct.voids;
                    for j = 1:length(voids)
                        void_data = voids{j};
                        if isempty(void_data)
                            continue;
                        end
                        % calculate the area fraction of voids in the current stack
                        if ~isfield(hist_data, sprintf('N%d', N))
                            hist_data.(sprintf('N%d', N)) = [];
                        end
                        boundaries = void_data.B;
                        % only keep boundaries that don't touch the edge of the image
                        boundaries = obj.remove_holes_on_edge(boundaries, 800, 800);
                        % remove the holes that are too small
                        boundaries = boundaries(cellfun(@(x) length(x) > 30, boundaries));
                        areas = zeros(length(boundaries), 1);
                        for k = 1:length(boundaries)
                            boundary = boundaries{k};
                            % calculate the area of the void
                            areas(k) = polyarea(boundary(:,2), boundary(:,1));
                        end
                        hist_data.(sprintf('N%d', N)) = [hist_data.(sprintf('N%d', N)); areas];
                        break;
                    end
                    fprintf('collected voids in %s\n', obj.app.path);
                end
                assignin('base', 'hist_data', hist_data);
                % save the voids data to results folder);
                save('F:\shake_table_data\Results\voids_hist_data.mat', 'hist_data');
            end
            % loop over the hist_data and plot the histogram of voids for each N
            unique_Ns = fieldnames(hist_data);
            for i = 1:length(unique_Ns)
                N = str2double(unique_Ns{i}(2:end));
                areas = hist_data.(sprintf('N%d', N));
                areas = areas(areas < 1500);
                % plot the histogram of voids for each N
                histogram(ax, areas, 100, 'Normalization', 'count', 'DisplayName', sprintf("N = %d",N), 'EdgeColor', obj.app.utils.get_color(N), 'FaceColor', obj.app.utils.get_color(N),'FaceAlpha',0.1);
            end
            legend(ax, 'Location', 'best');
            hold(ax, 'off');
        end
        function plot_void_size_distr_over_time(obj)
            % if hist_data already exists, load it
            % hist_data = struct();
            for i = 60%1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                % [N, fs] = obj.app.utils.get_info(obj.app.path);
                [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')                
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                % obj.app.utils.load_stack_info();
                voids_struct = load(sprintf('%s//stack_info_%s.mat', parentDir, iteration), '-mat', 'voids');
                % assignin('base', 'voids_struct', voids_struct);
                if ~isfield(voids_struct, 'voids')
                    fprintf('No voids found in %s\n', obj.app.path);
                    continue;
                end
                f = figure('Name', 'Void area Distribution');
                ax = axes(f);
                hold(ax, 'on');
                voids = voids_struct.voids;
                % remove empty voids
                voids = voids(~cellfun(@isempty, voids));
                % get length of non empty voids
                num_voids = length(voids);
                % make list of colors from gray colormap 
                colormap = summer(num_voids);
                colors = colormap(1:num_voids, :);
                for j = 1:num_voids
                    void_data = voids{j};
                    if isempty(void_data)
                        continue;
                    end
                    % % calculate the areas of voids in the current stack
                    % if ~isfield(hist_data, sprintf('N%d', N))
                    %     hist_data.(sprintf('N%d', N)) = [];
                    % end
                    boundaries = void_data.B;
                    areas = zeros(length(boundaries), 1);
                    for k = 1:length(boundaries)
                        boundary = boundaries{k};
                        % calculate the area of the void
                        areas(k) = polyarea(boundary(:,2), boundary(:,1));
                    end
                    areas = areas(areas < 1000);
                    % create a line plot of the histogram of voids
                    % [N,edges] = histcounts(areas);
                    % edges = edges(2:end) - (edges(2)-edges(1))/2;
                    % plot(edges, N, 'Color', 'b', 'DisplayName', 'None');
                    histogram(ax, areas, 200, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', colors(j,:));
                    % break;
                end
                hold(ax, 'off');
            end
        end
        function [anisotropy, eigenVectors, eigenValues] = get_eigen_vectors(obj, image_idx)
            % Assuming B, L, N, A are already computed using:
            % [B,L,N,A] = bwboundaries(binaryImage, 'noholes');
            % image_idx = obj.app.current_image_idx;
            particle_locations = obj.app.particle_locator.get_particle_locations(image_idx);
            binaryImage = obj.get_locations_image(particle_locations);
            data = obj.app.stack_info.voids{image_idx};
            data = obj.clean_void_data(data, size(binaryImage, 1), size(binaryImage, 2));
            B = data.B;
            N = length(B);
            % Initialize arrays to store results
            anisotropy = zeros(N, 1);
            eigenVectors = cell(N, 1);
            eigenValues = cell(N, 1);

            for k = 1:N
                % Extract boundary coordinates for void k
                boundary = B{k};
                x = boundary(:, 2); % column indices
                y = boundary(:, 1); % row indices
                
                % Center the data
                x_mean = sum(x)/polyarea(x, y);
                y_mean = sum(y)/polyarea(x, y);
                x_centered = x - x_mean;
                y_centered = y - y_mean;
                
                % Create a 2xN matrix of the centered coordinates
                coords = [x_centered, y_centered];
                
                % Compute the covariance matrix
                C = cov(coords);
                
                % Compute eigenvalues and eigenvectors
                [V, D] = eig(C);  % V: eigenvectors, D: diagonal eigenvalue matrix
                
                % Sort eigenvalues and vectors in descending order
                [eigvals, idx] = sort(diag(D), 'descend');
                V = V(:, idx); % sort eigenvectors accordingly
                D = diag(eigvals); % sorted eigenvalues
                
                % Store eigenvectors and eigenvalues
                eigenVectors{k} = V;
                eigenValues{k} = D;
                
                % Compute shape anisotropy
                lambda1 = eigvals(1);
                lambda2 = eigvals(2);
                anisotropy(k) = 1 - (lambda2 / lambda1);
            end
            obj.app.stack_info.voids{image_idx}.anisotropy = anisotropy;
            obj.app.stack_info.voids{image_idx}.eigenVectors = eigenVectors;
            % Display summary for each pore
            % for k = 1:N
            %     fprintf('Void %d:\n', k);
            %     fprintf('  Eigenvalues: %.3f, %.3f\n', eigenValues{k}(1,1), eigenValues{k}(2,2));
            %     fprintf('  Anisotropy: %.3f\n', anisotropy(k));
            %     fprintf('  Major axis direction: [%.3f %.3f]\n\n', eigenVectors{k}(:,1));
            % end
        end
        function overlay_eigen_vectors(obj)
            image_idx = obj.app.current_image_idx;
            particle_locations = obj.app.particle_locator.get_particle_locations(image_idx);
            binaryImage = obj.get_locations_image(particle_locations);
            data = obj.app.stack_info.voids{image_idx};
            data = obj.clean_void_data(data, size(binaryImage, 1), size(binaryImage, 2));
            B = data.B;
            N = length(B);
            cla(obj.app.ui.controls.ax2);
            % Show the binary image
            imshow(binaryImage, 'Parent', obj.app.ui.controls.ax2);
            obj.overlay_boundaries(B);
            hold on;

            % Loop through each void
            for k = 1:N
                boundary = B{k};
                x = boundary(:, 2); % columns
                y = boundary(:, 1); % rows

                % Compute centroid
                x_mean = mean(x);
                y_mean = mean(y);

                % Centered coordinates
                x_centered = x - x_mean;
                y_centered = y - y_mean;
                coords = [x_centered, y_centered];

                % Covariance and eigen decomposition
                C = cov(coords);
                [V, D] = eig(C);
                
                % Sort eigenvalues and vectors
                [~, idx] = sort(diag(D), 'descend');
                V = V(:, idx);
                
                % Scale eigenvectors for visualization
                scale = 10; % adjust this for visibility
                v1 = V(:,1) * scale; % major axis
                v2 = V(:,2) * scale; % minor axis

                % Plot major and minor axis as arrows
                quiver(x_mean, y_mean, v1(1), v1(2), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 2,'DisplayName', 'None');
                quiver(x_mean, y_mean, -v1(1), -v1(2), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 2,'DisplayName', 'None'); % opposite direction

                quiver(x_mean, y_mean, v2(1), v2(2), 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 2,'DisplayName', 'None');
                quiver(x_mean, y_mean, -v2(1), -v2(2), 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 2,'DisplayName', 'None'); % opposite direction

                % Optionally, mark the centroid
                plot(x_mean, y_mean, 'go', 'MarkerSize', 5, 'LineWidth', 2,'DisplayName', 'None');
            end

            hold off;
            title('Eigenvectors (red: major axis, blue: minor axis)');
            % save the ax2 axis as voids.eps file
            exportgraphics(obj.app.ui.controls.ax2, 'F://shake_table_data//Sauiba//voids_eigenvectors.eps', 'ContentType', 'vector');
        end
        function analyze_2d_fabric_frame(obj, image_idx)
            % Input: majorEigenVectors - Nx2 matrix of 2D major axis vectors for each void
            % Each row should be a unit vector: [vx, vy]
            [~, eigenVectors, ~] = obj.get_eigen_vectors(image_idx);
            N = length(eigenVectors);
            majorEigenVectors = zeros(N, 2);

            for i = 1:N
                majorEigenVectors(i, :) = eigenVectors{i}(:,1)';  % major axis eigenvector
            end

            % Normalize the vectors
            norms = sqrt(sum(majorEigenVectors.^2, 2));
            majorEigenVectors = majorEigenVectors ./ norms;

            % --- Fabric tensor computation ---
            fabric = zeros(2);
            for i = 1:size(majorEigenVectors, 1)
                v = majorEigenVectors(i, :)';
                fabric = fabric + (v * v');  % outer product dyadic product
            end
            fabric = fabric / size(majorEigenVectors, 1); % average

            % --- Deviatoric fabric ---
            trace_fabric = trace(fabric);
            hydro = trace_fabric / 2;
            dev_fabric = fabric - hydro * eye(2);
            f_dev = dev_fabric * (3/2);  % scaling factor (similar to 3D: 15/2)

            % --- Anisotropy norm (Frobenius norm of deviatoric tensor) ---
            norm_dev = sqrt((3/2) * sum(sum(f_dev .* f_dev)));
            % fprintf('Anisotropy norm (2D): %.4f\n', norm_dev);

            % --- Orientation distribution plot (polar plot) ---
            theta = linspace(0, pi, 180);  % angles in radians (0 to π for axis-aligned symmetry)
            radius = zeros(size(theta));
            for i = 1:length(theta)
                v = [cos(theta(i)); sin(theta(i))];
                radius(i) = (1 / pi) * (1 + v' * f_dev * v);  % 2D analogue of ODF
            end

            % --- Polar plot ---
            % figure;
            % polarplot(obj.app.ui.controls.ax2, [theta, theta + pi], [radius, radius], 'r-', 'LineWidth', 2);  % symmetric about π
            % title('2D Orientation Distribution (Major Axis)');
            % ax = gca;
            % ax.ThetaZeroLocation = 'top';
            % ax.ThetaDir = 'clockwise';
            % add to stack info
            obj.app.stack_info.voids{image_idx}.fabric = fabric;
            obj.app.stack_info.voids{image_idx}.dev_fabric = dev_fabric;
            obj.app.stack_info.voids{image_idx}.f_dev = f_dev;
            obj.app.stack_info.voids{image_idx}.norm_dev = norm_dev;
            obj.app.stack_info.voids{image_idx}.radius = radius;
            obj.app.stack_info.voids{image_idx}.theta = theta;
            % fprintf('Processed image %d/%d\n', image_idx, length(obj.app.stack_info.voids));
        end
        function analyze_2d_fabric_stack(obj)
            % Create a progress bar
            h = waitbar(0, 'Processing stack for 2D fabric analysis');
            n = length(obj.app.stack_info.voids);
            voids = obj.app.stack_info.voids;
            
            for i = 1:n
                % Check if the voids data is empty
                if isempty(voids{i})
                    progress = i / n;
                    waitbar(progress, h);
                    continue; % Skip this iteration if voids data is empty
                end
                % Call the analyze_2d_fabric_frame function for each image
                obj.analyze_2d_fabric_frame(i);
                % Update progress bar
                progress = i / n;
                waitbar(progress, h);
            end
            % save the stack info
            obj.app.utils.save_stack_callback();
            close(h);
        end
        function analyze_2d_fabric_stacks(obj)
            % Create a progress bar
            h = waitbar(0, 'Processing stacks for 2D fabric analysis');
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

                obj.analyze_2d_fabric_stack();
                % Update progress bar
                progress = i / length(obj.app.stack_paths);
                waitbar(progress, h);
            end
            close(h);
        end
        function plot_anisotropy_and_angle_evolution_stack(obj)
            % Load the pooled voids data from results folder
            if isempty(obj.voids_data)
                obj.voids_data = load('F:\shake_table_data\Results\voids_data.mat');
            end
            [N, fs] = obj.app.utils.get_info(obj.app.path);
            iteration = obj.app.stack_info.iteration;
            parentDir = obj.app.stack_info.parentDir;
            % Ensure save directory exists
            save_dir = fullfile(parentDir,sprintf("anisotropy_angle_frames_%s",iteration));
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            void_orientation_save_dir = fullfile(parentDir,sprintf("void_orientation_%s",iteration));
            if ~exist(void_orientation_save_dir, 'dir')
                mkdir(void_orientation_save_dir);
            end
            void_anisotropy_hist_save_dir = fullfile(parentDir,sprintf("anisotropy_hist_%s",iteration));
            if ~exist(void_anisotropy_hist_save_dir, 'dir')
                mkdir(void_anisotropy_hist_save_dir);
            end
            n_frames = length(obj.app.stack_info.voids);
            voids = obj.app.stack_info.voids;
            fprintf('Processing %d frames for anisotropy_and_angle_evolution\n', n_frames);
            total_area = sum(sum(~obj.app.stack_info.mask));
            % check if the non empty voids data have anisotropy, radius and theta fields
            % Find the first and last non-empty voids index
            first_valid = [];
            first_valid_index = [];
            % last_valid = [];
            last_valid_index = [];
            for idx = 1:length(voids)
                if ~isempty(voids{idx})
                    if isempty(first_valid)
                        first_valid = voids{idx};
                        first_valid_index = idx;
                    end
                    % last_valid = voids{idx};
                    last_valid_index = idx;
                end
            end
            if ~isempty(first_valid) && isfield(first_valid, 'anisotropy') && isfield(first_valid, 'theta') && isfield(first_valid, 'radius')
                % All fields exist
                disp('First non-empty voids data has all required fields.');
            else
                disp('First non-empty voids data is missing one or more required fields.');
                disp('Analyzing 2D fabric stack...');
                obj.analyze_2d_fabric_stack();
                % Re-check after analysis
                voids = obj.app.stack_info.voids;
                first_valid = voids{first_valid_index};
                if ~isfield(first_valid, 'anisotropy') || ~isfield(first_valid, 'theta') || ~isfield(first_valid, 'radius')
                    error('Required fields are still missing after analysis.');
                end
                disp('All required fields are now present.');
            end

            anisotropy_norm = [];anisotropy_norm_frame = [];
            start_frame_timestamp = obj.app.timer.time_2_sec(obj.app.stack_info.timestamps{obj.app.stack_info.start_index});
            data = obj.voids_data.all_data.(sprintf('N%d', N)).(sprintf('f%d', fs)).(sprintf('iter%s', iteration));
            normalized_x = ((data.image_indexes-min(data.image_indexes)) / max(data.image_indexes)) * 100;
            chain_area_frac = zeros(1, numel(data.image_indexes));
                        
            h1 = waitbar(0, 'Processing stack for anisotropy and angle evolution');
            for i = 1:n_frames
                void_data = voids{i};
                if isempty(void_data) || ~isfield(void_data, 'anisotropy') || ~isfield(void_data, 'theta') || ~isfield(void_data, 'radius')
                    continue;
                end
                % fprintf('Processing frame %d/%d\n', i, n_frames);
                void_data = obj.clean_void_data(void_data, 800, 800);
                timeStamp = obj.app.stack_info.timestamps{i};
                secs = obj.app.timer.time_2_sec(timeStamp);
                timestamp_text = sprintf("%.0f Sec", secs-start_frame_timestamp);
                particle_locations = obj.app.particle_locator.get_particle_locations(i);
                bw_image = obj.get_locations_image(particle_locations);
                % Prepare figure
                f = figure('Visible', 'off', 'Position', [100, 100, 1000, 800]);
                % add figure title
                sgtitle(sprintf('Anisotropy and Orientation distribution (fs = %d, N = %d, t = %s)', fs, N, timestamp_text), 'FontSize', 16);
                t = tiledlayout(f, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

                % --- Anisotropy histogram ---
                ax1 = nexttile(t, 1); 
                % Row 1, Col 1
                if isfield(void_data, 'anisotropy')
                    % Compute histogram data
                    [counts, edges] = histcounts(void_data.anisotropy, 20, 'Normalization', 'pdf');
                    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
                    plot(ax1, bin_centers, counts, 'Color', [0.2 0.6 0.8], 'LineWidth', 2);
                    xlabel(ax1, 'Anisotropy');
                    ylabel(ax1, 'PDF');
                    title(ax1, 'Anisotropy Distribution');
                    xlim(ax1, [0 1]);
                    exportgraphics(ax1, sprintf('%s//anisotropy_hist_%d.png', void_anisotropy_hist_save_dir, i));
                end

                % --- Void area fraction ---
                ax2 = nexttile(t, 9); % Row 3, Col 1
                end_index = find(data.image_indexes == i, 1, 'last');
                plot(ax2, normalized_x(1:end_index), data.void_area_frac(1:end_index), 'Color', [0.2 0.6 0.8], 'LineWidth', 2);
                ax2.XLim = [0 100];
                title(ax2, 'Void Area Fraction Evolution');
                ax2.XLabel.String = 'Percent Completion (%)';
                ax2.YLabel.String = 'Void Area Fraction';

                % --- Chains area fraction (big subplot, col 2) ---
                ax3 = nexttile(t, 5); % Span 3 rows, col 2
                end_index = find(data.image_indexes <= i, 1, 'last');
                chain_area = sum(bw_image(:) == 0);%1-white,0-black
                chain_area_frac(end_index) = chain_area / total_area; % Calculate area fraction
                plot(ax3, normalized_x(1:end_index), chain_area_frac(1:end_index), 'Color', [0.2 0.6 0.8], 'LineWidth', 2);
                ax3.XLim = [0 100];
                title(ax3, 'Chains Area Fraction Evolution');
                ax3.XLabel.String = 'Percent Completion (%)';
                ax3.YLabel.String = 'Chains Area Fraction';

                % --- orientation distribution (polar hist) ---
                nexttile(t, 2, [2 1]);
                eigenVectors = void_data.eigenVectors;
                if ~isempty(eigenVectors)
                    num_vec = length(eigenVectors);
                    majorEigenVectors = zeros(num_vec, 2);
                    for ii = 1:num_vec
                        majorEigenVectors(ii, :) = eigenVectors{ii}(:,1)';
                    end
                    cv = majorEigenVectors;
                    cv1 = cv;
                    cv2 = -cv;
                    cv3 = [cv1; cv2];
                    theta = atan2(cv3(:,2), cv3(:,1));
                    theta(theta < 0) = theta(theta < 0) + 2*pi;
                    polarhistogram(theta, 18);
                    title('Orientation Distribution');
                    pax = gca; % Get the polar axes
                    pax.ThetaZeroLocation = 'right';
                    pax.ThetaDir = 'counterclockwise';
                    pax.RAxis.Label.String = 'Frequency';
                    exportgraphics(pax, sprintf('%s//void_orientation_%d.png', void_orientation_save_dir, i));
                end
                % --- Scalar Anisotropy plot ---
                ax5 = nexttile(t, 10, [1 1]);
                if isfield(void_data, 'anisotropy')
                    anisotropy_norm_frame = [anisotropy_norm_frame; (i-obj.app.stack_info.start_index)/(obj.app.stack_info.end_index-obj.app.stack_info.start_index)];
                    anisotropy_norm = [anisotropy_norm; void_data.norm_dev];
                    plot(anisotropy_norm_frame * 100, anisotropy_norm, 'Color', [0.2 0.6 0.8], 'LineWidth', 2);
                    xlabel('Percent Completion (%)');
                    ylabel('Anisotropy Norm');
                    title('Anisotropy Norm Evolution');
                    xlim([1 100]);
                else
                    text(0.5, 0.5, 'No Anisotropy Data', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end
                % --- Voids boundaries with eigen vectors (big subplot, last 2 cols) ---
                ax6 = nexttile(t, 3, [3 2]); % Start at row 1, col 3, span 3 rows and 2 columns
                axes(ax6); hold on;
                voids_save_dir = sprintf("%s//voids_images_%s//", parentDir, iteration);
                if ~exist(voids_save_dir, 'dir')
                    mkdir(voids_save_dir);
                end 
                img_filename = sprintf('%s//voids_%d.png', voids_save_dir, i);
                if exist(img_filename, 'file')
                    imshow(img_filename, 'Parent', ax6);
                else
                    imshow(bw_image);
                    B = void_data.B;
                    colors = ['b' 'g' 'r' 'c' 'm' 'y'];
                    B_len = length(B);
                    centroids = zeros(B_len,2);
                    V1 = zeros(B_len,2);
                    V2 = zeros(B_len,2);

                    % Parallel computation of centroids and eigenvectors
                    parfor k = 1:B_len
                        boundary = B{k};
                        x = boundary(:,2);
                        y = boundary(:,1);
                        x_mean = mean(x);
                        y_mean = mean(y);
                        coords = [x - x_mean, y - y_mean];
                        C = cov(coords);
                        [V, D] = eig(C);
                        [~, idx_sort] = sort(diag(D), 'descend');
                        V = V(:, idx_sort);
                        centroids(k,:) = [x_mean, y_mean];
                        V1(k,:) = V(:,1)';
                        V2(k,:) = V(:,2)';
                    end

                    hold on
                    for k = 1:B_len
                        boundary = B{k};
                        cidx = mod(k-1, length(colors)) + 1;
                        plot(boundary(:,2), boundary(:,1), colors(cidx), 'LineWidth', 2, 'HandleVisibility', 'off');
                        scale = 10;
                        x_mean = centroids(k,1);
                        y_mean = centroids(k,2);
                        v1 = V1(k,:) * scale;
                        v2 = V2(k,:) * scale;
                        quiver(x_mean, y_mean, v1(1), v1(2), 0, 'r', 'LineWidth', 1, 'MaxHeadSize', 1, 'HandleVisibility', 'off');
                        quiver(x_mean, y_mean, -v1(1), -v1(2), 0, 'r', 'LineWidth', 1, 'MaxHeadSize', 1, 'HandleVisibility', 'off');
                        quiver(x_mean, y_mean, v2(1), v2(2), 0, 'b', 'LineWidth', 1, 'MaxHeadSize', 1, 'HandleVisibility', 'off');
                        quiver(x_mean, y_mean, -v2(1), -v2(2), 0, 'b', 'LineWidth', 1, 'MaxHeadSize', 1, 'HandleVisibility', 'off');
                        plot(x_mean, y_mean, 'ro', 'MarkerSize', 0.5, 'LineWidth', 0.2, 'HandleVisibility', 'off');
                    end
                    hold off
                    axis(ax6, 'image');
                    exportgraphics(ax6, sprintf('%s//voids_%d.png', voids_save_dir, i));
                end
                % save all the axes
                if(i==last_valid_index)
                    % Save the anisotropy histogram as an image
                    results_dir = fullfile(parentDir, sprintf('voids_results_%s', iteration));
                    if ~exist(results_dir, 'dir')
                        mkdir(results_dir);
                    end
                    % exportgraphics(ax1, fullfile(results_dir, sprintf('anisotropy_histogram_%d.png', i)));

                    exportgraphics(ax2, fullfile(results_dir, sprintf('void_area_fraction_%d.png', i)));

                    exportgraphics(ax3, fullfile(results_dir, sprintf('chain_area_fraction_%d.png', i)));

                    exportgraphics(ax5, fullfile(results_dir, sprintf('scalar_anisotropy_%d.png', i)));
                end

                % Save the figure
                annotation(f, 'textbox', [0 0.95 1 0.05], 'String', ...
                    sprintf('fs = %d, N = %d, t = %s', fs, N, timestamp_text), ...
                    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
                filename = fullfile(save_dir, sprintf('anisotropy_angle_%04d.png', i));
                set(f, 'Units', 'pixels', 'Position', [100, 100, 1200, 1000]);
                set(f, 'PaperUnits', 'inches');
                set(f, 'PaperPosition', [0 0 12 10]); % [left bottom width height] in inches
                set(f, 'PaperSize', [12 10]);         % [width height] in inches
                set(f, 'PaperPositionMode', 'manual');
                drawnow; % Ensure layout is updated % Force MATLAB to update layout
                
                print(f, filename, '-dpng', '-r100');
                close(f);
                progress = i / n_frames;
                waitbar(progress, h1);
                % break;
            end
            close(h1);
        end
        function plot_anisotropy_and_angle_evolution_stacks(obj)
            h1 = waitbar(0, 'Processing stacks for anisotropy and angle evolution');
            for i = 80:length(obj.app.stack_paths)%31,32
                % if i == 72
                %     continue;
                % end
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                % [N, fs] = obj.app.utils.get_info(obj.app.path);
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
                    fprintf('Skipping as empty %s\n', obj.app.path);
                    continue;
                end
                obj.app.utils.load_stack_info();
                obj.plot_anisotropy_and_angle_evolution_stack();
                progress = i / length(obj.app.stack_paths);
                waitbar(progress, h1);
            end
            close(h1);
        end
        function overlay_first_last_particle_locations(obj)

            % Ensure save directory exists
            save_dir = fullfile('F:', 'shake_table_data', 'Results', 'overlay_first_last');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
        
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                [N, fs] = obj.app.utils.get_info(obj.app.path);
                [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
        
                % Load stack_info
                obj.app.utils.load_stack_info();
                stack_info = obj.app.stack_info;
        
                % Find first and last image with non-empty particle locations
                first_idx = [];
                last_idx = [];
                for idx = stack_info.start_index:stack_info.end_index
                    ploc = obj.app.stack_info.particle_locations{idx};

                    if ~isempty(ploc)
                        if isempty(first_idx)
                            first_idx = idx;
                        end
                        last_idx = idx;
                    end
                end
                if isempty(first_idx) || isempty(last_idx)
                    fprintf('No valid particle locations in %s\n', obj.app.path);
                    continue;
                end
                fprintf('First index: %d, Last index: %d\n', first_idx, last_idx);
                % Load images
                img_files = stack_info.img_data.img_files;
                img1 = imread(fullfile(img_files(first_idx).folder, img_files(first_idx).name));
                img2 = imread(fullfile(img_files(last_idx).folder, img_files(last_idx).name));
                fprintf('Found images for first index %d and last index %d\n', first_idx, last_idx);
                % Get particle locations
                ploc1 = obj.app.particle_locator.get_particle_locations(first_idx);
                ploc2 = obj.app.particle_locator.get_particle_locations(last_idx);
                fprintf('Found particle locations for first index %d and last index %d\n', first_idx, last_idx);
                % Create figure
                f = figure('Visible', 'off', 'Position', [100, 100, 1000, 500]);
                subplot(1,2,1);
                imshow(img1, []);
                hold on;
                plot(ploc1.x, ploc1.y, 'ro', 'MarkerSize', 3, 'LineWidth', 0.5);
                title(sprintf('First (idx=%d)', first_idx));
                hold off;
                
                subplot(1,2,2);
                imshow(img2, []);
                hold on;
                plot(ploc2.x, ploc2.y, 'go', 'MarkerSize', 3, 'LineWidth', 0.5);
                title(sprintf('Last (idx=%d)', last_idx));
                hold off;
                fprintf('Overlaying particle locations for first and last images in %s\n', obj.app.path);
                % Save overlay image
                save_name = sprintf('overlay_first_last_N%d_f%d_iter%s.png', N, fs, iteration);
                exportgraphics(f, fullfile(save_dir, save_name));
                close(f);
            end
        end
        function convert_image_dirs_to_videos(obj)
            h1 = waitbar(0, 'Processing stacks for video conversion');
            % Converts directories of saved images for each stack into videos
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                [~, ~] = obj.app.utils.get_info(obj.app.path);
                [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
        
                % List of subdirectories to convert (add/remove as needed)
                subdirs = { ...
                    % sprintf('voids_images_%s', iteration), ...
                    sprintf('anisotropy_hist_%s', iteration), ...
                    sprintf('void_orientation_%s', iteration), ...
                    sprintf('anisotropy_angle_frames_%s', iteration) ...
                };
        
                for s = 1:length(subdirs)
                    img_dir = fullfile(parentDir, subdirs{s});
                    if ~exist(img_dir, 'dir')
                        continue;
                    end
                    % Find all PNG images in the directory
                    img_files = dir(fullfile(img_dir, '*.png'));
                    if isempty(img_files)
                        continue;
                    end
                    % Sort files by name (assumes zero-padded numbering)
                    [~, idx] = sort({img_files.name});
                    img_files = img_files(idx);
        
                    % Read the first image to get frame size
                    first_img = imread(fullfile(img_dir, img_files(1).name));
                    [height, width, ~] = size(first_img);
        
                    % Create VideoWriter object
                    video_filename = fullfile(img_dir, [subdirs{s} '.mp4']);
                    v = VideoWriter(video_filename, 'MPEG-4');
                    v.FrameRate = 10; % Adjust as needed
                    open(v);
        
                    % Write each image as a frame
                    for k = 1:length(img_files)
                        img = imread(fullfile(img_dir, img_files(k).name));
                        % Ensure all frames are the same size
                        if size(img,1) ~= height || size(img,2) ~= width
                            img = imresize(img, [height, width]);
                        end
                        writeVideo(v, img);
                    end
                    close(v);
                    fprintf('Saved video: %s\n', video_filename);
                end
                progress = i / length(obj.app.stack_paths);
                waitbar(progress, h1);
            end
            close(h1);
        end
        % plot chain and void area fraction vs absolute time, 
        % superpose all chains together for each frequency.
        function plot_chain_and_void_area_fraction_vs_time_combined(obj)
            % TODO: avg trials together for each frequency
            % Load the pooled voids data from results folder
            if isempty(obj.voids_data)
                obj.voids_data = load('F:\shake_table_data\Results\voids_data.mat');
            end
            % Save dir for the figures
            save_dir = fullfile('F:', 'shake_table_data', 'Results', 'area_fracs');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            % pool chain area fractions in to void_data
            obj.pool_chain_area_fraction_into_voidData();
            all_data = obj.voids_data.all_data;
            for ns = fieldnames(all_data)'
                n = ns{1};
                if n == "N4"
                    continue; % skip N4 for now
                end
                n_num = str2double(n(2:end)); % Extract stack number
                for fs = fieldnames(all_data.(n))'
                    f = fs{1};
                    f_num = str2double(f(2:end)); % Extract frequency number
                    for iters = fieldnames(all_data.(n).(f))'
                        iter = iters{1};
                        % if iter contains _cont skip it
                        if contains(iter, '_cont')
                            continue;
                        end
                        iter_num = str2double(iter(5:end)); % Extract iteration number
                        data = all_data.(n).(f).(iter);
                        fprintf('Processing frequency %s, stack %s, iteration %s\n', f, n, iter);
                        % assignin('base', 'temp', data); % For debugging
                        empty = load(sprintf('F://shake_table_data//N%d//%dhz_hopperflow//60deg//10cm//stack_info_%d.mat',n_num,f_num,iter_num), '-mat', 'empty');
                        if ~empty.empty                  
                            fprintf('Skipping %s\n', obj.app.path);           
                            continue;
                        end

                        figure(f_num);ax = gca;hold(ax, 'on');
                        % Plot void area fraction
                        plot(ax, data.normalized_x, data.void_area_frac,'color',obj.app.utils.get_color(n_num), 'LineWidth', 2, 'DisplayName', sprintf('%s %s %s', f, n, iter));

                        figure(f_num+21);ax2 = gca;hold(ax2, 'on');                     
                        % % Plot chain area fraction
                        plot(ax2, data.normalized_x, data.chain_area_frac,'color',obj.app.utils.get_color(n_num), 'LineWidth', 2, 'DisplayName', sprintf('%s %s %s', f, n, iter));
                    end
                end                
            end

            Ns = [12, 24, 48];
            legend_labels = arrayfun(@(n) sprintf('N = %d', n), Ns, 'UniformOutput', false);
            for ns = fieldnames(all_data)'
                n = ns{1};
                if n == "N4"
                    continue; % skip N4 for now
                end
                for fs = fieldnames(all_data.(n))'
                    f = fs{1};
                    f_num = str2double(f(2:end)); % Extract frequency number
                    figure(f_num);ax = gca;hold(ax, 'on');

                    figure(f_num+21);ax2 = gca;hold(ax2, 'on');                     

                    title(ax, 'Void Area Fraction Evolution');
                    ax.XLabel.String = 'Percent Completion (%)';
                    ax.YLabel.String = 'Void Area Fraction';
                    ax.XLim = [0 100];
                    % ax.YLim = [0 0.4];
                    legend_handles = gobjects(1, numel(Ns));
                    for k = 1:numel(Ns)
                        legend_handles(k) = plot(ax, NaN, NaN, 'LineWidth', 2, ...
                            'Color', obj.app.utils.get_color(Ns(k)), 'DisplayName', legend_labels{k});
                    end
                    legend(ax, legend_handles, legend_labels, 'Location', 'best');

                    title(ax2, 'Chains Area Fraction Evolution');
                    ax2.XLabel.String = 'Percent Completion (%)';
                    ax2.YLabel.String = 'Chains Area Fraction';
                    ax2.XLim = [0 100];
                    % ax2.YLim = [0 1];
                    legend_handles2 = gobjects(1, numel(Ns));
                    for k = 1:numel(Ns)
                        legend_handles2(k) = plot(ax2, NaN, NaN, 'LineWidth', 2, ...
                            'Color', obj.app.utils.get_color(Ns(k)), 'DisplayName', legend_labels{k});
                    end
                    legend(ax2, legend_handles2, legend_labels, 'Location', 'best');

                    exportgraphics(ax, fullfile(save_dir, sprintf('void_area_fraction_f%d.png', f_num)));
                    exportgraphics(ax2, fullfile(save_dir, sprintf('chain_area_fraction_f%d.png', f_num)));
                end                
            end          
        end
        function plot_chain_and_void_area_fraction_vs_time_averaged(obj)
            if isempty(obj.voids_data)
                obj.voids_data = load('F:\shake_table_data\Results\voids_data.mat');
            end
            save_dir = fullfile('F:', 'shake_table_data', 'Results', 'area_fracs','averaged');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            obj.pool_chain_area_fraction_into_voidData();
            all_data = obj.voids_data.all_data;
            Ns = [12, 24, 48];
            freq_list = [];
            % Find all unique frequencies
            for ns = fieldnames(all_data)'
                n = ns{1};
                for fs = fieldnames(all_data.(n))'
                    f = fs{1};
                    f_num = str2double(f(2:end));
                    freq_list = [freq_list, f_num];
                end
            end
            unique_fs = unique(freq_list);
            disp(sprintf('Unique frequencies found: %s', num2str(unique_fs)));

            % Common x-axis for interpolation
            xq = 0:1:100;

            for n_idx = 1:numel(Ns)
                N_val = Ns(n_idx);
                n_field = sprintf('N%d', N_val);
                if ~isfield(all_data, n_field)
                    fprintf('No data for N = %d\n', N_val);
                    continue;
                end
                for f_idx = 1:numel(unique_fs)
                    f_val = unique_fs(f_idx);
                    f_field = sprintf('f%d', f_val);
                    if ~isfield(all_data.(n_field), f_field)
                        fprintf('No data for frequency %d in N = %d\n', f_val, N_val);
                        continue;
                    end
                    % Gather all iterations for this (N, f)
                    iters = fieldnames(all_data.(n_field).(f_field));
                    void_mat = [];
                    chain_mat = [];
                    for it = 1:numel(iters)
                        iter = iters{it};
                        data = all_data.(n_field).(f_field).(iter);
                        col_names = data.Properties.VariableNames;
                        % Interpolate to common xq
                        if ismember('normalized_x', col_names) && ismember('void_area_frac', col_names) && ismember('chain_area_frac', col_names)
                            y_void = interp1(data.normalized_x, data.void_area_frac, xq, 'linear', 'extrap');
                            y_chain = interp1(data.normalized_x, data.chain_area_frac, xq, 'linear', 'extrap');
                            void_mat = [void_mat; y_void];
                            chain_mat = [chain_mat; y_chain];
                        end
                    end
                    if isempty(void_mat)
                        fprintf('No valid data for N = %d, f = %d\n', N_val, f_val);
                        continue;
                    end
                    % Compute mean and std
                    void_mean = mean(void_mat, 1, 'omitnan');
                    void_std = std(void_mat, 0, 1, 'omitnan');
                    chain_mean = mean(chain_mat, 1, 'omitnan');
                    chain_std = std(chain_mat, 0, 1, 'omitnan');
                    
                    indxs = 1:4:100;
                    % Plot
                    figure(100+f_val);
                    ax = gca; hold(ax, 'on');
                    % fill([xq fliplr(xq)], [void_mean+void_std fliplr(void_mean-void_std)], ...
                    %     obj.app.utils.get_color(N_val), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                    plot(xq, void_mean, 'LineWidth', 2, 'Color', obj.app.utils.get_color(N_val), ...
                        'DisplayName', sprintf('N = %d', N_val));
                    % add errorbars with handlevisibility off linestyle 'none' and marker 'none'
                    errorbar(xq(indxs), void_mean(indxs), void_std(indxs), 'LineStyle', 'none', 'Marker', 'none', ...
                        'Color', obj.app.utils.get_color(N_val), 'HandleVisibility', 'off');
                    xlabel('Percent Completion (%)');
                    ylabel('Void Area Fraction');
                    title(sprintf('Void Area Fraction (Averaged) for f = %d', f_val));
                    legend('show');
                    hold(ax, 'off');
                    exportgraphics(ax, fullfile(save_dir, sprintf('void_area_fraction_avg_N%d_f%d.png', N_val, f_val)));

                    figure(200+f_val);
                    ax2 = gca; hold(ax2, 'on');
                    % fill([xq fliplr(xq)], [chain_mean+chain_std fliplr(chain_mean-chain_std)], ...
                    %     obj.app.utils.get_color(N_val), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                    plot(xq, chain_mean, 'LineWidth', 2, 'Color', obj.app.utils.get_color(N_val), ...
                        'DisplayName', sprintf('N = %d', N_val));
                    % add errorbars with handlevisibility off linestyle 'none' and marker 'none'
                    errorbar(xq(indxs), chain_mean(indxs), chain_std(indxs), 'LineStyle', 'none', 'Marker', 'none', ...
                        'Color', obj.app.utils.get_color(N_val), 'HandleVisibility', 'off');
                    xlabel('Percent Completion (%)');
                    ylabel('Chain Area Fraction');
                    title(sprintf('Chain Area Fraction (Averaged) for f = %d', f_val));
                    legend('show');
                    hold(ax2, 'off');
                    exportgraphics(ax2, fullfile(save_dir, sprintf('chain_area_fraction_avg_N%d_f%d.png', N_val, f_val)));
                end
            end
        end
        function plot_anisotropy_norm_vs_time_combined(obj)
            % Load the pooled voids data from results folder
            if isempty(obj.voids_data)
                obj.voids_data = load('F:\shake_table_data\Results\voids_data.mat');
            end
            % Save dir for the figures
            save_dir = fullfile('F:', 'shake_table_data', 'Results', 'anisotropy_norm');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            obj.pool_anisotropy_norm_into_voidData();
            all_data = obj.voids_data.all_data;
            for ns = fieldnames(all_data)'
                n = ns{1};
                if n == "N4"
                    continue; % skip N4 for now
                end
                n_num = str2double(n(2:end)); % Extract stack number
                for fs = fieldnames(all_data.(n))'
                    f = fs{1};
                    f_num = str2double(f(2:end)); % Extract frequency number
                    for iters = fieldnames(all_data.(n).(f))'
                        iter = iters{1};
                        % if iter contains _cont skip it
                        if contains(iter, '_cont')
                            continue;
                        end
                        iter_num = str2double(iter(5:end)); % Extract iteration number
                        data = all_data.(n).(f).(iter);
                        fprintf('Processing frequency %s, stack %s, iteration %s\n', f, n, iter);
                        % assignin('base', 'temp', data); % For debugging
                        empty = load(sprintf('F://shake_table_data//N%d//%dhz_hopperflow//60deg//10cm//stack_info_%d.mat',n_num,f_num,iter_num), '-mat', 'empty');
                        if ~empty.empty                  
                            fprintf('Skipping %s\n', obj.app.path);           
                            continue;
                        end

                        figure(f_num);ax = gca;hold(ax, 'on');
                        % Plot void area fraction
                        plot(ax, data.normalized_x, data.anisotropy_norm,'color',obj.app.utils.get_color(n_num), 'LineWidth', 2, 'DisplayName', sprintf('%s %s %s', f, n, iter));
                    end
                end                
            end

            Ns = [12, 24, 48];
            legend_labels = arrayfun(@(n) sprintf('N = %d', n), Ns, 'UniformOutput', false);
            for ns = fieldnames(all_data)'
                n = ns{1};
                if n == "N4"
                    continue; % skip N4 for now
                end
                for fs = fieldnames(all_data.(n))'
                    f = fs{1};
                    f_num = str2double(f(2:end)); % Extract frequency number
                    figure(f_num);ax = gca;hold(ax, 'on');

                    figure(f_num+21);ax2 = gca;hold(ax2, 'on');                     

                    title(ax, 'Scaler anisotropy Evolution');
                    ax.XLabel.String = 'Percent Completion (%)';
                    ax.YLabel.String = 'Anisotropy';
                    ax.XLim = [0 100];
                    % ax.YLim = [0 0.4];
                    legend_handles = gobjects(1, numel(Ns));
                    for k = 1:numel(Ns)
                        legend_handles(k) = plot(ax, NaN, NaN, 'LineWidth', 2, ...
                            'Color', obj.app.utils.get_color(Ns(k)), 'DisplayName', legend_labels{k});
                    end
                    legend(ax, legend_handles, legend_labels, 'Location', 'best');

                    exportgraphics(ax, fullfile(save_dir, sprintf('anisotropy_f%d.png', f_num)));
                end                
            end
        end
        function plot_anisotropy_norm_vs_time_averaged(obj)
            if isempty(obj.voids_data)
                obj.voids_data = load('F:\shake_table_data\Results\voids_data.mat');
            end
            save_dir = fullfile('F:', 'shake_table_data', 'Results', 'anisotropy_norm', 'averaged');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            all_data = obj.voids_data.all_data;
            Ns = [12, 24, 48];
            freq_list = [];
            % Find all unique frequencies
            for ns = fieldnames(all_data)'
                n = ns{1};
                for fs = fieldnames(all_data.(n))'
                    f = fs{1};
                    f_num = str2double(f(2:end));
                    freq_list = [freq_list, f_num];
                end
            end
            unique_fs = unique(freq_list);

            % Common x-axis for interpolation
            xq = 0:1:100;

            for n_idx = 1:numel(Ns)
                N_val = Ns(n_idx);
                n_field = sprintf('N%d', N_val);
                if ~isfield(all_data, n_field)
                    fprintf('No data for N = %d\n', N_val);
                    continue;
                end
                for f_idx = 1:numel(unique_fs)
                    f_val = unique_fs(f_idx);
                    f_field = sprintf('f%d', f_val);
                    if ~isfield(all_data.(n_field), f_field)
                        fprintf('No data for frequency %d in N = %d\n', f_val, N_val);
                        continue;
                    end
                    % Gather all iterations for this (N, f)
                    iters = fieldnames(all_data.(n_field).(f_field));
                    norm_mat = [];
                    for it = 1:numel(iters)
                        iter = iters{it};
                        data = all_data.(n_field).(f_field).(iter);
                        % Check for required fields
                        if isfield(data, 'normalized_x') && isfield(data, 'anisotropy_norm')
                            y_norm = interp1(data.normalized_x, data.anisotropy_norm, xq, 'linear', 'extrap');
                            norm_mat = [norm_mat; y_norm];
                        end
                    end
                    if isempty(norm_mat)
                        fprintf('No valid anisotropy norm data for N = %d, f = %d\n', N_val, f_val);
                        continue;
                    end
                    % Compute mean and std
                    norm_mean = mean(norm_mat, 1, 'omitnan');
                    norm_std = std(norm_mat, 0, 1, 'omitnan');

                    % Plot
                    figure(f_val);
                    ax = gca; hold(ax, 'on');
                    indxs = 1:4:100;
                    % fill([xq fliplr(xq)], [norm_mean+norm_std fliplr(norm_mean-norm_std)], ...
                    %     obj.app.utils.get_color(N_val), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                    plot(xq, norm_mean, 'LineWidth', 2, 'Color', obj.app.utils.get_color(N_val), ...
                        'DisplayName', sprintf('N = %d', N_val));
                    % add errorbars with handlevisibility off linestyle 'none' and marker 'none'
                    errorbar(xq(indxs), norm_mean(indxs), norm_std(indxs), 'LineStyle', 'none', 'Marker', 'none', ...
                        'Color', obj.app.utils.get_color(N_val), 'HandleVisibility', 'off');
                    xlabel('Percent Completion (%)');
                    ylabel('Anisotropy Norm');
                    title(sprintf('Anisotropy Norm (Averaged) for f = %d', f_val));
                    legend('show');
                    hold(ax, 'off');
                    exportgraphics(ax, fullfile(save_dir, sprintf('anisotropy_norm_avg_N%d_f%d.png', N_val, f_val)));
                end
            end
        end
        function pool_chain_area_fraction_into_voidData(obj)
            % loop over all the stacks
            h1 = waitbar(0, 'Processing stacks for chain area fraction');
            for i = 30:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                [N, fs] = obj.app.utils.get_info(obj.app.path);
                [iteration, ~] = obj.app.utils.getIteration(obj.app.path);

                tbl = obj.voids_data.all_data.(sprintf('N%d', N)).(sprintf('f%d', fs)).(sprintf('iter%s', iteration));
                if ismember('chain_area_frac', tbl.Properties.VariableNames)
                    continue; % chain_area_frac already exists
                else
                    obj.app.utils.load_stack_info();
                    voids = obj.app.stack_info.voids;                
                    total_area = sum(sum(~obj.app.stack_info.mask));
                    data = obj.voids_data.all_data.(sprintf('N%d', N)).(sprintf('f%d', fs)).(sprintf('iter%s', iteration));
                    normalized_x = ((data.image_indexes-min(data.image_indexes)) / max(data.image_indexes)) * 100;
                
                    % chain_area_frac is not empty
                    if ~isfield(obj.app.stack_info, 'chain_area_frac')
                        n_frames = length(obj.app.stack_info.voids);
                        fprintf('Processing %d frames for anisotropy_and_angle_evolution\n', n_frames);
                        chain_area_frac = zeros(1, numel(data.image_indexes));
                        h2 = waitbar(0, 'Processing stack for anisotropy and angle evolution');
                        for ii = 1:n_frames
                            void_data = voids{ii};
                            if isempty(void_data) || ~isfield(void_data, 'anisotropy') || ~isfield(void_data, 'theta') || ~isfield(void_data, 'radius')
                                continue;
                            end

                            particle_locations = obj.app.particle_locator.get_particle_locations(ii);
                            radius = 5;
                            num_particles = numel(particle_locations.x);
                            circle_area = pi * radius^2;
                            chain_area = num_particles * circle_area;
                            chain_area_frac(end_index) = chain_area / total_area; % Calculate area fraction
                            
                            progress = ii / n_frames;
                            waitbar(progress, h2);
                        end
                        % add the chain area fraction to the data structure
                        obj.app.stack_info.chain_area_frac = chain_area_frac;
                        % save the stack info
                        obj.app.utils.save_stack_callback();
                        close(h2);
                    else
                        chain_area_frac = obj.app.stack_info.chain_area_frac;
                    end 
                    % --- Void area fraction ---
                    % obj.voids_data.(sprintf('f%d', fs)).(sprintf('N%d', N)).(sprintf('iter%s', iteration)).void_area_frac = data.void_area_frac;
                    obj.voids_data.all_data.(sprintf('N%d', N)).(sprintf('f%d', fs)).(sprintf('iter%s', iteration)).normalized_x = normalized_x;

                    % --- Chains area fraction ---
                    obj.voids_data.all_data.(sprintf('N%d', N)).(sprintf('f%d', fs)).(sprintf('iter%s', iteration)).chain_area_frac = chain_area_frac';
                    % obj.voids_data.(sprintf('f%d', fs)).(sprintf('N%d', N)).(sprintf('iter%s', iteration)).image_indexes = data.image_indexes;
                end
                progress = i / numel(obj.app.stack_paths);
                waitbar(progress, h1);
            end
            close(h1);
            % save the plot data
            all_data = obj.voids_data.all_data;
            save('F:\shake_table_data\Results\voids_data.mat', 'all_data');
        end
        function pool_anisotropy_norm_into_voidData(obj)
            % loop over all the stacks
            h1 = waitbar(0, 'Processing stacks for anisotropy norm vs time');
            for j = 30:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', j);
                obj.app.path = obj.app.stack_paths{j};
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end

                [N, fs] = obj.app.utils.get_info(obj.app.path);
                [iteration, ~] = obj.app.utils.getIteration(obj.app.path);

                tbl = obj.voids_data.all_data.(sprintf('N%d', N)).(sprintf('f%d', fs)).(sprintf('iter%s', iteration));
                if ismember('anisotropy_norm', tbl.Properties.VariableNames)
                    continue; % chain_area_frac already exists
                else
                    obj.app.utils.load_stack_info();
                    n_frames = length(obj.app.stack_info.voids);
                    voids = obj.app.stack_info.voids;
                    fprintf('Processing %d frames for anisotropy_and_angle_evolution\n', n_frames);

                    anisotropy_norm = [];
                                
                    h1 = waitbar(0, 'Processing stack for anisotropy and angle evolution');
                    for i = 1:n_frames
                        void_data = voids{i};
                        if isempty(void_data) || ~isfield(void_data, 'anisotropy') || ~isfield(void_data, 'theta') || ~isfield(void_data, 'radius')
                            continue;
                        end
                        void_data = obj.clean_void_data(void_data, 800, 800);
                        % --- Scalar Anisotropy plot ---
                        if isfield(void_data, 'norm_dev')
                            anisotropy_norm = [anisotropy_norm; void_data.norm_dev];
                        else
                            text(0.5, 0.5, 'No Anisotropy Data', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                        end
                        progress = i / n_frames;
                        waitbar(progress, h1);                    
                    end
                    obj.voids_data.all_data.(sprintf('N%d', N)).(sprintf('f%d', fs)).(sprintf('iter%s', iteration)).anisotropy_norm = anisotropy_norm;
                end
            end
            close(h1);
            % save the plot data
            all_data = obj.voids_data.all_data;
            save('F:\shake_table_data\Results\voids_data.mat', 'all_data');
        end
    end
end