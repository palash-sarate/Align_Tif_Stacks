classdef RDF < handle
    properties
        app
    end
    methods
        function obj = RDF(app)
            obj.app = app;
        end
        %%%%%%%%%%%%%%%%%%%%%% GR %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Gr_all_stacks_callback(obj, ~,~)
            % iterate over all the stacks
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp') || contains(obj.app.path, 'cont')
                    continue;
                end
                obj.app.load_images_callback();
                if isfield(obj.app.stack_info, 'gr')
                    fprintf('Gr already exists for stack %s\n', obj.app.path);
                    if ~obj.app.forced
                        fprintf('Skipping stack %s\n', obj.app.path);
                        continue;
                    end
                end
                obj.get_Gr('mode', 'auto');
            end
        end
        function get_Gr(obj, ~,~)
            radius = 15 * 5;
            img_idx = obj.app.stack_info.start_index;
            [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
            particle_locations = obj.app.particle_locator.get_particle_locations(img_idx);

            % check if gr exists in stack_info
            if ~isfield(obj.app.stack_info, 'gr') || obj.app.forced
                disp('calculating gr');
                % calculate gr for the first image
                bin_width = 3;
                obj.calculate_gr(radius, bin_width, particle_locations);
                % assignin('base', 'gr', gr);
            end
            % plot the gr
            obj.plot_gr(iter, parentDir);
            fprintf('Gr calculated for stack %s\n', obj.app.path);
        end
        function gr = calculate_gr(obj, r_max, dr, particle_locations)
            % calculate the radial distribution function
            % bins
            bin_centers = dr:0.1:r_max-dr;
            % calculate the histogram
            disp('calculating radial distribution function');

            % Option 2: For precise window counting (if dr != bin width/2)
            counts = get_counts(bin_centers, dr, particle_locations);
            % fprintf('Counts before normalizing: %d\n', counts);
            gr = normalize_count(counts, particle_locations);

            obj.app.stack_info.gr = gr;
            obj.app.stack_info.gr_bins = bin_centers;
            % assignin('base', 'stack_info', stack_info);
            obj.app.utils.save_stack_callback();
            function counts = get_counts(bins_centers, dr, particle_locations)
                counts = zeros(1, numel(bins_centers));
                wait = waitbar(0, 'Calculating Count for each Particle');
                for i = 1:numel(particle_locations.x)       
                    % get the particle location
                    x = particle_locations.x(i);
                    y = particle_locations.y(i);
                    if obj.app.stack_info.mask(floor(y), floor(x)) == 1
                        continue;
                    end
                    % get the distances of the particle from all other particles
                    distances = sqrt((particle_locations.x - x).^2 + (particle_locations.y - y).^2);
                    for j = 1:numel(bins_centers)
                        bin_center = bins_centers(j);
                        % get the distances in the range [bin_center-dr, bin_center+dr]
                        in_range = (distances >= bin_center-dr) & (distances < bin_center+dr);
                        % normalize by the area of the bin with correction for border
                        in_range = in_range / get_sliceArea(x, y, bin_center, dr);
                        % add the count to the counts array
                        counts(j) = counts(j) + sum(in_range);                
                    end
                    waitbar(i/numel(particle_locations.x));
                end
                close(wait);
            end
            function area = get_sliceArea(x, y, bin_center, dr)
                % depending on the particle location, bin_center and dr, 
                % get the area of the bin considering the border and the mask
                % usually area = pi*(bin_center+dr)^2 - pi*(bin_center-dr)^2;
                % but if the bin is close to the border or the mask, the area will be less
                percent = calculate_unmasked_ring_percentage([y, x], bin_center, dr, obj.app.stack_info.mask);
                % if percent == 0
                %     fprintf('Percentage for x: %d, y: %d, bin: %d, percent: %d\n',x,y, bin_center, percent);
                % end
                area = (pi*(bin_center+dr)^2 - pi*(bin_center-dr)^2) * percent;
            end
            function percentage = calculate_unmasked_ring_percentage(center, radius, thickness, masked_image, mask_value)
                % Calculate the percentage of a complete circular ring that is both:
                % 1) Inside the image boundaries
                % 2) Not masked
                %
                % Parameters:
                % -----------
                % center : [row, col]
                %     Center coordinates of the circle in [row, column] format
                % radius : double
                %     Radius of the circle
                % thickness : double
                %     Thickness of the ring (dr), the ring extends from r-dr to r+dr
                % masked_image : 2D matrix
                %     Binary mask where mask_value indicates masked areas
                % mask_value : scalar, optional
                %     Value in masked_image that indicates masked areas (default: 0)
                %
                % Returns:
                % --------
                % percentage : double
                %     Percentage of the complete circular ring that is inside the image
                %     and not masked (0-100%)
                
                % Set default mask value if not provided
                if nargin < 5
                    mask_value = 1;
                end
                
                % Get image dimensions
                [height, width] = size(masked_image);
                
                % Extract center coordinates
                center_row = center(1);
                center_col = center(2);
                
                % Calculate inner and outer radius
                inner_radius = max(0, radius - thickness);
                outer_radius = radius + thickness;
                
                % Calculate the theoretical area of the complete circular ring
                % Area = π(R₂² - R₁²)
                % theoretical_ring_area = pi * (outer_radius^2 - inner_radius^2);
                
                % Create a large enough grid to cover the entire ring (ignoring image boundaries)
                row_min = floor(center_row - outer_radius);
                row_max = ceil(center_row + outer_radius);
                col_min = floor(center_col - outer_radius);
                col_max = ceil(center_col + outer_radius);
                
                % Create coordinate matrices for the full theoretical ring
                [cols, rows] = meshgrid(col_min:col_max, row_min:row_max);
                
                % Calculate distances from the center for all points
                dist_from_center = sqrt((rows - center_row).^2 + (cols - center_col).^2);
                
                % Create a mask for the full circular ring
                ring_mask = (dist_from_center >= inner_radius) & (dist_from_center <= outer_radius);
                
                % Create a mask for points that are inside the image boundaries
                valid_rows = (rows >= 1) & (rows <= height);
                valid_cols = (cols >= 1) & (cols <= width);
                inside_image_mask = valid_rows & valid_cols;
                
                % Count pixels in the theoretical complete ring
                total_ring_pixels = sum(ring_mask(:));
                
                % For pixels inside both the ring and image boundaries, check if they're masked
                valid_points = ring_mask & inside_image_mask;
                
                % Initialize a counter for unmasked pixels
                unmasked_count = 0;
                
                % Get linear indices of points inside both ring and image
                [valid_row_indices, valid_col_indices] = find(valid_points);
                
                % Check each valid point against the mask
                for i = 1:length(valid_row_indices)
                    img_row = valid_row_indices(i) + row_min - 1;  % Convert back to image coordinates
                    img_col = valid_col_indices(i) + col_min - 1;
                    
                    % Check if this point is not masked
                    if masked_image(img_row, img_col) ~= mask_value
                        unmasked_count = unmasked_count + 1;
                    end
                end
                
                % Calculate percentage: (unmasked pixels) / (total theoretical ring pixels) * 100
                if total_ring_pixels > 0
                    percentage = (unmasked_count / total_ring_pixels);
                else
                    percentage = 0;
                end
            end
            function norm_count = normalize_count(counts, particle_locs)
                % 1.Divide your total count by N, the number of reference particles you considered
                %  -- probably the total number of particles in your data.
                norm_count = counts/numel(particle_locs.x);
                % 2. Divide by the area of the bin
                %  -- probably the area of the ring between r and r+dr. 
                % NOTE: The area normalization is being done in the get_count function
                % norm_count = norm_count ./ get_slice_area(bin_centers, dr, particle_locs);
                % 3. Divide by the particle number density, i.e number/V
                if ~isfield(obj.app.stack_info, 'number_density')
                    calculate_number_density(particle_locs);
                end
                norm_count = norm_count ./ obj.app.stack_info.number_density;
            end
            function calculate_number_density(particle_locations)
                % calculate the number density of the particles
                
                % get the box dimensions
                xmin = min(particle_locations.x);
                xmax = max(particle_locations.x);
                ymin = min(particle_locations.y);
                ymax = max(particle_locations.y);
                % get the mask area
                mask_area = sum(obj.app.stack_info.mask(:));
                % calculate the area of the box
                area = (xmax - xmin) * (ymax - ymin) - mask_area;
                % calculate the number density
                number_density = numel(particle_locations.x) / area;
                obj.app.stack_info.number_density = number_density;
            end
        end
        function plot_local_number_density(obj)
            particle_locations = obj.app.particle_locator.get_particle_locations(obj.app.stack_info.start_index);
            % plot the local number density as a colormap only for 1st image
            if ~isfield(obj.app.stack_info, 'local_number_density')
                obj.calculate_local_number_density(particle_locations);
            end
            local_number_density = obj.app.stack_info.local_number_density;
            % obj.app.utils.save_stack_callback();
            % plot the local number density
            cla(obj.app.ui.controls.ax2);
            imagesc(obj.app.ui.controls.ax2, local_number_density);
            set(obj.app.ui.controls.ax2, 'YDir', 'normal'); % Ensure the y-axis is not flipped
            title('Local number density');
            colorbar;
            clim([-1, 1]); % Set colorbar limits to 0 and 1
        end
        function calculate_local_number_density(obj, particle_locations)
            bd = obj.app.stack_info.bd;
            % number density at each pixel is the...
            % number of particles in a box of size 20bd x 20bd
            % divided by the area of the box 20bd x 20bd - mask area
            obj.app.utils.setFrame(1);
            [image_height,image_width] = size(obj.app.stack_info.img_data.imgs{obj.app.stack_info.start_index});
            % fprintf('height: %d, width: %d\n', image_height, image_width);
            fprintf('Calculating local number density\n');
            local_number_density = zeros(image_height, image_width);
            wait = waitbar(0, 'Calculating local number density');
            for i = 1:image_height
                for j = 1:image_width
                    % get the box corners
                    xmin = max(1, j - floor(10*bd));
                    xmax = min(image_width, j + floor(10*bd));
                    ymin = max(1, i - floor(10*bd));
                    ymax = min(image_height, i + floor(10*bd));
                    % fprintf('xmin: %d, xmax: %d, ymin: %d, ymax: %d\n', xmin, xmax, ymin, ymax);
                    % get the particles in the box
                    particles_in_box = particle_locations.x >= xmin & particle_locations.x <= xmax & ...
                        particle_locations.y >= ymin & particle_locations.y <= ymax;
                    % calculate the number density
                    mask_area = sum(obj.app.stack_info.mask(xmin:xmax, ymin:ymax), 'all');
                    area = (xmax - xmin) * (ymax - ymin) - mask_area;
                    local_number_density(i,j) = sum(particles_in_box) / area;
                    waitbar(i/image_height);
                end
            end
            wait.close();
            obj.app.stack_info.local_number_density = local_number_density;
            assignin('base', 'stack_info', obj.app.stack_info);
            % obj.app.utils.save_stack_callback();
        end
        function plot_gr(obj, iter, parentDir)
            % clear axis
            cla(obj.app.ui.controls.ax2);
            % plot the gr on obj.app.ui.controls.ax2        
            plot(obj.app.ui.controls.ax2, obj.app.stack_info.gr_bins(1:end)/obj.app.stack_info.bd, obj.app.stack_info.gr);
            title('Radial distribution function');
            xlabel('r/bd');
            ylabel('g(r)');
            axis(obj.app.ui.controls.ax2, 'tight');
            % save the plot
            temp = figure(visible='off');
            plot(obj.app.stack_info.gr_bins(1:end)/obj.app.stack_info.bd, obj.app.stack_info.gr);
            title('Radial distribution function');
            xlabel('r/bd');
            ylabel('g(r)');
            saveas(temp, fullfile(parentDir, sprintf('gr_%s.png', iter)));
            close(temp);
        end
        function plot_all_gr(obj, ~,~)
            Ns = [];    
            % check if struct with gr of all stacks exists at F:\shake_table_data\Results
            if exist('F:\shake_table_data\Results\gr_all_stacks.mat', 'file')
                gr_all_stacks = load('F:\shake_table_data\Results\gr_all_stacks.mat');
                while isfield(gr_all_stacks, 'gr_all_stacks')
                    gr_all_stacks = gr_all_stacks.gr_all_stacks;
                end
            else
                gr_all_stacks = struct();
            end
            % clear axis
            cla(obj.app.ui.controls.ax2);hold on;
            % iterate over all the stacks
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                current_idx = get(obj.app.ui.controls.stackDropdown, 'Value');
                obj.app.path = obj.app.stack_paths{current_idx};
                [iteration, parentDir] = obj.app.utils.getIteration(obj.app.path);
                [N, fs] = obj.app.utils.get_info(obj.app.path);
                Ns = [Ns, N];
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp') || contains(obj.app.path, 'cont')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                if obj.app.trial.is_rec_start_late(N,fs,iteration)
                    fprintf('Ignoring %d,%d,%s\n as rec started late',N,fs,iteration);
                    continue;
                end
                % check if gr_all_stacks has gr, gr_bins data for N, fs, iteration
                if isfield(gr_all_stacks, sprintf('N%d', N)) && ...
                        isfield(gr_all_stacks.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                        isfield(gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration)) && ~obj.app.forced
                    fprintf('gr data found in gr_all_stacks for %s\n', obj.app.path);
                    gr = gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr');
                    gr_bins = gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr_bins');
                else
                    if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                        obj.app.utils.load_stack_info();
                        if isfield(obj.app.stack_info, 'gr')
                            fprintf('gr data found in stack_info for %s\n', obj.app.path);
                            % add the data to gr_all_stacks
                            gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr') = obj.app.stack_info.gr;
                            gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr_bins') = obj.app.stack_info.gr_bins;
                            gr = obj.app.stack_info.gr;
                            gr_bins = obj.app.stack_info.gr_bins;
                        else
                            fprintf('gr data not found in stack_info for %s calculating now\n', obj.app.path);
                            obj.get_Gr('mode', 'auto');
                        end
                    end
                end
                % plot the gr on obj.app.ui.controls.ax2
                plot(obj.app.ui.controls.ax2, gr_bins(1:end)/obj.app.stack_info.bd, gr, "DisplayName", "None", "Color", obj.app.utils.get_color(N));
            end
            % add legend with different colors for each N
            unique_Ns = unique(Ns);
            legend_handles = zeros(1, numel(unique_Ns));
            legend_entries = cell(1, numel(unique_Ns));
            
            for j = 1:numel(unique_Ns)
                N_val = unique_Ns(j);
                % Create a "dummy" line just for the legend with the right color
                legend_handles(j) = plot(obj.app.ui.controls.ax2, NaN, NaN, '-', 'Color', obj.app.utils.get_color(N_val));
                legend_entries{j} = sprintf('N = %d', N_val);
            end
            
            % Create legend using only our dummy lines
            legend(obj.app.ui.controls.ax2, legend_handles, legend_entries, 'Location', 'best');
            hold off;
            axis(obj.app.ui.controls.ax2, 'tight');
            title('Radial distribution function');
            xlabel('r/bd');
            ylabel('g(r)');
            % save the gr_all_stacks
            save('F:\shake_table_data\Results\gr_all_stacks.mat', 'gr_all_stacks');
            % save the plot
            exportgraphics(obj.app.ui.controls.ax2, 'F:\shake_table_data\Results\gr_all_stacks.png');
        end
    end
end
