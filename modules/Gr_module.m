classdef Gr_module < handle
    properties
        stack_info
    end
    methods
        % constructor
        function obj = Gr_module(stack_info)
            obj.stack_info = stack_info;
        end
        function all_distances = calculate_all_distances(~, locations)
            x_gpu = gpuArray(locations.x);
            y_gpu = gpuArray(locations.y);
            n_particles = size(locations, 1);
            % Pre-allocate distances array in GPU memory
            all_distances = zeros(n_particles*(n_particles-1)/2, 1, 'gpuArray');
            % calculate all distances if not already calculated
            disp('calculating distances');
            % Create indices for vectorized distance calculation
            idx = 1;
            for i = 1:n_particles-1
                % Calculate distances between particle i and all particles j>i
                dx = x_gpu(i) - x_gpu(i+1:end);
                dy = y_gpu(i) - y_gpu(i+1:end);
                
                % Calculate Euclidean distances
                d = sqrt(dx.^2 + dy.^2);
                
                % Store in the pre-allocated array
                n_dists = length(d);
                all_distances(idx:idx+n_dists-1) = d;
                idx = idx + n_dists;
            end
        end
        function gr = calculate_gr(obj, r_max, dr)
            % calculate the radial distribution function
            % load the particle locations
            particle_locations = obj.stack_info.particle_locations;
            % calculate all distances if not already calculated
            if ~isfield(obj.stack_info, 'distances')
                disp('calculating distances');
                distances = calculate_all_distances(particle_locations);
                obj.stack_info.distances = distances;
            else 
                distances = obj.stack_info.distances;
            end
            % bins
            bin_centers = dr:1:r_max-dr;
            % calculate the histogram
            disp('calculating radial distribution function');
    
            % Option 1: Fast histogram approach (slightly less precise)
            % uniform_locations = get_uniform_distribution(particle_locations, 4000);
            % uniform_distances = calculate_all_distances(uniform_locations);
            % [counts, ~] = histcounts(distances, bins);
            % [uniform_counts, ~] = histcounts(uniform_distances, bins);
            % gr = counts ./ uniform_counts;
    
            % Option 2: For precise window counting (if dr != bin width/2)
            count = get_count(bin_centers, dr, distances);
            gr = normalize_count(count, particle_locations);
    
            obj.stack_info.gr = gr;
            obj.stack_info.gr_bins = bins;
            save_stack_callback();
        end
        function norm_count = normalize_count(~, counts, particle_locs)
            % 1.Divide your total count by N, the number of reference particles you considered
            %  -- probably the total number of particles in your data.
            norm_count = counts/numel(particle_locs.x);
            % 2. Divide by the area of the bin
            %  -- probably the area of the ring between r and r+dr.
            norm_count = norm_count ./ (2*pi*bin_centers*dr);
            % 3. Divide by the particle number density, i.e max possible count in the bin
            norm_count = norm_count ./ (2*pi*bin_centers*dr/pi*(bd/2)^2);
        end
        function counts = get_count(~, bins_centers, dr, distances)
            counts = zeros(1, numel(bins_centers));
            for i = 1:numel(bins_centers)
                bin_center = bin_centers(i);
                % Count elements in range [bin_center-dr, bin_center+dr]
                in_range = (distances >= bin_center-dr) & (distances < bin_center+dr);
                counts(i) = sum(in_range);
            end
        end
        function uniform_locations = get_uniform_distribution(~, particle_locations, n)
            % # find box corners
            xmin = int32(min(particle_locations.x));
            xmax = int32(max(particle_locations.x));
            ymin = int32(min(particle_locations.y));
            ymax = int32(max(particle_locations.y));
            x = randi([xmin, xmax], n, 1);
            y = randi([ymin, ymax], n, 1);
            uniform_locations = table(x, y);
        end
        function get_particle_locations(obj, image_path, save_path)
            % get the particle locations from the image
            py.track.find_particle_locations(image_path=image_path, diam=int32(5), max_iterations=int32(10), minmass=int32(1), separation=int32(5), save_path=save_path);
            % load the saved csv from save_path
            particle_locations = readtable(save_path);
            obj.stack_info.particle_locations = particle_locations;
            save_stack_callback();
            % assignin('base', 'particle_locations', particle_locations);
        end
        function get_Gr(~,~)
            radius = 15 * 5;
            start = stack_info.start_index;
            image_path = fullfile(stack_info.img_data.img_files(start).folder, stack_info.img_data.img_files(start).name);
            [iter, parentDir] = getIteration(path);
            save_path = fullfile(parentDir, sprintf('particle_locations_%s.csv', iter));
    
            % if save_path doesn't exist, get the particle locations from the first image
            if ~isfile(save_path)
                disp('getting particle locations');
                get_particle_locations(image_path, save_path);
            end
    
            % check if gr exists in stack_info
            disp('calculating gr');
            % calculate gr for the first image
            bin_width = 10;
            calculate_gr(radius, bin_width);
            % plot the gr
            plot_gr(iter, parentDir);
            fprintf('Gr calculated for stack %s\n', path);
        end
    
        function plot_gr(iter, parentDir)
            % clear axis
            cla(ui.ax2);
            % plot the gr on ui.ax2
            plot(ui.ax2, stack_info.gr_bins(1:end-1), stack_info.gr);
            title('Radial distribution function');
            xlabel('r');
            ylabel('g(r)');
            axis(ui.ax2, 'tight');
            % save the plot
            temp = figure(visible='off');
            plot(stack_info.gr_bins(1:end-1), stack_info.gr);
            title('Radial distribution function');
            xlabel('r');
            ylabel('g(r)');
            saveas(temp, fullfile(parentDir, sprintf('gr_%s.png', iter)));
            close(temp);
        end
        function plot_all_gr(~,~)
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
            cla(ui.ax2);hold on;
            % iterate over all the stacks
            for i = 1:4%length(stack_paths)
                set(ui.stack_dropdown, 'Value', i);
                current_idx = get(ui.stack_dropdown, 'Value');
                path = stack_paths{current_idx};
                [iteration, parentDir] = getIteration(path);
                [N, fs] = get_info(path);
                Ns = [Ns, N];
                if contains(path, 'time_control') || contains(path, 'temp') || contains(path, 'cont')
                    fprintf('Skipping %s\n', path);           
                    continue;
                end
                % check if gr_all_stacks has gr, gr_bins data for N, fs, iteration
                if isfield(gr_all_stacks, sprintf('N%d', N)) && ...
                        isfield(gr_all_stacks.(sprintf('N%d', N)), sprintf('F%d', fs)) && ...
                        isfield(gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)), sprintf('Iter%s', iteration))
                    fprintf('gr data found in gr_all_stacks for %s\n', path);
                    gr = gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr');
                    gr_bins = gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr_bins');
                else
                    if exist(sprintf('%s//stack_info_%s.mat', parentDir, iteration), 'file')
                        load_stack_info();
                        if isfield(stack_info, 'gr')
                            fprintf('gr data found in stack_info for %s\n', path);
                            % add the data to gr_all_stacks
                            gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr') = stack_info.gr;
                            gr_all_stacks.(sprintf('N%d', N)).(sprintf('F%d', fs)).(sprintf('Iter%s', iteration)).('gr_bins') = stack_info.gr_bins;
                            gr = stack_info.gr;
                            gr_bins = stack_info.gr_bins;
                        else
                            fprintf('gr data not found in stack_info for %s calculating now\n', path);
                            get_Gr('mode', 'auto');
                        end
                    end
                end
                % plot the gr on ui.ax2
                plot(ui.ax2, gr_bins(1:end-1)/7.5, gr, "DisplayName", "None", "Color", get_color(N, [4,12,24,48]));                        
            end
            % add legend with different colors for each N
            unique_Ns = unique(Ns);
            legend_handles = zeros(1, numel(unique_Ns));
            legend_entries = cell(1, numel(unique_Ns));
            
            for j = 1:numel(unique_Ns)
                N_val = unique_Ns(j);
                % Create a "dummy" line just for the legend with the right color
                legend_handles(j) = plot(ui.ax2, NaN, NaN, '-', 'Color', get_color(N_val, unique_Ns));
                legend_entries{j} = sprintf('N = %d', N_val);
            end
            
            % Create legend using only our dummy lines
            legend(ui.ax2, legend_handles, legend_entries, 'Location', 'best');
            hold off;
            axis(ui.ax2, 'tight');
            title('Radial distribution function');
            xlabel('r/bd');
            ylabel('g(r)');
            % save the gr_all_stacks
            save('F:\shake_table_data\Results\gr_all_stacks.mat', 'gr_all_stacks');
        end
        function Gr_all_stacks_callback(~,~)
            % iterate over all the stacks
            for i = 1:length(stack_paths)
                set(ui.stack_dropdown, 'Value', i);
                load_images_callback();
                if isfield(stack_info, 'gr') || contains(path, 'time_control')
                    fprintf('Gr already exists for stack %s\n', path);
                    if ~forced
                        continue;
                    end
                end
                get_Gr('mode', 'auto');
            end
        end
    end
end