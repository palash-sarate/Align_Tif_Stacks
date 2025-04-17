classdef Steinhardt < handle
    properties
        app
    end
    methods
        function obj = Steinhardt(app)
            obj.app = app;
        end
        %%%%%%%%%%%%%%%%%%%%%% LOCAL ORDER PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%%
        function psi = localOrderParameter(~, pos, m, r_cut)
            % obj.localOrderParameter computes the local bond-orientational order parameter.
            %
            %   psi = obj.localOrderParameter(pos, m) computes the order parameter for
            %   each particle given in pos (an N-by-2 array of [x,y] positions) using the
            %   symmetry index m. For instance, use m = 4 for square symmetry or m = 6 for
            %   triangular (hexagonal) symmetry.
            %
            %   psi = obj.localOrderParameter(pos, m, r_cut) only considers neighbors that are
            %   within the distance r_cut.
            %
            %   The local order parameter for particle i is defined as:
            %       psi_m(i) = | (1/N_i) * sum_{j in neighbors(i)} exp( i*m*theta_ij ) |
            %   where theta_ij is the angle between the vector (pos(j,:) - pos(i,:)) and the x-axis.
            %
            %   Example:
            %       % pos: N-by-2 array of particle coordinates
            %       psi_top = obj.localOrderParameter(pos(pos(:,2)>y_thresh,:), 4);
            %       psi_bottom = obj.localOrderParameter(pos(pos(:,2)<=y_thresh,:), 6);
            %
            %   See also delaunayTriangulation, atan2.
            %
            
            if nargin < 2
                error('You must provide at least the positions and symmetry parameter m.');
            end
            
            % Optional neighbor cutoff: if not provided, use all neighbors from the triangulation.
            if nargin < 3
                useCutoff = false;
            else
                useCutoff = true;
            end
            
            % Create a Delaunay triangulation from the particle positions.
            dt = delaunayTriangulation(pos(:,1), pos(:,2));
            
            % For each particle, find the indices of attached triangles (neighbors)
            attachList = vertexAttachments(dt);
            
            N = size(pos,1);
            psi = zeros(N,1);
            
            for i = 1:N
                % Extract the indices of triangles attached to particle i
                triIndices = attachList{i};
                % Get all vertices from these triangles
                nb = unique(dt.ConnectivityList(triIndices,:));
                % Remove the particle itself from its neighbor list
                nb(nb == i) = [];
                
                % If a cutoff distance is provided, filter the neighbors by distance.
                if useCutoff && ~isempty(nb)
                    distances = sqrt(sum((pos(nb,:) - pos(i,:)).^2, 2));
                    nb = nb(distances <= r_cut);
                end
                
                % If no neighbors are found, set order parameter to NaN.
                if isempty(nb)
                    psi(i) = NaN;
                    continue;
                end
                
                % Calculate the angle between particle i and each of its neighbors.
                angles = atan2(pos(nb,2) - pos(i,2), pos(nb,1) - pos(i,1));
                
                % Compute the local order parameter for particle i.
                psi(i) = abs(sum(exp(1i*m*angles)) / numel(angles));
            end 
        end
        function psi = freudLocalOrderParameter(~, particles_path, bead_dia, l, nhood)
            psi = py.track.compute_ql(particles_path, bead_dia, int32(l), int32(nhood));
            % convert psi from python ndarray to matlab array
            psi = double(py.array.array('d', py.numpy.nditer(psi)));
        end
        % function to find the nhood where psi meand and std stabilise
        function find_nhood(obj)
            L = 6;
            nhood_values = [4, 6, 8, 10 , 12, 14, 16];
            colors = mat2cell(jet(length(nhood_values)), ones(1, length(nhood_values)), 3);
            cla(obj.app.ui.controls.ax2);
            for i = 1:length(nhood_values)
                nhood = nhood_values(i);
                % get the psi for the current nhood
                [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
                particles_path = fullfile(parentDir, sprintf('particle_locations_%s.csv', iter));
                if ~isfile(particles_path)
                    fprintf('Particle locations not found for stack %s\n', obj.app.path);
                    continue;
                end
                psi = obj.freudLocalOrderParameter(particles_path, obj.app.stack_info.bd, L, nhood);
                % plot the mean and std of psi
                mean_psi = mean(psi);
                std_psi = std(psi);
                % check if mean and std are stable
                if i > 1 && abs(mean_psi - mean_psi_prev) < 0.01 && abs(std_psi - std_psi_prev) < 0.01
                    fprintf('Stable psi found for nhood = %d\n', nhood);
                    break;
                end
                mean_psi_prev = mean_psi;
                std_psi_prev = std_psi;
                % plot the psi histogram on obj.app.ui.controls.ax2
                histogram(obj.app.ui.controls.ax2, psi, 'Normalization', 'pdf', 'FaceColor', colors{i}, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                hold(obj.app.ui.controls.ax2, 'on');
            end
            obj.app.stack_info.nhood = nhood;
            obj.app.utils.save_stack_callback();
            % save the obj.app.ui.controls.ax2
            exportgraphics(obj.app.ui.controls.ax2, fullfile(parentDir, sprintf('nhood_vs_psi_%s.png', iter)));
        end
        function find_all_nhood(obj)
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                obj.app.utils.load_stack_info();
                obj.find_nhood();
            end
        end
        function find_all_psi(obj)
            % iterate over all the stacks
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                obj.app.utils.load_stack_info();
                [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
                particles_path = fullfile(parentDir, sprintf('particle_locations_%s.csv', iter));
                if ~isfile(particles_path)
                    fprintf('Particle locations not found for stack %s\n', obj.app.path);
                    continue;
                end
                bead_dia = obj.app.stack_info.bd;
                nhood = obj.app.stack_info.nhood;
                l = 6;
                psi = obj.freudLocalOrderParameter(particles_path, bead_dia, l, nhood);
                obj.app.stack_info.(sprintf("psi_%d",l)) = psi;
                % save the psi to the stack_info
                % assignin('base', 'stack_info', stack_info);
                obj.app.utils.save_stack_callback();
                fprintf('Psi calculated for stack %s\n', obj.app.path);
            end
        end
        function get_LBOOP(obj)
            num_nebors = 6;
            l = 6;
            [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
            particles_path = fullfile(parentDir, sprintf('particle_locations_%s.csv', iter));
            if ~isfile(particles_path)
                fprintf('Particle locations not found for stack %s\n', obj.app.path);
                return;
            end
            bond_orders = py.track.compute_bond_orientation_order(particles_path, int32(l), int32(num_nebors));
            % convert bond_orders from python ndarray to matlab array
            bond_orders = double(py.array.array('d', py.numpy.nditer(bond_orders)));
            % save the bond_orders to the stack_info
            obj.app.stack_info.bond_orders = bond_orders;
            obj.app.utils.save_stack_callback();
        end
        function get_all_LBOOP(obj)
            % iterate over all the stacks
            for i = 1:length(obj.app.stack_paths)
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                obj.app.path = obj.app.stack_paths{i};
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp') || contains(obj.app.path, 'cont')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                fprintf('Calculating LBOOP for stack %s\n', obj.app.path);
                obj.app.utils.load_stack_info();
                obj.get_LBOOP();
            end
        end
        function plot_LBOOP(obj)
            % function to plot the current stack LBOOP on ax2
            % get the current stack info
            % [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
            data = obj.app.stack_info.bond_orders;
            % check if data is empty
            if isempty(data)
                fprintf('No data found for stack %s\n', obj.app.path);
                return;
            end
            % plot the data as a histogram on ax2
            cla(obj.app.ui.controls.ax2);
            histogram(obj.app.ui.controls.ax2, data, 'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            % set the x and y limits to be the same as the current axis
            xlim(obj.app.ui.controls.ax2, [0 1]);
            ylim(obj.app.ui.controls.ax2, [0 1]);
            % set the title and labels
            title(obj.app.ui.controls.ax2, 'LBOOP');
            xlabel(obj.app.ui.controls.ax2, 'LBOOP');
            ylabel(obj.app.ui.controls.ax2, 'Probability Density');
            % save the plot to the current stack folder
            [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
            exportgraphics(obj.app.ui.controls.ax2, fullfile(parentDir, sprintf('LBOOP_hist_%s.png', iter)));
        end
        function plot_LBOOP_per_chain(obj)
            cla(obj.app.ui.controls.ax2);
            hold(obj.app.ui.controls.ax2, 'on');
            % struct to hold the data for each chain
            % check if all_data struct exists in the base workspace
            if evalin('base', 'exist(''all_data'', ''var'')')
                all_data = evalin('base', 'all_data');
                all_data_exist = true;
            else
                % check if all_data struct at F:\\shake_table_data\\Results\\lboop_all_stacks.mat
                if isfile('F:\\shake_table_data\\Results\\lboop_all_stacks.mat')
                    all_data = load('F:\\shake_table_data\\Results\\lboop_all_stacks.mat');
                    while isfield(all_data, 'all_data')
                        all_data = all_data.all_data;
                    end
                    all_data_exist = true;
                else
                    all_data = struct();
                    all_data_exist = false;
                end
            end
            % get stack paths for the current N
            for i = 1:1:length(obj.app.stack_paths)
                obj.app.path = obj.app.stack_paths{i};
                % set the current stack in the dropdown
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                [N, ~] = obj.app.utils.get_info(obj.app.path); 
                % [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp') || contains(obj.app.path, 'cont')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                fprintf('Plotting LBOOP for stack %s\n', obj.app.path);
                if all_data_exist
                    data = all_data.(sprintf('i_%d', i));
                else
                    % load the stack info
                    obj.app.utils.load_stack_info();
                    data = obj.app.stack_info.bond_orders;
                    % store the data in the struct
                    all_data.(sprintf('i_%d',i)) = data;    
                end
                % check if data is empty
                if isempty(data)
                    fprintf('No data found for stack %s\n', obj.app.path);
                    return;
                end
                % Define common bin edges for all histograms
                common_edges = linspace(0, 1, 20); % Example: 20 bins between 0 and 1
                % plot the data as a histogram on ax2
                [counts, edges] = histcounts(data,common_edges, 'Normalization', 'pdf');
                counts = counts / sum(counts); % Normalize counts to get probability density
                edges = edges(2:end) - (edges(2)-edges(1))/2;
                plot(obj.app.ui.controls.ax2, edges, counts, '-o', 'Color',...
                 obj.app.utils.get_color(N), 'LineWidth', 2, 'DisplayName', 'None');
            end
            assignin('base', 'all_data', all_data);
            % save the all_data struct
            save('F:\\shake_table_data\\Results\\lboop_all_stacks.mat', 'all_data');
            % add legend with different colors for each N
            unique_Ns = [4,12,24,48];
            legend_handles = zeros(1, numel(unique_Ns));
            legend_entries = cell(1, numel(unique_Ns));
            
            for j = 1:numel(unique_Ns)
                N_val = unique_Ns(j);
                % Create a "dummy" line just for the legend with the right color
                legend_handles(j) = plot(obj.app.ui.controls.ax2, NaN, NaN, '-o', 'LineWidth', 2, 'Color', obj.app.utils.get_color(N_val));
                legend_entries{j} = sprintf('N = %d', N_val);
            end
            
            % Create legend using only our dummy lines
            legend(obj.app.ui.controls.ax2, legend_handles, legend_entries, 'Location', 'best');

            hold(obj.app.ui.controls.ax2, 'off');
            % set the x and y limits to be the same as the current axis
            xlim(obj.app.ui.controls.ax2, [0 1]);
            % ylim(obj.app.ui.controls.ax2, [0 1]);
            % set the title and labels
            title(obj.app.ui.controls.ax2, 'LBOOP');
            xlabel(obj.app.ui.controls.ax2, 'LBOOP');
            ylabel(obj.app.ui.controls.ax2, 'Probability Density');
            % save the plot to the current results folder F:\shake_table_data\Results\LBOOP
            exportgraphics(obj.app.ui.controls.ax2, fullfile('F:\shake_table_data\Results\', 'LBOOP_hist.png'));
        end
        function plot_avg_LBOOP_per_chain(obj)
            cla(obj.app.ui.controls.ax2);
            hold(obj.app.ui.controls.ax2, 'on');
            % struct to hold the data for each chain
            % check if all_data struct at F:\\shake_table_data\\Results\\lboop_all_stacks.mat
            if isfile('F:\\shake_table_data\\Results\\lboop_all_stacks.mat')
                all_data = load('F:\\shake_table_data\\Results\\lboop_all_stacks.mat');
                while isfield(all_data, 'all_data')
                    all_data = all_data.all_data;
                end
            else
                disp('No LBOOP data found for all stacks');
                return;
            end
            avg_data = struct();
            % get stack paths for the current N
            for i = 1:1:length(obj.app.stack_paths)
                obj.app.path = obj.app.stack_paths{i};
                % set the current stack in the dropdown
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                [N, ~] = obj.app.utils.get_info(obj.app.path); 
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp') || contains(obj.app.path, 'cont')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                fprintf('Processing LBOOP for stack %s\n', obj.app.path);

                data = all_data.(sprintf('i_%d', i));

                % check if data is empty
                if isempty(data)
                    fprintf('No data found for stack %s\n', obj.app.path);
                    continue;
                end
                % Define common bin edges for all histograms
                common_edges = linspace(0, 1, 20); % Example: 20 bins between 0 and 1
                [counts, edges] = histcounts(data, common_edges);
                counts = counts / sum(counts); % Normalize counts to get probability density
                edges = edges(2:end) - (edges(2)-edges(1))/2;

                % Initialize the structure for this N if it doesn't exist
                if ~isfield(avg_data, sprintf('N_%d', N))
                    avg_data.(sprintf('N_%d', N)).edges = edges;
                    avg_data.(sprintf('N_%d', N)).counts = [];
                end

                % Concatenate counts as columns for later mean and std calculation
                avg_data.(sprintf('N_%d', N)).counts = [avg_data.(sprintf('N_%d', N)).counts, counts'];
            end
        
            % Normalize the counts for each N
            for j = [4,12,24,48]
                % Calculate mean and standard deviation for the counts
                mean_counts = mean(avg_data.(sprintf('N_%d', j)).counts, 2);
                std_counts = std(avg_data.(sprintf('N_%d', j)).counts, 0, 2);% / sqrt(size(avg_data.(sprintf('N_%d', j)).counts, 2));
                edges = avg_data.(sprintf('N_%d', j)).edges;
                % Plot the mean counts with error bars
                errorbar(obj.app.ui.controls.ax2, edges, mean_counts, std_counts, '-o', ...
                    'Color', obj.app.utils.get_color(j), 'LineWidth', 2, 'DisplayName', sprintf('N = %d', j));
            end
            % save the all_data struct
            save('F:\\shake_table_data\\Results\\lboop_avg.mat', 'avg_data');
            % add legend with different colors for each N
            legend(obj.app.ui.controls.ax2, 'Location', 'best');
            hold(obj.app.ui.controls.ax2, 'off');
            % set the x and y limits to be the same as the current axis
            xlim(obj.app.ui.controls.ax2, [0 1]);
            % ylim(obj.app.ui.controls.ax2, [0 1]);
            % set the title and labels
            title(obj.app.ui.controls.ax2, 'LBOOP');
            xlabel(obj.app.ui.controls.ax2, 'LBOOP');
            ylabel(obj.app.ui.controls.ax2, 'Probability Density');
            % save the plot to the current results folder F:\shake_table_data\Results\LBOOP
            exportgraphics(obj.app.ui.controls.ax2, fullfile('F:\shake_table_data\Results\', 'LBOOP_hist_avg.png'));
        end
        function plot_avg_BOOP(obj)
            % load all data
            if isfile('F:\\shake_table_data\\Results\\lboop_all_stacks.mat')
                all_data = load('F:\\shake_table_data\\Results\\lboop_all_stacks.mat');
                while isfield(all_data, 'all_data')
                    all_data = all_data.all_data;
                end
            else
                disp('No LBOOP data found for all stacks');
                return;
            end
            avg_data = struct();
            % avg the lboop for each trial 
            for i = 1:1:length(obj.app.stack_paths)
                obj.app.path = obj.app.stack_paths{i};
                % set the current stack in the dropdown
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                [N, ~] = obj.app.utils.get_info(obj.app.path); 
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp') || contains(obj.app.path, 'cont')
                    fprintf('Skipping %s\n', obj.app.path);           
                    continue;
                end
                fprintf('Processing LBOOP for stack %s\n', obj.app.path);

                data = all_data.(sprintf('i_%d', i));

                % check if data is empty
                if isempty(data)
                    fprintf('No data found for stack %s\n', obj.app.path);
                    continue;
                end
            
                % Initialize the structure for this N if it doesn't exist
                if ~isfield(avg_data, sprintf('N_%d', N))
                    % avg_data.(sprintf('N_%d', N)).edges = edges;
                    avg_data.(sprintf('N_%d', N)).data = [];
                end

                % Concatenate counts as columns for later mean and std calculation
                avg_data.(sprintf('N_%d', N)).data = [avg_data.(sprintf('N_%d', N)).data, mean(data)];
            end
            assignin('base', 'avg_data', avg_data);
            % then avg the boop from each trial for each N
            cla(obj.app.ui.controls.ax2);
            hold(obj.app.ui.controls.ax2, 'on');
            for j = [4,12,24,48]
                % Calculate mean and standard deviation for the counts
                mean_boop = mean(avg_data.(sprintf('N_%d', j)).data, 2);
                std_boop = std(avg_data.(sprintf('N_%d', j)).data, 0, 2);
                % edges = avg_data.(sprintf('N_%d', j)).edges;
                % Plot the mean counts with error bars
                errorbar(obj.app.ui.controls.ax2, j, mean_boop^2, std_boop, '-o', ...
                    'Color', obj.app.utils.get_color(j), 'LineWidth', 2, 'DisplayName', sprintf('N = %d', j));
            end
            % add legend with different colors for each N
            legend(obj.app.ui.controls.ax2, 'Location', 'best');
            hold(obj.app.ui.controls.ax2, 'off');
            % set the x and y limits to be the same as the current axis
            % xlim(obj.app.ui.controls.ax2, [0 50]);
            % ylim(obj.app.ui.controls.ax2, [0.4 0.6]);
            % set the title and labels
            title(obj.app.ui.controls.ax2, 'BOOP');
            xlabel(obj.app.ui.controls.ax2, 'Chain length, N');
            ylabel(obj.app.ui.controls.ax2, 'BOOP^2');
            % save the plot to the current results folder F:\shake_table_data\Results\LBOOP
            exportgraphics(obj.app.ui.controls.ax2, fullfile('F:\shake_table_data\Results\', 'BOOP.png'));
        end
        function superpose_all_LBOOP(obj)
            % load all data
            if isfile('F:\\shake_table_data\\Results\\lboop_all_stacks.mat')
                all_data = load('F:\\shake_table_data\\Results\\lboop_all_stacks.mat');
                while isfield(all_data, 'all_data')
                    all_data = all_data.all_data;
                end
            else
                disp('No LBOOP data found for all stacks');
                return;
            end

            for i = 1:1:length(obj.app.stack_paths)
                obj.app.path = obj.app.stack_paths{i};
                % set the current stack in the dropdown
                set(obj.app.ui.controls.stackDropdown, 'Value', i);
                % [N, ~] = obj.app.utils.get_info(obj.app.path);
                if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp') || contains(obj.app.path, 'cont')
                    fprintf('Skipping %s\n', obj.app.path);
                    continue;
                end
                fprintf('Processing LBOOP for stack %s\n', obj.app.path);

                data = all_data.(sprintf('i_%d', i));
                
                % check if data is empty
                if isempty(data)
                    fprintf('No data found for stack %s\n', obj.app.path);
                    continue;
                end
            
                % load the stack info
                obj.app.utils.load_stack_info();
                particle_locations = obj.app.stack_info.particle_locations;
                x = particle_locations.x;
                y = particle_locations.y;
                y = 800-y;
                % plot a hexagon for each particle location where the color is the lboop value i.e data(i)
                temp_min = min(data);
                temp_max = max(data);
                normalized_temp = (data - temp_min) / (temp_max - temp_min);
                
                % Get colors from the jet colormap
                colormap_jet = jet(256); % 256 colors in the colormap
                color_indices = round(normalized_temp * 255) + 1; % Map normalized values to colormap indices
                colors = colormap_jet(color_indices, :); 
                cla(obj.app.ui.controls.ax2);
                hold(obj.app.ui.controls.ax2, 'on');
                % for j = 1:size(particle_locations, 1)
                    scatter(obj.app.ui.controls.ax2,...
                     x, y,...
                      25, colors, 'filled');
                % end
                hold(obj.app.ui.controls.ax2, 'off');
                % Add colorbar for reference
                colormap(jet);
                colorbar;
                clim([0, 1]);

                % save the plot to the stack folder
                [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
                exportgraphics(obj.app.ui.controls.ax2, fullfile(parentDir, sprintf('LBOOP_%s.png', iter)));
            end
        end
    end
end