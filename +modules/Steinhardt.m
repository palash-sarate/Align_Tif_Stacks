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
            for N  = [4]
                cla(obj.app.ui.controls.ax2);
                hold(obj.app.ui.controls.ax2, 'on');
                % get stack paths for the current N
                stack_paths = obj.app.stack_paths(contains(obj.app.stack_paths, sprintf('\N%d\', N)));
                for i = 1:3%length(stack_paths)
                    obj.app.path = stack_paths{i};
                    [~, f] = obj.app.utils.get_info(obj.app.path); 
                    % [iter, parentDir] = obj.app.utils.getIteration(obj.app.path);
                    if contains(obj.app.path, 'time_control') || contains(obj.app.path, 'temp') || contains(obj.app.path, 'cont')
                        fprintf('Skipping %s\n', obj.app.path);           
                        continue;
                    end
                    fprintf('Plotting LBOOP for stack %s\n', obj.app.path);
                    % load the stack info
                    obj.app.utils.load_stack_info();
                    data = obj.app.stack_info.bond_orders;
                    % check if data is empty
                    if isempty(data)
                        fprintf('No data found for stack %s\n', obj.app.path);
                        return;
                    end
                    % plot the data as a histogram on ax2
                    [counts, edges] = histcounts(data, 'Normalization','pdf');
                    edges = edges(2:end) - (edges(2)-edges(1))/2;
                    plot(edges, counts, 'Color', obj.app.utils.get_color(f), 'LineWidth', 2, 'DisplayName', sprintf('N = %d', N));
                end
                hold(obj.app.ui.controls.ax2, 'off');
                % set the x and y limits to be the same as the current axis
                xlim(obj.app.ui.controls.ax2, [0 1]);
                ylim(obj.app.ui.controls.ax2, [0 1]);
                % set the title and labels
                title(obj.app.ui.controls.ax2, 'LBOOP');
                xlabel(obj.app.ui.controls.ax2, 'LBOOP');
                ylabel(obj.app.ui.controls.ax2, 'Probability Density');
                % save the plot to the current results folder F:\shake_table_data\Results\LBOOP
                exportgraphics(obj.app.ui.controls.ax2, fullfile('F:\shake_table_data\Results\LBOOP', sprintf('LBOOP_hist_%d.png', N)));
            end
        end
    end
end