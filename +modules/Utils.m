classdef Utils < handle
    methods
        % function all_distances = calculate_all_distances(locations)
        %     x_gpu = gpuArray(locations.x);
        %     y_gpu = gpuArray(locations.y);
        %     n_particles = size(locations, 1);
        %     % Pre-allocate distances array in GPU memory
        %     all_distances = zeros(n_particles*(n_particles-1)/2, 1, 'gpuArray');
        %     % calculate all distances if not already calculated
        %     disp('calculating distances');
        %     % Create indices for vectorized distance calculation
        %     idx = 1;
        %     for i = 1:n_particles-1
        %         % Calculate distances between particle i and all particles j>i
        %         dx = x_gpu(i) - x_gpu(i+1:end);
        %         dy = y_gpu(i) - y_gpu(i+1:end);
                
        %         % Calculate Euclidean distances
        %         d = sqrt(dx.^2 + dy.^2);
                
        %         % Store in the pre-allocated array
        %         n_dists = length(d);
        %         all_distances(idx:idx+n_dists-1) = d;
        %         idx = idx + n_dists;
        %     end
        % end
        % function uniform_locations = get_uniform_distribution(~, particle_locations, n)
        %     % # find box corners
        %     xmin = int32(min(particle_locations.x));
        %     xmax = int32(max(particle_locations.x));
        %     ymin = int32(min(particle_locations.y));
        %     ymax = int32(max(particle_locations.y));
        %     x = randi([xmin, xmax], n, 1);
        %     y = randi([ymin, ymax], n, 1);
        %     uniform_locations = table(x, y);
        % end
        % function mark_dead_zone_callback(~,~)
        %     % get the start and end frame of the dead timeline and save it in the stack_info
        %     display_warning("Select the start frame of the dead timeline");
        %     wait_for_keypress("n");
        %     start_frame = round(get(slider, 'Value'));
        %     display_warning("Select the end frame of the dead timeline");
        %     wait_for_keypress("n");
        %     end_frame = round(get(slider, 'Value'));   
        %     stack_info.dead_zone = [start_frame, end_frame];
        %     save_stack_callback();
        % end
        % function lastNonConstIndex = findLastNonConstIndex(arr)
        %     displacements = stack_info.displacements;
        %     displacements(:,1) = displacements(:,1) - displacements(end,1);
        %     displacements(:,2) = displacements(:,2) - displacements(end,2);
        %     % find the last non zero index
        %     lastNonConstIndex = find(displacements(:,2), 1, 'last');
        % end
    end
end