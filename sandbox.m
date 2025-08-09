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
        x_mean = mean(x);
        y_mean = mean(y);
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
        anisotropy(k) = 1 - (lambda2 / lambda1); %shape anisotropy
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
        fabric = fabric + (v * v');  % outer product
    end
    fabric = fabric / size(majorEigenVectors, 1);

    % --- Deviatoric fabric ---
    trace_fabric = trace(fabric);
    hydro = trace_fabric / 2;
    dev_fabric = fabric - hydro * eye(2);
    f_dev = dev_fabric * (4);  % scaling factor (similar to 3D: 15/2)

    % --- Anisotropy norm (Frobenius norm of deviatoric tensor) ---
    norm_dev = sqrt((2) * sum(sum(f_dev .* f_dev)));
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

figure
void_data = obj.app.stack_info.voids{image_index}; % assuming image_index is defined
 % --- Anisotropy plot ---
subplot(1,3,1);
if isfield(void_data, 'anisotropy')
    histogram(void_data.anisotropy, 20, 'FaceColor', [0.2 0.6 0.8], 'Normalization', 'pdf');
    xlabel('Anisotropy');
    ylabel('PDF');
    title(sprintf('Anisotropy Distribution (Frame %d)', i));
    xlim([0 1]);
end

% --- Angle distribution (polar plot) ---
subplot(1,3,2);
if isfield(void_data, 'theta') && isfield(void_data, 'radius')
    theta = [void_data.theta, void_data.theta+pi];
    radius = [void_data.radius, void_data.radius];
    polarplot(theta, radius, 'r-', 'LineWidth', 2);
    title('Angle Distribution');
end

% plot the anisotropy norm
if isfield(void_data, 'norm_dev')
    anisotropy_norm = [anisotropy_norm; void_data.norm_dev];
    anisotropy_norm_frame = [anisotropy_norm_frame; image_index];
end
% --- Anisotropy norm plot ---
subplot(1,3,3);
if ~isempty(anisotropy_norm)
    plot(anisotropy_norm_frame, anisotropy_norm, 'b-', 'LineWidth', 2);
    xlabel('Frame Index');
    ylabel('Anisotropy Norm');
    title('Anisotropy Norm Evolution');
    xlim([1 n_frames]);
end
            % directory of the tif stacks
            % folder = 'F:\shake_table_data\';
            % populate the list of paths to the tiff stacks
            Ns = [4,12,24,48];
            fs = [4,6,8,10,12,14,16,18,20];
            deg = 60;
            wd = 10;
            stack_paths = [];
            % fps = 47;

            for n = 1:length(Ns)
                for freq = 1:length(fs)
                    for w = 1:length(wd)
                        % Construct the path with escaped backslashes
                        c_path = sprintf("F:\\shake_table_data\\N%d\\%dhz_hopperflow\\%ddeg\\%dcm\\",Ns(n),fs(freq),deg,wd);
                        % find folders in the path directory
                        subFolders = dir(c_path);
                        subFolders = subFolders([subFolders.isdir]);  % Keep only directories
                        % remove directories that have voids_images in them
                        subFolders = subFolders(~contains({subFolders.name}, 'voids_images'));  % Remove directories that contain 'voids_images'
                        subFolders = subFolders(~contains({subFolders.name}, 'anisotropy_angle_frames'));  % Remove directories that contain 'voids_images'
                        subFolders = subFolders(~contains({subFolders.name}, 'void_orientation'));  % Remove directories that contain 'voids_images'
                        subFolders = subFolders(~contains({subFolders.name}, 'voids_results'));  % Remove directories that contain 'voids_images'
                        subFolders = subFolders(~contains({subFolders.name}, 'anisotropy_hist'));  % Remove directories that contain 'voids_images'
                        subFolders = subFolders(~contains({subFolders.name}, 'area_fractions'));  % Remove directories that contain 'voids_images'
                        subFolders = subFolders(~contains({subFolders.name}, 'largest_void_area'));  % Remove directories that contain 'voids_images'
                        subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % Remove '.' and '..' directories

                        for k = 1:length(subFolders)
                            stack_paths = [stack_paths; fullfile(c_path, subFolders(k).name)];
                        end
                    end
                end
            end

            stack_paths = [stack_paths; "F:\shake_table_data\time_control\\1"];
            stack_paths = [stack_paths; "F:\shake_table_data\time_control\\2"];
            stack_paths = [stack_paths; "F:\shake_table_data\time_control\\3"];
        end