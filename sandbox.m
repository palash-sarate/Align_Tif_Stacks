% stack_info = load("F:\shake_table_data\N4\4hz_hopperflow\60deg\10cm\stack_info_1.mat");
% while isfield(stack_info, 'stack_info')
%     stack_info = stack_info.stack_info;
% end
% %%
% bd = stack_info.bd;
% r_cut = 2 * bd;
% m = 6;
% pos = [stack_info.particle_locations.x, stack_info.particle_locations.y];

% % localOrderParameter computes the local bond-orientational order parameter.
% %
% %   psi = localOrderParameter(pos, m) computes the order parameter for
% %   each particle given in pos (an N-by-2 array of [x,y] positions) using the
% %   symmetry index m. For instance, use m = 4 for square symmetry or m = 6 for
% %   triangular (hexagonal) symmetry.
% %
% %   psi = localOrderParameter(pos, m, r_cut) only considers neighbors that are
% %   within the distance r_cut.
% %
% %   The local order parameter for particle i is defined as:
% %       psi_m(i) = | (1/N_i) * sum_{j in neighbors(i)} exp( i*m*theta_ij ) |
% %   where theta_ij is the angle between the vector (pos(j,:) - pos(i,:)) and the x-axis.
% %
% %   Example:
% %       % pos: N-by-2 array of particle coordinates
% %       psi_top = localOrderParameter(pos(pos(:,2)>y_thresh,:), 4);
% %       psi_bottom = localOrderParameter(pos(pos(:,2)<=y_thresh,:), 6);
% %
% %   See also delaunayTriangulation, atan2.
% %

% useCutoff = true;

% % Create a Delaunay triangulation from the particle positions.
% dt = delaunayTriangulation(pos(:,1), pos(:,2));

% % For each particle, find the indices of attached triangles (neighbors)
% attachList = vertexAttachments(dt);

% N = size(pos,1);
% psi = zeros(N,1);

% for i = 1:N
%     % Extract the indices of triangles attached to particle i
%     triIndices = attachList{i};
%     % Get all vertices from these triangles
%     nb = unique(dt.ConnectivityList(triIndices,:));
%     % Remove the particle itself from its neighbor list
%     nb(nb == i) = [];
    
%     % If a cutoff distance is provided, filter the neighbors by distance.
%     if useCutoff && ~isempty(nb)
%         distances = sqrt(sum((pos(nb,:) - pos(i,:)).^2, 2));
%         nb = nb(distances <= r_cut);
%     end
    
%     % If no neighbors are found, set order parameter to NaN.
%     if isempty(nb)
%         psi(i) = NaN;
%         continue;
%     end
    
%     % Calculate the angle between particle i and each of its neighbors.
%     angles = atan2(pos(nb,2) - pos(i,2), pos(nb,1) - pos(i,1));
    
%     % Compute the local order parameter for particle i.
%     psi(i) = abs(sum(exp(1i*m*angles)) / numel(angles));
% end

% %% Plot histogram of psi
% histogram(psi, 'Normalization', 'probability');
% xlabel('\psi');
% ylabel('Probability');
% title('Histogram of \psi');

% %% get diameter and center of boundaries in ff_voids
% ff_voids = voids{start_index};
% boundaries = ff_voids.B;
% % calculate the diameter of each boundary
% % and the center of each boundary
% ff_voids.diameters = zeros(length(boundaries), 1);
% ff_voids.centers = zeros(length(boundaries), 2);

% for i = 1:length(boundaries)
%     % get the coordinates of the boundary
%     boundary = boundaries{i};
%     % calculate the diameter of the boundary
%     % ff_voids.diameters(i) = max(pdist(boundary));
%     % Fit a minimum enclosing circle to the boundary
%     [center, radius] = minboundcircle(boundary(:,1), boundary(:,2));
%     ff_voids.diameters(i) = 2 * radius; % Diameter is twice the radius
%     % calculate the center of the boundary
%     ff_voids.centers(i, :) = mean(boundary, 1);
% end
% ff_voids.diameters = ff_voids.diameters/ 7;

%% code to analyse voids in any BW image and overlay their eigen vectors on the image

image_path = "F:\shake_table_data\Results\voids_images\N48_voids.png";
image = imread(image_path);
if ndims(image) == 3
    image = rgb2gray(image);
end
% check if the image is binary 
if ~islogical(image)
    disp('Image is not binary, converting to binary...');
    image = imbinarize(image);
end

overlay_eigen_vectors(image);

function overlay_eigen_vectors(binaryImage)
    % get the size of the image
    % [height, width] = size(binaryImage);
    % Get the boundaries of the voids
    [B,~,~,~] = bwboundaries(binaryImage, 'noholes');
    % only keep boundaries that don't touch the edge of the image
    % B = remove_holes_on_edge(B, height, width);
    % remove the holes that are too small
    B = B(cellfun(@(x) length(x) > 30, B)); % change 30 to a suitable threshold for your application %
    N = length(B);
    % Show the binary image
    imshow(binaryImage);
    overlay_boundaries(B);
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
    ax = gca;
    exportgraphics(ax, 'voids_eigenvectors.eps', 'ContentType', 'vector');
end
function overlay_boundaries(B)
    colors=['b' 'g' 'r' 'c' 'm' 'y'];
    hold('on');
    for k=1:length(B)
        boundary = B{k};
        cidx = mod(k,length(colors))+1;
        plot( boundary(:,2), boundary(:,1),...
            colors(cidx),'LineWidth', 2, 'HandleVisibility', 'off');
    
        %randomize text position for better visibility
    %   rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
    %   col = boundary(rndRow,2); row = boundary(rndRow,1);
    %   h = text(col+1, row-1, num2str(L(row,col)));
    %   set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
    end
    hold('off');
end