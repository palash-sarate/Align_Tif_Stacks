stack_info = load("F:\shake_table_data\N4\4hz_hopperflow\60deg\10cm\stack_info_1.mat");
while isfield(stack_info, 'stack_info')
    stack_info = stack_info.stack_info;
end
%%
bd = stack_info.bd;
r_cut = 2 * bd;
m = 6;
pos = [stack_info.particle_locations.x, stack_info.particle_locations.y];

% localOrderParameter computes the local bond-orientational order parameter.
%
%   psi = localOrderParameter(pos, m) computes the order parameter for
%   each particle given in pos (an N-by-2 array of [x,y] positions) using the
%   symmetry index m. For instance, use m = 4 for square symmetry or m = 6 for
%   triangular (hexagonal) symmetry.
%
%   psi = localOrderParameter(pos, m, r_cut) only considers neighbors that are
%   within the distance r_cut.
%
%   The local order parameter for particle i is defined as:
%       psi_m(i) = | (1/N_i) * sum_{j in neighbors(i)} exp( i*m*theta_ij ) |
%   where theta_ij is the angle between the vector (pos(j,:) - pos(i,:)) and the x-axis.
%
%   Example:
%       % pos: N-by-2 array of particle coordinates
%       psi_top = localOrderParameter(pos(pos(:,2)>y_thresh,:), 4);
%       psi_bottom = localOrderParameter(pos(pos(:,2)<=y_thresh,:), 6);
%
%   See also delaunayTriangulation, atan2.
%

useCutoff = true;

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

%% Plot histogram of psi
histogram(psi, 'Normalization', 'probability');
xlabel('\psi');
ylabel('Probability');
title('Histogram of \psi');

%% get diameter and center of boundaries in ff_voids
ff_voids = voids{start_index};
boundaries = ff_voids.B;
% calculate the diameter of each boundary
% and the center of each boundary
ff_voids.diameters = zeros(length(boundaries), 1);
ff_voids.centers = zeros(length(boundaries), 2);

for i = 1:length(boundaries)
    % get the coordinates of the boundary
    boundary = boundaries{i};
    % calculate the diameter of the boundary
    % ff_voids.diameters(i) = max(pdist(boundary));
    % Fit a minimum enclosing circle to the boundary
    [center, radius] = minboundcircle(boundary(:,1), boundary(:,2));
    ff_voids.diameters(i) = 2 * radius; % Diameter is twice the radius
    % calculate the center of the boundary
    ff_voids.centers(i, :) = mean(boundary, 1);
end
ff_voids.diameters = ff_voids.diameters/ 7;
