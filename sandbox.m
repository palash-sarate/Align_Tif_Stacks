
center = [1,1];
radius = 10;
thickness = 2;
masked_image = stack_info.mask;
mask_value = 1;

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
theoretical_ring_area = pi * (outer_radius^2 - inner_radius^2);

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
    percentage = (unmasked_count / total_ring_pixels) * 100;
else
    percentage = 0;
end

 

