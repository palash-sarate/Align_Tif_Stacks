classdef Voids < handle
    properties
        app
    end
    methods
        function obj = Voids(app)
            obj.app = app;
        end

        function detect_voids_callback(obj)
            % Find voids in the image using morphological operations
            particle_locations = obj.app.particle_locator.get_particle_locations(obj.app.current_image_idx);
            hold(obj.app.ui.controls.ax2, 'on');
            plot(obj.app.ui.controls.ax2, particle_locations.x, particle_locations.y, 'b*');
            hold(obj.app.ui.controls.ax2, 'off');
            
        end
    end
end