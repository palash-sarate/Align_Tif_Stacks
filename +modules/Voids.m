classdef Voids < handle
    properties
        app
    end
    methods
        function obj = Voids(app)
            obj.app = app;
            obj.addUi();
        end
        function addUi(obj)
            % Add UI elements for voids detection
            obj.app.ui.controls.voidsButton = uicontrol('Style', 'pushbutton', 'String', 'Detect Voids', ...
                'Position', [20, 20, 100, 30], 'Callback', @(src, event) obj.detect_voids_callback(src, event));
        end
        function voids = find_voids(obj, image)
            % Find voids in the image using morphological operations

            % plot the image on obj.app.ui.controls.ax2
            imshow(image, 'Parent', obj.app.ui.controls.ax2);
        end
    end
end