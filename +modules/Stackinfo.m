classdef Stackinfo < handle
    properties (Access = private)
        matfileObj % Handle to the matfile object
        cache
    end
    % Preloaded fields (small fields that are loaded into memory at initialization)
    properties
        start_index
        end_index
        parentDir
        iteration
        aligned
        shortened
        masked
        bd
        number_density
        empty
        nhood
        jammed
        mask
        mask_vertices
    end
    properties (Dependent)
        % Dependent properties (loaded on demand)
        displacements
        img_data
        particle_locations
        distances
        gr
        gr_bins
        timestamps
        psi_6
        bond_orders
        voids
    end
    methods
        % Constructor
        function obj = Stackinfo(matfilePath)
            if exist(matfilePath, 'file') ~= 2
                error('The specified .mat file does not exist: %s', matfilePath);
            end
            obj.matfileObj = matfile(matfilePath, 'Writable', true); % Load as matfile object
            obj.cache = struct();
            % print the stack info fields
            % fprintf('Stack info fields:\n');
            % fields = fieldnames(obj.matfileObj);
            % for i = 1:length(fields)
            %     fprintf('%s\n', fields{i});
            % end

            % Preload small fields into memory
            obj.start_index = obj.matfileObj.start_index;
            obj.end_index = obj.matfileObj.end_index;
            obj.parentDir = obj.matfileObj.parentDir;
            obj.iteration = obj.matfileObj.iteration;
            obj.aligned = obj.matfileObj.aligned;
            obj.shortened = obj.matfileObj.shortened;
            obj.masked = obj.matfileObj.masked;
            obj.bd = obj.matfileObj.bd;
            obj.number_density = obj.matfileObj.number_density;
            obj.empty = obj.matfileObj.empty;
            obj.nhood = obj.matfileObj.nhood;
            obj.jammed = obj.matfileObj.jammed;
            obj.mask = obj.matfileObj.mask;
            obj.mask_vertices = obj.matfileObj.mask_vertices;
            fprintf('loaded small fields into memory\n');
        end
        % Getter for dependent properties
        function value = get.displacements(obj)
            fprintf('Loading displacements from matfile\n');
            if isfield(obj.cache, 'displacements')
                fprintf('Loading displacements from cache\n');
                value = obj.cache.displacements;
            else
                value = obj.matfileObj.displacements;
                obj.cache.displacements = value;
            end
        end
        % Setter for displacements
        function set.displacements(obj, value)
            obj.matfileObj.displacements = value; % Update matfile
            obj.cache.displacements = value; % Update cache
        end
        function value = get.img_data(obj)
            fprintf('Loading img_data from matfile\n');
            if isfield(obj.cache, 'img_data')
                fprintf('Loading img_data from cache\n');
                value = obj.cache.img_data;
            else
                value = obj.matfileObj.img_data;
                obj.cache.img_data = value;
            end
        end
        % Setter for img_data
        function set.img_data(obj, value)
            obj.matfileObj.img_data = value; % Update matfile
            obj.cache.img_data = value; % Update cache
        end
        function value = get.particle_locations(obj)
            if isfield(obj.cache, 'particle_locations')
                value = obj.cache.particle_locations;
            else
                value = obj.matfileObj.particle_locations;
                obj.cache.particle_locations = value;
            end
        end
                % Setter for particle_locations
                function set.particle_locations(obj, value)
                    obj.matfileObj.particle_locations = value; % Update matfile
                    obj.cache.particle_locations = value; % Update cache
                end
        function value = get.distances(obj)
            if isfield(obj.cache, 'distances')
                value = obj.cache.distances;
            else
                value = obj.matfileObj.distances;
                obj.cache.distances = value;
            end
        end
        function set.distances(obj, value)
            obj.matfileObj.distances = value; % Update matfile
            obj.cache.distances = value; % Update cache
        end
        function value = get.gr(obj)
            if isfield(obj.cache, 'gr')
                value = obj.cache.gr;
            else
                value = obj.matfileObj.gr;
                obj.cache.gr = value;
            end
        end
        function set.gr(obj, value)
            obj.matfileObj.gr = value; % Update matfile
            obj.cache.gr = value; % Update cache
        end
        function value = get.gr_bins(obj)
            if isfield(obj.cache, 'gr_bins')
                value = obj.cache.gr_bins;
            else
                value = obj.matfileObj.gr_bins;
                obj.cache.gr_bins = value;
            end
        end
        function set.gr_bins(obj, value)
            obj.matfileObj.gr_bins = value; % Update matfile
            obj.cache.gr_bins = value; % Update cache
        end
        function value = get.timestamps(obj)
            if isfield(obj.cache, 'timestamps')
                value = obj.cache.timestamps;
            else
                value = obj.matfileObj.timestamps;
                obj.cache.timestamps = value;
            end
        end
        function set.timestamps(obj, value)
            obj.matfileObj.timestamps = value; % Update matfile
            obj.cache.timestamps = value; % Update cache
        end
        function value = get.psi_6(obj)
            if isfield(obj.cache, 'psi_6')
                value = obj.cache.psi_6;
            else
                value = obj.matfileObj.psi_6;
                obj.cache.psi_6 = value;
            end
        end
        function set.psi_6(obj, value)
            obj.matfileObj.psi_6 = value; % Update matfile
            obj.cache.psi_6 = value; % Update cache
        end
        function value = get.bond_orders(obj)
            if isfield(obj.cache, 'bond_orders')
                value = obj.cache.bond_orders;
            else
                value = obj.matfileObj.bond_orders;
                obj.cache.bond_orders = value;
            end
        end
        function set.bond_orders(obj, value)
            obj.matfileObj.bond_orders = value; % Update matfile
            obj.cache.bond_orders = value; % Update cache
        end
        function value = get.voids(obj)
            if isfield(obj.cache, 'voids')
                value = obj.cache.voids;
            else
                value = obj.matfileObj.voids;
                obj.cache.voids = value;
            end
        end
        function set.voids(obj, value)
            obj.matfileObj.voids = value; % Update matfile
            obj.cache.voids = value; % Update cache
        end

    end
end