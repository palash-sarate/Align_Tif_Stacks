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

            % Preload small fields into memory
            obj.start_index = obj.matfileObj.stack_info.start_index;
            obj.end_index = obj.matfileObj.stack_info.end_index;
            obj.parentDir = obj.matfileObj.stack_info.parentDir;
            obj.iteration = obj.matfileObj.stack_info.iteration;
            obj.aligned = obj.matfileObj.stack_info.aligned;
            obj.shortened = obj.matfileObj.stack_info.shortened;
            obj.masked = obj.matfileObj.stack_info.masked;
            obj.bd = obj.matfileObj.stack_info.bd;
            obj.number_density = obj.matfileObj.stack_info.number_density;
            obj.empty = obj.matfileObj.stack_info.empty;
            obj.nhood = obj.matfileObj.stack_info.nhood;
            obj.jammed = obj.matfileObj.stack_info.jammed;
            obj.mask = obj.matfileObj.stack_info.mask;
            obj.mask_vertices = obj.matfileObj.stack_info.mask_vertices;
        end
        % Getter for dependent properties
        function value = get.displacements(obj)
            if isfield(obj.cache, 'displacements')
                value = obj.cache.displacements;
            else
                value = obj.matfileObj.stack_info.displacements;
                obj.cache.displacements = value;
            end
        end
        function value = get.img_data(obj)
            if isfield(obj.cache, 'img_data')
                value = obj.cache.img_data;
            else
                value = obj.matfileObj.stack_info.img_data;
                obj.cache.img_data = value;
            end
        end
        function value = get.particle_locations(obj)
            if isfield(obj.cache, 'particle_locations')
                value = obj.cache.particle_locations;
            else
                value = obj.matfileObj.stack_info.particle_locations;
                obj.cache.particle_locations = value;
            end
        end
        function value = get.distances(obj)
            if isfield(obj.cache, 'distances')
                value = obj.cache.distances;
            else
                value = obj.matfileObj.stack_info.distances;
                obj.cache.distances = value;
            end
        end
        function value = get.gr(obj)
            if isfield(obj.cache, 'gr')
                value = obj.cache.gr;
            else
                value = obj.matfileObj.stack_info.gr;
                obj.cache.gr = value;
            end
        end
        function value = get.gr_bins(obj)
            if isfield(obj.cache, 'gr_bins')
                value = obj.cache.gr_bins;
            else
                value = obj.matfileObj.stack_info.gr_bins;
                obj.cache.gr_bins = value;
            end
        end
        function value = get.timestamps(obj)
            if isfield(obj.cache, 'timestamps')
                value = obj.cache.timestamps;
            else
                value = obj.matfileObj.stack_info.timestamps;
                obj.cache.timestamps = value;
            end
        end
        function value = get.psi_6(obj)
            if isfield(obj.cache, 'psi_6')
                value = obj.cache.psi_6;
            else
                value = obj.matfileObj.stack_info.psi_6;
                obj.cache.psi_6 = value;
            end
        end
        function value = get.bond_orders(obj)
            if isfield(obj.cache, 'bond_orders')
                value = obj.cache.bond_orders;
            else
                value = obj.matfileObj.stack_info.bond_orders;
                obj.cache.bond_orders = value;
            end
        end
        function value = get.voids(obj)
            if isfield(obj.cache, 'voids')
                value = obj.cache.voids;
            else
                value = obj.matfileObj.stack_info.voids;
                obj.cache.voids = value;
            end
        end

    end
end