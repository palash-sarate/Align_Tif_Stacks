classdef Stackinfo < handle
    properties (Access = private)
        matfileObj % Handle to the matfile object
        matfilePath
        varList
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
        displacements
        img_data
    end
    properties (Dependent)
        % Dependent properties (loaded on demand)
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
            obj.matfilePath = matfilePath;
            if exist(obj.matfilePath, 'file') ~= 2
                error('The specified .mat file does not exist: %s', obj.matfilePath);
            end
            obj.matfileObj = matfile(obj.matfilePath, 'Writable', true); % Load as matfile object
            obj.varList = who(obj.matfileObj);
            obj.cache = struct();

            % Preload small fields into memory
            obj.start_index = loadVar(obj, 'start_index');
            obj.end_index = loadVar(obj, 'end_index');
            obj.parentDir = loadVar(obj, 'parentDir');
            obj.iteration = loadVar(obj, 'iteration');
            obj.aligned = loadVar(obj, 'aligned');
            obj.shortened = loadVar(obj, 'shortened');
            obj.masked = loadVar(obj, 'masked');
            obj.bd = loadVar(obj, 'bd');
            obj.number_density = loadVar(obj, 'number_density');
            obj.empty = loadVar(obj, 'empty');
            obj.nhood = loadVar(obj, 'nhood');
            obj.jammed = loadVar(obj, 'jammed');
            obj.mask = loadVar(obj, 'mask');
            obj.mask_vertices = loadVar(obj, 'mask_vertices');
            obj.displacements = loadVar(obj, 'displacements');
            obj.img_data = loadVar(obj, 'img_data');
            fprintf('loaded small fields into memory\n');
        end
        function value = loadVar(obj, varName)
            if ismember(varName, obj.varList)
                % tic; % Start timing
                temp = load(obj.matfilePath, varName);
                value = temp.(varName);
                % value = obj.matfileObj.(varName); % extremely SLOW
                % elapsedTime = toc; % Stop timing
                % fprintf('Time taken to load variable "%s": %.6f seconds\n', varName, elapsedTime);
            else
                value = [];
                warning('Variable "%s" does not exist in the matfile. Setting to null.', varName);
            end
            obj.cache.(varName) = value; % Cache the loaded variable
        end
        % Getter for dependent properties
        % function value = get.displacements(obj)
        %     fprintf('Loading displacements from matfile\n');
        %     if isfield(obj.cache, 'displacements')
        %         fprintf('Loading displacements from cache\n');
        %         value = obj.cache.displacements;
        %     else
        %         value = loadVar(obj, 'displacements');
        %     end
        % end
        % % Setter for displacements
        % function set.displacements(obj, value)
        %     obj.matfileObj.displacements = value; % Update matfile
        %     obj.cache.displacements = value; % Update cache
        % end
        % function value = get.img_data(obj)
        %     fprintf('Loading img_data from matfile\n');
        %     if isfield(obj.cache, 'img_data')
        %         fprintf('Loading img_data from cache\n');
        %         value = obj.cache.img_data;
        %     else
        %         value = loadVar(obj, 'img_data');
        %     end
        % end
        % % Setter for img_data
        % function set.img_data(obj, value)
        %     obj.matfileObj.img_data = value; % Update matfile
        %     obj.cache.img_data = value; % Update cache
        % end
        function value = get.particle_locations(obj)
            if isfield(obj.cache, 'particle_locations')
                value = obj.cache.particle_locations;
            else
                value = loadVar(obj, 'particle_locations');
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
                value = loadVar(obj, 'distances');
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
                value = loadVar(obj, 'gr');
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
                value = loadVar(obj, 'gr_bins');
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
                value = loadVar(obj, 'timestamps');
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
                value = loadVar(obj, 'psi_6');
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
                value = loadVar(obj, 'bond_orders');
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
                value = loadVar(obj, 'voids');
            end
        end
        function set.voids(obj, value)
            obj.matfileObj.voids = value; % Update matfile
            obj.cache.voids = value; % Update cache
        end

    end
end