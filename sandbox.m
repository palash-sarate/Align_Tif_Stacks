clc;
% list all tif files > 5mb
files = dir('E:\shake_table_data\**\*.tif');
upper_limit = 20 * 1024 * 1024;
for k = 1:length(files)
    if files(k).bytes > upper_limit
        % confirm that files(k).name only has single number in it
        if sum(isstrprop(files(k).name, 'digit')) == 1
            sprintf("%s - %s",files(k).name, files(k).folder)
            % delete file
            delete(fullfile(files(k).folder, files(k).name))
        end
    end
end

files = dir('E:\shake_table_data\**\*.avi');
upper_limit = 20 * 1024 * 1024;
for k = 1:length(files)
    if files(k).bytes > upper_limit
        % confirm that files(k).name only has single number in it
        if sum(isstrprop(files(k).name, 'digit')) == 1
            sprintf("%s - %s",files(k).name, files(k).folder)
            % delete file
            delete(fullfile(files(k).folder, files(k).name))
        end
    end
end