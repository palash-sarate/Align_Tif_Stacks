wd = 'D:\HCA\DATA\Radius distribution';
before_folder = '\before';before_binarised_folder = '\before_binarised';
before_final_folder = '\before_final';before_onlylargest_folder = '\before_onlyLargest';
after_folder = '\after';after_binarised_folder = '\after_binarised';
after_final_folder = '\after_final';after_onlylargest_folder = '\after_onlyLargest';
%% Binarise
% binarise(wd,before_folder,before_binarised_folder,before_onlylargest_folder,before_final_folder);
% binarise(wd,after_folder,after_binarised_folder,after_onlylargest_folder,after_final_folder);
before_data = analyse(wd,before_final_folder);
% after_data = analyse(wd,after_final_folder);
%% find area of white
% close all;
% cd([wd,before_final_folder]);
% final_files = dir('**/*.tiff');
% whiteArea = [];
% for j = 1:numel(final_files)
%     file = [final_files(j).folder,'\',final_files(j).name];
%     img = Tiff(file,'r');
%     img_data = read(img);img_data = ~(~img_data);
%     
%     whiteArea(j) = bwarea(img_data);%nnz(img_data);%
%     
% end
% plot(whiteArea)


%% funcitons
%% find boundaries, center(centroid),radius vs theta,
%this was made into a different function instead of inserting in 2nd for
%loop of binarise funciton because few cases where final image is not as
%required were manually deleted
function table_ = analyse(wd,final_folder)
    cd([wd,final_folder]);
    files = dir('**/*.tif');
    area = zeros(numel(files),1);
    center = zeros(numel(files),2);
    radius = zeros(numel(files),2);
    rounded_ness = zeros(numel(files),1);
    for i = 1:numel(files)
        file = [files(i).folder,'\',files(i).name];
        img = Tiff(file,'r');
        img_data = read(img);
        img_data = uint8(img_data)/255;
        img_data = ~(~img_data);
        stats = regionprops('table',img_data,'Centroid',...
            'MajorAxisLength','MinorAxisLength');
%         disp(stats.Centroid)
        current_area = bwarea(img_data);
        area(i) = current_area;
        center(i,:) = stats.Centroid;
        radius(i,:) = [stats.MajorAxisLength/2 stats.MinorAxisLength/2];
        rounded_ness(i) = roundedness(current_area,stats.MajorAxisLength/2);
%         boundary = bwboundaries(img_data)
    end
    % imshow(img_data)
    % hold on;
    % viscircles(centers,radii);
    % hold off
    table_ = table(area,center,radius,rounded_ness);
    writetable(table_,'data.csv');
end

function r_ness = roundedness(a,rmax)
%sphericity is the ratio of the volume of particle and volume of sphere
%having equal surface area as the particle
% Roundness is defined as the ratio of the surface area of an object to the
% area of the circle whose diameter is equal to the maximum diameter of the object
    r_ness = a/(pi*(rmax^2));
end

function binarise(wd,folder,binarised_folder,onlyLargest_folder,final_folder)
%Binarise -> find largest area and fill all other areas ->
%erode-fill-dilate to remove connected areas
cd([wd,folder]);
files = dir('**/*.tif');
    for i = 1:numel(files)
        file = [files(i).folder,'\',files(i).name];
        img = Tiff(file,'r');
        img_data = read(img);
        img_data = uint8(img_data)/255;
        img_data = ~img_data;
        img_data = imfill(img_data,'holes');
        n = 1;
        padded = padarray(img_data,[n n],0,'both'); % pad
        largestBlob = bwareafilt(padded, 1);
        onlyLargest = padded & largestBlob;
        onlyLargest = onlyLargest(n+1:end-n,n+1:end-n); % unpad
    %     imshow(onlyLargest);
        imwrite(img_data,[wd,binarised_folder,'\',files(i).name]);
        imwrite(onlyLargest*255,[wd,onlyLargest_folder,'\',files(i).name]);
    end
cd([wd,onlyLargest_folder]);
binarised_files = dir('**/*.tif');
    for j = 1:numel(binarised_files)%31%
        file = [binarised_files(j).folder,'\',binarised_files(j).name];
        img = Tiff(file,'r');
        img_data = read(img);img_data = ~(~img_data);

    %     img_data(1:end,1)

        se = strel('disk',10);
        m=10;
        eroded = img_data;
        for i = 1:m
            eroded = imerode(eroded,se);
        end
            largestBlob = bwareafilt(eroded, 1);
        %     onlyLargest = eroded & largestBlob;
        dilated = largestBlob;
        for i = 1:m
            dilated = imdilate(dilated,se);
        end
    %     imshowpair(eroded,dilated,'montage');
        imwrite(dilated*255,[wd,final_folder,'\',binarised_files(j).name]);
    end
end