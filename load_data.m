
% Source and Data Paths
% currently searches for data and source folder within the repo that are
% gitignore. You can add the folders or change source and data directories.
src_dir = '.\src\'
data_dir = '.\data\'

% The number of images loaded - 1500 is the maximum number of images
img_number = 1500;
% img_number = 300;

% make the functions in the source avaialble to use
addpath(strcat(src_dir,'NIfTI_20140122\'));
addpath(strcat(src_dir,'display\'));
addpath(strcat(src_dir,'transforms\'));

% get data path
image_path = strcat(data_dir,'images\');
registration_path = strcat(data_dir,'registrations\');
segmentation_path = strcat(data_dir,'segmentation\');

% get list of all files in the folders
image_listing           = dir(image_path);
registration_listing    = dir(registration_path);
segmentation_listing    = dir(segmentation_path);

% load initial nifti structure
images = load_untouch_nii(strcat(image_path,'0001.nii.gz'));
cpg1 = load_untouch_nii(strcat(registration_path, '0001_cpp_region1.nii.gz'));
cpg2 = load_untouch_nii(strcat(registration_path, '0001_cpp_region2.nii.gz'));
segmentations = load_untouch_nii(strcat(segmentation_path,'0007_mask.nii.gz'));

disp('loading images')

% load images
count1 = 1;
for i=1:size(image_listing,1) 
    if contains(image_listing(i).name,'nii')

        if count1 > img_number
            break;
        end
        
        filename = strcat(image_path, image_listing(i).name);
        images(count1) = load_untouch_nii(filename);
        count1 = count1 + 1;
    end
end

disp('loading registrations')

%load registration
count1 = 1;
count2 = 1;
for i=1:size(registration_listing,1) 
    if contains(registration_listing(i).name,'nii')
        
        if count2 > img_number
            break;
        end
        
        filename = strcat(registration_path, registration_listing(i).name);      
        if contains(registration_listing(i).name, 'region1')
            cpg1(count1) = load_untouch_nii(filename);
            count1 = count1 + 1;
        else
            cpg2(count2) = load_untouch_nii(filename);
            count2 = count2 + 1;
        end  
    end
end

disp('loading segmentations')

% load masks
count1 = 1;
for i=1:size(segmentation_listing,1) 
    if contains(segmentation_listing(i).name,'nii')
        filename = strcat(segmentation_path, segmentation_listing(i).name);
        segmentations(count1) = load_untouch_nii(filename);
        count1 = count1+1;
    end
end

addpath(strcat(data_dir,'images\'));
addpath(strcat(data_dir,'registrations\'));
addpath(strcat(data_dir,'segmentation\'));

disp('loading data completed')