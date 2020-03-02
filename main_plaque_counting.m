%This script calls iteratively the analysis for single tif files and save the results in a csv file

%%%%%%%%%%%%%%%%%%%%%%%%%% SETTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if method == 0 use global Otsu, otherwise 1 use adaptive threshold of Bradley and Roth.
method = 1;
% If you want to use the active contours segmentation to remove artefacts at the border of the pictures (as signal in the meningis)
% '1' means  use active contours, '0' don't. Consider that depending on the size of the volume this can seriously slow down the analysis
use_snake = 0;
load_partially = 1; % 0 to load one volume intererly, 1 to load parts

%%%%%%%%%%%%%%%%%%%%%%%%%% LOADING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If all files in one folder use the list
%fname ='Q:\Data\Stitching_AnnaMaria\Stitched\Full_resolution_515LP_filter\FullRes_em515LP__2x__left_x2064_y3879_z1353.tif';
%names = {fname};

working_dir = 'Y:\data\a.reuss\200108_MB_APPPS1_Tg_6mo_vessels\output\test\';

files = dir(fullfile(working_dir, '*.tif'));
filenames = arrayfun(@(x) x.name, files, 'UniformOutput', false);
filenames = filenames';
filepaths = arrayfun(@(x) x.folder, files, 'UniformOutput', false);
filepaths = filepaths';
names = fullfile(filepaths, filenames);
 
% Go to this dir to save data (and where the code is)
%cd('C:/Users/Administrator/Documents/cell_tools');


%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CALL %%%%%%%%%%%%%%%%%

res = zeros( length(names),1);
vol_plaque = zeros( length(names),1);
convex_hull = zeros( length(names),1);

% Split the big volume in n pieces of specific length
tot_block = 4;

for kk = 1 : length(names)
    disp(names(kk))
    
    if (load_partially)
        [data,num_images] = imread_big( names{kk});
        length_block = floor(num_images/tot_block);
        for zz = 1 : tot_block
            subdata =  data(:,:, 1 + ( (zz-1)*length_block )   :   length_block + ( (zz-1)*length_block ) );
            [res_t, dens_z, vol_t] = extract_feat_par( names(kk), subdata, method, use_snake,zz,tot_block);
                                                    
            temp_res(zz) = res_t;
            temp_vol(zz) = vol_t;
            temp_dens(zz) = dens_z;
        end
        res(kk) = sum(temp_res);
        vol_plaque(kk) = sum(temp_vol);
        dens(kk) = sum(temp_dens);
    else
        [res(kk),dens(kk),vol_plaque(kk)] = extract_feat_big( names{kk} , method, use_snake);
    end
    
    filename = [working_dir  'results.csv'];
    T = table( names',res,vol_plaque,dens);
    T.Properties.VariableNames={'Filename', 'CellCount','Volume','Density'};
    writetable(T,filename,'Delimiter',',');
end

%Initalize all labels as Placebo
labels = repmat('Placebo',length(names),1);
for uu = 1 : length(names)
         k = strfind(names{uu},'high');
         if (~isempty(k))
            labels(uu) = 'HD';
         end
         k = strfind(names{uu},'low');
         if (~isempty(k))
            labels(uu) = 'LD';
         end    
end

save %This will save variables into a .mat file
%split results
data_placebo = [];
data_hd = [];
data_ld = [];
data_placebo_vol = [];
data_hd_vol = [];
data_ld_vol = [];
for uu = 1 : length(names)
         k = strfind(names{uu},'high');
         if (~isempty(k))
             data_hd(end+1) = res(uu);
             data_hd_vol(end+1) = vol(uu);
         end
         k = strfind(names{uu},'low');
         if (~isempty(k))
            data_ld(end+1) = res(uu);
            data_ld_vol(end+1) = vol(uu);
         end    
         k = strfind(names{uu},'Placebo');
         if (~isempty(k))
            data_placebo(end+1)  = res(uu);
            data_placebo_vol(end+1)  = vol(uu);
         end    
end
hold on
boxplot(res,labels)
scatter(ones(size(data_hd))  ,data_hd,'r','filled')
scatter(ones(size(data_placebo))*3,data_placebo,'r','filled')
scatter(ones(size(data_ld))*2,data_ld,'r','filled')
hold off
hold on
figure; boxplot(vol_plaque,res)
scatter(ones(size(data_hd))  ,data_hd,'r','filled')
scatter(ones(size(data_placebo_vol))*3,data_placebo_vol,'r','filled')
scatter(ones(size(data_ld_vol))*2,data_ld_vol,'r','filled')
hold off

% This will export everything into a csv file
filename = 'dati_wren_90days.csv'; 
T = table( names',res,vol_plaque);
T.Properties.VariableNames={'Filename', 'CellCount','Volume'};
writetable(T,filename,'Delimiter',',');
