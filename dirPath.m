clc;
clear;
ImgDirPath = 'D:\Car Side View\';
OutDirPath = 'D:\sun\OUT\';
ImgFiles = dir(strcat(ImgDirPath, '*.jpg'));
OutFiles = dir(strcat(OutDirPath, '*.jpg'));
length = size(ImgFiles,1);

for i = 1 : length
    ImgFiles(i).isdir = 0;
end

for j = 1 : size(OutFiles, 1)
    for k = 1 : length
        if strcmp(ImgFiles(k).name, OutFiles(j).name) == 1
            ImgFiles(k).isdir = 1;
        end
    end
end