clc;
clear;
ImgDirPath = 'd:\\';
TxtDirPath = 'd:\\';
OutDirPath = 'd:\\OUT\\';
OutResizeDirPath = 'd:\\OUTResize\\';
ImgFiles = dir(strcat(ImgDirPath, '*.jpg'));
LengthFiles = length(ImgFiles);
for k = 1:LengthFiles
    [pathstr, TxtName, ext] = fileparts(strcat(ImgDirPath,ImgFiles(k).name));
    [PointMatrix, GroundLine] = ReadImgFile(strcat(TxtDirPath, TxtName, '.TMPLT'));
    [xv, yv]  = CreateBezier(PointMatrix);
    [left, right, top, bottom]  = GetBox(xv, GroundLine);
    img = imread(strcat(ImgDirPath,ImgFiles(k).name));
    width = size(img,2);
    height = size(img,1);
    for i= 1:width
        for j = 1:height
            if i > left && j > top && i < right && j < bottom
                if inpolygon(i,j,xv,yv) ~= 1
                    img(j,i,1) = 255;
                    img(j,i,2) = 255;
                    img(j,i,3) = 255;
                end
            else
                img(j,i,1) = 255;
                img(j,i,2) = 255;
                img(j,i,3) = 255;
            end
        end
        if mod(i,100) == 0
            [i, width, k]
        end
    end
    RGB = imcrop(img,[left, top, int32(right - left), int32(bottom - top)]);
    imwrite(img, strcat(OutDirPath,ImgFiles(k).name),'jpg');
    imwrite(RGB, strcat(OutResizeDirPath,ImgFiles(k).name),'jpg');
    clear img width height i j PointMatrix xv yv RGB;
end