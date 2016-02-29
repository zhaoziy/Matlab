function [xv yv]  = GetVec(PointMatrix)
length =size(PointMatrix,1);
xv = [];
yv = [];
for i = 1:length
    xv = [xv PointMatrix(i,1)];
    yv = [yv PointMatrix(i,2)];
end
