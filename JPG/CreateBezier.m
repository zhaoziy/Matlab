function [x,y]  = CreateBezier(PointMatrix)
length = size(PointMatrix,1)/4;
x = [];
y = [];
for i =1:length
    p0 = PointMatrix((i-1)*4 + 1, 1);
    p1 = PointMatrix((i-1)*4 + 2, 1);
    p2 = PointMatrix((i-1)*4 + 3, 1);
    p3 = PointMatrix((i-1)*4 + 4, 1);
    
    p0 = [p0 PointMatrix((i-1)*4 + 1, 2)];
    p1 = [p1 PointMatrix((i-1)*4 + 2, 2)];
    p2 = [p2 PointMatrix((i-1)*4 + 3, 2)];
    p3 = [p3 PointMatrix((i-1)*4 + 4, 2)];
    
    t=0:0.1:1;
    x=[x (1-t).^3*p0(1)+3*t.*(1-t).^2*p1(1)+3*t.^2.*(1-t)*p2(1)+t.^3*p3(1)];
    y=[y (1-t).^3*p0(2)+3*t.*(1-t).^2*p1(2)+3*t.^2.*(1-t)*p2(2)+t.^3*p3(2)];
end
end
