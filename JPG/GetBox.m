function [left, right, top, bottom]  = GetBox(xv, GroundLine)
bottom = GroundLine(1,2);
for i= 2:4
    if GroundLine(i,2) > bottom
        bottom = GroundLine(i,2);
    end
end

left = xv(1);
right = xv(1);
length = size(xv);
for i= 2:length(2)
    if xv(i) < left
        left = xv(i);
    end
    if xv(i) > right
        right = xv(i);
    end
end

left = left - 3;
right = right + 3;

top = bottom - (right - left) * 3 / 7;
end
