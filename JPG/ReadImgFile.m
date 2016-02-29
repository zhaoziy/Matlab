function [PointMatrix, GroundLine]  = ReadImgFile(Path)
File = fopen(Path,'r+');% 需要读取的文件
i = 0;
PointMatrix = [];
GroundLine =[];
while i<100
    tline=fgetl(File);%读取一行
    i=i+1;
    if mod(i,4) ~= 1
        S = regexp(tline, '\s+', 'split');
        n1 = str2double(S(2));
        n2 = str2double(S(3));
        PointMatrix = [PointMatrix; n1 n2];
    else
        S = regexp(tline, '\s+', 'split');
        n1 = str2double(S(5));
        n2 = str2double(S(6));
        PointMatrix = [PointMatrix; n1 n2];
    end
end

while i<168
    tline=fgetl(File);%读取一行
    i=i+1;
    if i > 164
        if mod(i,4) ~= 1
            S = regexp(tline, '\s+', 'split');
            n1 = str2double(S(2));
            n2 = str2double(S(3));
            GroundLine = [GroundLine; n1 n2];
        else
            S = regexp(tline, '\s+', 'split');
            n1 = str2double(S(5));
            n2 = str2double(S(6));
            GroundLine = [GroundLine; n1 n2];
        end
    end
end

i=0;
fclose(File);