function [t p] = readVTK(fn)

fn = add_extension(fn,'.vtk');

fid = OpenFile(fn,'rt');

s = fgetl(fid);
if isempty(regexp(s,'Version 3\.0', 'once'))
    error('Can only read version 2.0 if vtk file format.')
end
s = fgetl(fid);
s = fgetl(fid);

if ~strcmp(s,'ASCII')
    error('Can only read ASCII form of vtk files.')
end

s = fgetl(fid);
tmp = sscanf(s,'%*s%s');
if ~strcmp(tmp,'POLYDATA')
    error(['Can not read dataset of type: ' tmp])
end

s = fgetl(fid);
keyword = sscanf(s,'%s',1);
if ~strcmp(keyword,'POINTS')
    error('Expecting POINTS data in vtk file')
end
nn = str2double(sscanf(s,'%*s%s'));
p=zeros(nn*3,1);

s = fgetl(fid);
sidx = 1;

while ~strcmp(sscanf(s,'%s',1),'POLYGONS')
    tmp = textscan(s,'%f');
    tmp=tmp{1};
    eidx = sidx + length(tmp) - 1;
    p(sidx:eidx) = tmp;
    sidx = eidx + 1;
    s = fgetl(fid);
    while isempty(s)
        s = fgetl(fid);
    end
end
p = (reshape(p,3,nn))';

ne = str2double(sscanf(s,'%*s%s'));

data = textscan(fid,'%d %d %d %d');
if length(data{1}) ~= ne
    error('Can not parse the vtk file.')
end
t = [data{2} data{3} data{4}];
t=t+1;
fclose(fid);





