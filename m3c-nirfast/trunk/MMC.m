function MMC(maskloc, xypixelsize, zpixelsize, edgesize, inploc)
% Use either a set of bmp files or a 3D matlab matrix (mask) to create a
% multiple material surface by using Ziji's m3c algorithm
% 

if ischar(maskloc)
    mask = GetBMPStack(maskloc);
else
    mask = maskloc;
end

[outputdir, outputfn, ext] = fileparts(inploc);

if isempty(outputfn)
    outputfn='m3c-surface';
end
if isempty(outputdir)
    outputdir=pwd;
end

[nrow ncol nslice]=size(mask);
[t p]=run_mmc(mask,nrow,ncol,nslice,xypixelsize,zpixelsize,edgesize,outputdir,outputfn);

writenodelm_nod_elm([outputdir filesep outputfn],t,p,[],2);





function mask = GetBMPStack(maskloc)

foo = dir([maskloc '*.bmp']);
if isempty(foo)
    errordlg('Can not find BMP files','Meshing Error');
    error('Can not find BMP files');
end

a = imread([maskloc '1.bmp']);
if ndims(a)==3
    a=rgb2gray(a);
end
[nrow ncol]=size(a);
mask = zeros(nrow,ncol,length(foo),'int16');

for i=1:length(foo)
    a = imread([maskloc num2str(i) '.bmp']);
    if ndims(a)==3
        a=rgb2gray(a);
    end
    mask(:,:,i) = a;
end
