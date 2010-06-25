function [t p] = MMC(maskloc, xypixelsize, zpixelsize, edgesize, inploc)
% Use either a set of bmp files or a 3D matlab matrix (mask) to create a
% multiple material surface by using Ziji's m3c algorithm
% 
% Output:
% t : is a n by 5 matrix (first 3 are nodes and the rest are material properties.
% p : is the list of coordinates (m by 3)
% 
% The generated surface will also be saved in a pair of files as tetgen format:
% 'inploc.node' and 'inploc.ele'

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
writenodelm_surface_medit([outputdir filesep outputfn '.mesh'],t(:,1:3),p(:,1:3));

if nargout == 0
    clear t p
end





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
