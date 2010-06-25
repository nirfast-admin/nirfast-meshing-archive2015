function [e p] = run_mmc(img,nrow,ncol,nslice,xypixelsize,zpixelsize,edgesize,outputdir,outputfn)

mypwd = pwd;
cd(outputdir);


props.dim = 3;
props.mysize = [nrow ncol nslice];
props.fov = [[nrow ncol]*xypixelsize zpixelsize*nslice];
props.interval = [xypixelsize xypixelsize zpixelsize];
props.orient = 'axis';
props.datatype = 'WORD';
props.endian = 'ieee-le';

mmcfn = 'mmc_imgstack';
WriteSPRHeaderInfo(mmcfn,props);

fid = fopen([mmcfn '.sdt'],'wb');
for i=1:nslice
    count = fwrite(fid,(img(:,:,i))','int16');
    assert(count==nrow*ncol);
end



fclose(fid);
%% Run M3C
% set up the input file to M3C
fid = fopen('m3cinput.txt','wt');
fprintf(fid,'%s\n',outputdir);
fprintf(fid,'%s\n',mmcfn);
fprintf(fid,'0\n1\n1\n1\n');
% X-Y and Z spacing
xyspacing = edgesize/xypixelsize;
zspacing  = edgesize/ zpixelsize;
fprintf(fid,'%d\n%d\n',ceil(xyspacing),ceil(zspacing));
fprintf(fid,'%s\n',outputfn);
fclose(fid);

systemcommand=GetSystemCommand('m3c');
m3c_cmd=['! "' systemcommand '" < m3cinput.txt > junk.txt'];

fprintf('\n-----------> Running Marching Cube, please wait...\n\n');
eval(m3c_cmd);
fprintf('\n-----------> done.\n\n');

%% Read in the generated surface file and call checkerboard3d
[e p] = read_nod_elm([outputfn '_tetgen'],1);

%% Delete extra files
delfiles = {'tmp.elm','tmp.nod','mc_layer.elm','mc_layer.nod','layer.nod',...
             'm3cinput.txt',[mmcfn '.sdt'],[mmcfn '.spr'],'junk.txt',...
             [outputfn '_nod.dat'],[outputfn '_elm.dat'],...
             [outputfn '_tetgen.node'],[outputfn '_tetgen.ele']};
for i=1:length(delfiles)
    delete(delfiles{i});
end

cd(mypwd);