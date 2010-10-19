function [vol,q,q_area,status]=CheckMesh3D(e,p,nodemap,input_flags)
% CheckMesh3D(type,e,p,nodemap) checks the integrity of a 3D mesh
% User can provide list of elements and points or provide a filename that
% contains the mesh. If 'type' is zero then he will be prompted for a
% filename in form of .nod .elm, provided that 'fn' is empty
% The format of .nod header is:
% <# of points> <dimension> <# of attributes> <# of boundary markers (0 or 1)>
% # Remaining lines list # of points:
% <point #> <x> <y> <z> [attributes] [boundary marker]
% The format for .elm is:
% <# of elements> <nodes per element> <# of attributes>
% # Remaining lines list of # of elements:
% <elm #> <node> <node> <node> <node> ... [attributes]
% type:
% 0: do all diagnostics on a mesh read from file
% 1: do all the diagnostics (edge count, void element, quality)
% 2: do a mesh quality only
% 3: only check if there is any void element
% Return valuse: (status is a 'struct')
% 1: Provided surface is multiple material
% 2: Provided closed surface is not valid! 
% 3: surface is both multip material and open
% 4: Provided surface has faces with low quality < q_area_threshold

if nargin>3 && ~isempty(input_flags)
    if isfield(input_flags,'type')
        type = input_flags.type;
    else
        type=0;
    end
    if isfield(input_flags,'verbose')
        verbose=input_flags.verbose;
    else
        verbose=0;
    end
elseif nargin<3
    type = 1;
    verbose = 1;
    input_flags=[];
end
vol=0;
q=0;
q_area=0;
edgeflag=0;
q_area_threshold=0.1;
status.a=0; status.b=0; status.c=0;
global TetrahedronFailQuality
if isempty(TetrahedronFailQuality)
    TetrahedronFailQuality=0.03;
end
if type~=0
    fn='Diagnostic-3DMesh';
end
if type==0
    [e,p,nodemap,elemap,dim,nnpe,fn]=ui_read_nod_elm(0);
    if dim~=3
        disp('Wrong type of mesh. It is not a 3D mesh!');
        return;
    end
else
    np=size(p,1);
    ne=size(e,1);
    nnpe = size(e,2);
end
if type==0
    if nnpe==3
        if verbose, disp([fn ' should be a surface mesh.']); end
        [q,q_area,area,zeroflag,edges,edgeflag]=checkmesh3d_surface(e,p,nodemap);
    else
        if verbose, disp([fn ' should be a solid mesh.']); end
        [vol,q,zeroflag]=checkmesh3d_solid(e,p,1,nodemap);
    end
else
    nnpe=size(e,2);
    if nnpe==3
        if verbose, disp('Checking surface mesh:');end
        if nargin==4 && ~isempty(nodemap)
            [q,q_area,area,zeroflag,edges,edgeflag]=checkmesh3d_surface(e,p,nodemap);
        else
            [q,q_area,area,zeroflag,edges,edgeflag]=checkmesh3d_surface(e,p);
        end
    else
        if verbose, disp('Checking solid mesh:');end
        if nargin==4
            [vol,q,zeroflag]=checkmesh3d_solid(e,p,1,nodemap);
        else
            [vol,q,zeroflag]=checkmesh3d_solid(e,p,1);
        end
    end
end
if nnpe==3
    if edgeflag~=0
        switch edgeflag
            case 1
                status.b=1;
                if verbose, 
                    fprintf('\n Some of mesh edges are shared by more than two triangles! (multiple material mesh ?!)\n');
                    fprintf(' Please check all the diagnostic results stored in %s.txt files.\n',fn);
                end
            case 2
                status.b=2;
                if verbose, 
                    fprintf('\n  Provided surface is not closed:\n');
                    fprintf('  At least one of the edges is only shared by only one triagnle (it should be two, at least)\n');
                end
            otherwise
                if verbose, fprintf('\n Surface mesh is open AND has edges shared by more than 2 triangles!\n');end
                status.b=3;
        end
    end
    qcheck=q_area<q_area_threshold;
    a=sum(qcheck);
    if a~=0
        if verbose, fprintf(' There are %d faces with low quality (q<%f).\n', a, q_area_threshold);end
        status.b=4;
    end
elseif nnpe==4
    if zeroflag~=0
    end
    qcheck=q<TetrahedronFailQuality;
    a=sum(qcheck);
    if a~=0
        if verbose, fprintf('There are %d tetrahedrons that have a very low quality (q<%f).\n', a, TetrahedronFailQuality);end
        status.b=4;
    end
end

if isfield(input_flags,'writefiles') && input_flags.writefiles==1
    fn1=[fn '-quality-radius.txt'];
    fn2=[fn '-area.txt'];
    fn3=[fn '-edges.txt'];
    fn4=[fn '-edge-conn.txt'];
    fn5=[fn '-volume.txt'];
    os=computer;
    if ~isempty(strfind(os,'PCWIN')) % Windows
        newlinech ='pc';
    elseif ~isempty(strfind(os,'MAC')) ||  ~isempty(strfind(os,'GLNX86')) % Mac OS or Linux
        newlinech ='unix';
    end
    fprintf('\n');
    if ~isempty(q) && nnpe==3
        fprintf(' Avg Quality (radius ratio): %f\n', mean(q));
        fprintf(' Avg Quality (area   ratio): %f\n', mean(q_area));
        dlmwrite(fn1,q,'delimiter',' ','newline',newlinech);
        fn1=[fn '-quality-area.txt'];
        dlmwrite(fn1,q_area,'delimiter',' ','newline',newlinech);
    end
    if nnpe==3 && ~isempty(area)
        fprintf(' Avg Area: %f\n', mean(area));
        dlmwrite(fn2,area,'delimiter',' ','newline',newlinech);
    end
    if nnpe==3 && ~isempty(edges)
        dlmwrite(fn3,edges,'delimiter',' ','newline',newlinech);
    end
    if nnpe==4 && ~isempty(vol)
        [v,idx]=sort(vol);
        dlmwrite(fn5,[idx v],'delimiter',' ','newline',newlinech);
    end
    if nnpe==4 && ~isempty(q)
        fprintf(' Avg Quality (vol ratio): %f\n', mean(q));
        dlmwrite([fn 'quality-vol-ratio.txt'],q,'delimiter',' ','newline',newlinech);
    end
end
    













