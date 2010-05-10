function [tc,nodes,status]=FixPatchOrientation(p,sfaces,meshnodemap,silentflag)
% [tc,nodes,status]=FixPatchOrientation(p,sfaces,meshnodemap)
% Assuming the 3D shell defined by p and sfaces is a closed one, this
% routine will make sure  that all the facets are pointing outward.
global logfilename
% st=CheckEulerFormula(sfaces);
% if isempty(logfilename), logfilename='FixPOrient.log';end
% if ~st
%     logmessage(logfilename,'FPO: Input surface should be closed Euler surface!',1,1);
% end 
status=0;
if nargin<=3 || silentflag==0
    fprintf('%s\n','===========================================')
    fprintf('%s\n','Entered FixPatchOrientation() function')
end

nodes=[sfaces(:,1);sfaces(:,2);sfaces(:,3)];
nodes=unique(nodes);
if nargin==2 || (nargin>=3 && isempty(meshnodemap))
    ren_nodes=nodes;
elseif nargin>=3 && ~isempty(meshnodemap)
    [tflag1 ren_nodes]=ismember(nodes,meshnodemap);
    MyAssert(sum(tflag1)==size(nodes,1),'Error in meshnodemap!');
end
pp=p(ren_nodes,:);
[tf ee]=ismember(sfaces,nodes);
% writenodelm_2dm('fix_facets_orientation_inp_temp.2dm',ee,pp);
% 
ee=double(ee);
list=GetListOfConnTri2Tri_mex(ee,pp); 
[status tc1]=orient_surface_mex(ee,pp,list);
if status~=0
%     MyAssert(status==0,['Could not orient the given surface. Error Code: ' num2str(status)]);
    tc=sfaces;
    return
else
    tc=nodes(tc1);
    new_nodes=unique([tc(:,1);tc(:,2);tc(:,3)]);
    MyAssert(length(new_nodes)==length(nodes),sprintf(' FixPatchOrientation: sfchk.exe has altered the number of input points.\n Probably your input surface had a dangling triangle/edge/node!'))
    if nargin<=3 || silentflag==0
        fprintf('%s\n','Exitting FixPatchOrientation() function.')
        fprintf('%s\n','===========================================')
    end
end


% % Prepare an input text file for sfck
% fid=fopen('sfchk.txt','wt');
% fprintf(fid,'fix_facets_orientation_inp_temp.2dm\n');
% fprintf(fid,'fix_facets_orientation_out_temp.2dm\n');
% fclose(fid);

% os=computer;
% if ~isempty(strfind(os,'PCWIN')) % Windows
%     systemcommand = 'sfchk.exe < sfchk.txt';
% elseif ~isempty(strfind(os,'MAC')) % Mac OS
%     systemcommand = './sfchk-mac < sfchk.txt';
% elseif ~isempty(strfind(os,'GLNX86')) % Linux
%     systemcommand = './sfchk-linux < sfchk.txt';
% end
% 
% disp(sprintf('\t%s','Executing sfchk...'))
% [status,msg]=system(systemcommand);
% disp(sprintf('\b%s\n',' done.'))
% 
% if status==1
%     error('Can not run sfchk.exe!!!');
% end
% [tc,pc]=read2dm_3d('fix_facets_orientation_out_temp.2dm');
% tc=uint32(tc);

% xx=tc1-tc;
% MyAssert(sum(xx(:))==0,'Different Orientation from
% orient_surface_mex()!');
% tc=[tc(:,1) tc(:,2) tc(:,3)];
