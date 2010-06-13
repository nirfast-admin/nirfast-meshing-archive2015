function [P]=ExpandBoundaryBufferZone(t,p,P,shell_normals,density,DELTA,llc)
% Expands each triangular facet of the boundary by a factor (depending on
% local mesh density) and then extrudes a prism out of the expanded
% triangle. Then checks for all the pixels that are within this prism and
% tags them as buffer pixels. This way we can make sure that there is no
% 'hole' in the boundary buffer.
global tiny NA boundary_node_code 
if isempty(tiny), tiny=1e-9; end

fprintf('%s\n','===========================================')

% Mesh template properties
[nrow ncol npln]=size(P);
dx=DELTA(1); dy=DELTA(2); dz=DELTA(3);
xmin=llc(1); ymin=llc(2); zmin=llc(3);

nf = size(t,1);
sqroot2=sqrt(2)/1.4;
bethafactor=sqroot2;

tmp=load('prism_topology');
prism=tmp.prism;
list=tmp.list;
indexing=tmp.indexing;
clear tmp;

cprintf([0 0 1],'\tCalculating desired length at boundary nodes...')
tp1=p(t(:,1),:); tp2=p(t(:,2),:); tp3=p(t(:,3),:);
[pc]= incircle(tp1,tp2,tp3);
v=[tp1 tp2 tp3]-repmat(pc,1,3);
L=[v_magn(v(:,1:3)) v_magn(v(:,4:6)) v_magn(v(:,7:9))];
sf=ones(nf,3,'single')*density;
% % Or use this one when Size function is available
% for i=1:nf
%     for k=1:3
%         sf(i,j)=GetEdgeSizeAtXY(p(t(i,j),:,density);
%     end
% end
betha=(L+bethafactor*sf)./L;
clear L sf;
tp1=v(:,1:3).*repmat(betha(:,1),1,3)+pc;
tp2=v(:,4:6).*repmat(betha(:,2),1,3)+pc;
tp3=v(:,7:9).*repmat(betha(:,3),1,3)+pc;
clear v pc;
fprintf('\b%s\n',' done.')

% Calculate each face's prism coordinates
pp=zeros(nf,3,6);
pp(:,:,1)=tp1+shell_normals*bethafactor*density;
pp(:,:,2)=tp2+shell_normals*bethafactor*density;
pp(:,:,3)=tp3+shell_normals*bethafactor*density;
pp(:,:,4)=tp1-shell_normals*bethafactor*density;
pp(:,:,5)=tp2-shell_normals*bethafactor*density;
pp(:,:,6)=tp3-shell_normals*bethafactor*density;
% Perturb the above nodes
for i=1:6
    pp(:,:,i)=RandomizeBoundaryPoints(pp(:,:,i),llc,DELTA);
end
    
% Get bbx of prisms' facets
prism_facets_bbx=zeros(nf,8,6,'single');
prism_normals=zeros(nf,3,8,'double');
cprintf([0 0 1],'\tCalculating prism normals and bounding boxes...')
for i=1:nf
    tpp=(reshape(pp(i,:,:),3,6))';
    n1=tpp(prism(:,1),:); n2=tpp(prism(:,2),:); n3=tpp(prism(:,3),:);
    prism_facets_bbx(i,:,1)=min([n1(:,1) n2(:,1) n3(:,1)],[],2);
    prism_facets_bbx(i,:,2)=min([n1(:,2) n2(:,2) n3(:,2)],[],2);
    prism_facets_bbx(i,:,3)=min([n1(:,3) n2(:,3) n3(:,3)],[],2);
    prism_facets_bbx(i,:,4)=max([n1(:,1) n2(:,1) n3(:,1)],[],2);
    prism_facets_bbx(i,:,5)=max([n1(:,2) n2(:,2) n3(:,2)],[],2);
    prism_facets_bbx(i,:,6)=max([n1(:,3) n2(:,3) n3(:,3)],[],2);
    v1=tpp(prism(:,2),:)-tpp(prism(:,1),:);
    v2=tpp(prism(:,3),:)-tpp(prism(:,2),:);
    tmp=cross(v1,v2);
    norm_len=v_magn(tmp);
    prism_normals(i,:,:)=(tmp./repmat(norm_len,1,3))';
end
fprintf('\b%s\n',' done.')
% Get bbx of prisms
prisms_bbx=zeros(nf,6);
for i=1:3
    prisms_bbx(:,i)   = min(pp(:,i,:),[],3);
    prisms_bbx(:,i+3) = max(pp(:,i,:),[],3);
end
mydir=[1 0 0];

% First we tag all the pixels around the boundary nodes based on desired
% length
np=size(p,1);
for i=1:np
    [I J K]=mapxyz2ijk(p(i,1:3),[dx dy dz],[nrow ncol npln],[xmin ymin zmin]);
    if P(I,J,K)==boundary_node_code % already tagged as boundary
        continue;
    end
    dd = density; 
    % Or use this one when density function is available
    % dd = GetDensityAtXY(p(i,1:3),density);
    dI = round(dd/dx) - 1;
    dI = round(dI/2); % To address node jumps at boundary
    istart=max(1,I-dI); iend=min(nrow,I+dI);
    for ii=istart:iend
        jstart=max(1,J-dI); jend=min(ncol,J+dI);
        for jj=jstart:jend
            kstart=max(1,K-dI); kend=min(npln,K+dI);
            for kk=kstart:kend
                if P(ii,jj,kk)==0
                    P(ii,jj,kk)=NA;
                end
            end
        end
    end
    P(I,J,K) = boundary_node_code;
end


% Get the maximum edge length of each triangle
edges=[t(:,[1 2]); t(:,[1 3]); t(:,[2 3])];
[edges m n] = unique(sort(edges,2),'rows');
edgelength=sqrt(sum((p(edges(:,2),:)-p(edges(:,1),:)).^2,2));


cprintf([0 0 1],'\tSealing boundary buffer zone...')
h = waitbar(0,'Sealing boundary buffer zone. Please wait!') ;
div=200;
ff=false;
totaltime=cputime;
total=nf/div;

t2=0;t4=0;
infX=10*abs(max(p(:,1))-min(p(:,1)))+max(p(:,1));

for i=1:nf
    if mod(i,div)==0
        total=total-1;
        minutes=fix(total*(cputime-totaltime)/60);
        seconds=round(mod(total*totaltime,60));
        waitbar(i/nf,h,['Sealing boundary buffer zone. Please wait! Left: '...
                        sprintf('%u:%2.0f min',minutes,seconds)]);
        totaltime=cputime;
    end
%     i
    % Check to see if half of the maximum edge length of the current face
    % is less the the desired density/distance
    eidx=[i i+nf i+2*nf];
    if max(edgelength(n(eidx)))/sqroot2 < density
        continue
    end
    [imin jmin kmin]=mapxyz2ijk(prisms_bbx(i,1:3),[dx dy dz],[nrow ncol npln],[xmin ymin zmin]);
    [imax jmax kmax]=mapxyz2ijk(prisms_bbx(i,4:6),[dx dy dz],[nrow ncol npln],[xmin ymin zmin]);
    istart=min(imin,imax); iend=max(imin,imax); % iend=istart+dI; % iend=max(imin,imax);
    jstart=min(jmin,jmax); jend=max(jmin,jmax); % jend=jstart+dJ; % jend=max(jmin,jmax);
    kstart=min(kmin,kmax); kend=max(kmin,kmax); % kend=kstart+dK; % kend=max(kmin,kmax);
    
    istart=max(istart,1); jstart=max(jstart,1); kstart=max(kstart,1);
    iend=min(iend,nrow); jend=min(jend,ncol); kend=min(kend,npln);

    sub_P=zeros(iend-istart+1,jend-jstart+1,kend-kstart+1,'int8');
    matrix_dim=size(sub_P);
    
    prism_p=(reshape(pp(i,:,:),3,6))';
    bbx=reshape(prism_facets_bbx(i,:,:),8,6);
    prism_n=reshape(prism_normals(i,:,:),3,8)';
    sub_P=P(istart:iend,jstart:jend,kstart:kend);
    for ii=istart:iend
        for kk=kstart:kend
            for jj=jstart:jend
                if P(ii,jj,kk)==0
%                     rp1 = mapijk2xyz(ii,jj,kk,[dx dy dz],[nrow ncol npln],[xmin ymin zmin]);
                    rp1 = [(jj-2)*dx+xmin (nrow-ii-1)*dy+ymin (npln-kk-1)*dz+zmin];
                    rp2 = [infX rp1(2) rp1(3)];
                    t1=tic;
                    [st intpnts all_ok intersection_status]=...
                        intersect_ray_shell_mex(rp1,rp2,prism_p,double(prism),tiny,bbx,prism_n,list,indexing,mydir);    
                        t2=t2+toc(t1);
%                     maxX = max(prism_p(:,1)); minX = min(prism_p(:,1));
%                     st = involume_mex(rp1, double(prism), prism_p, 200, minX, maxX, tiny);
                    if st==1
                        P(ii,jj,kk)=NA;
                    end
                    t3=tic;
                    sub_P=tag_row_subzone(rp1,intpnts,intersection_status,sub_P,st,...
                                      ii-istart+1,jj-jstart+1,kk-kstart+1,...
                                      [dx dy dz],matrix_dim,prisms_bbx(i,1:3),prisms_bbx(i,4:6));
                    t4=t4+toc(t3);
                    break
                end
            end
        end
    end
    P(istart:iend,jstart:jend,kstart:kend)=sub_P;
end
close(h);
cprintf([0 0 1],'\b%s\n\n',' done.')
% fprintf('  Time spent in intersect_ray_shell_mex: %4.6f\n  Time spent in tag_row_subzone: %4.6f\n',t2,t4); 
% fprintf('  Avg time spent per face for intersect_ray_shell_mex: %4.6f\n',t2/nf);
% fprintf('  Avg time spent per face for tag_row_subzone: %4.6f\n\n',t4/nf);

fprintf('\n%s\n','Exitting ExpandBoundaryBufferZone() function.')
fprintf('%s\n','===========================================')


















function P=tag_row_subzone(p0,intpnts,intersection_status,P,st,ii,jj,kk,resolution,matrix_dim,llc,urc)
global NA outside
dx=resolution(1); dy=resolution(2); dz=resolution(3);
nrow=matrix_dim(1); ncol=matrix_dim(2); npln=matrix_dim(3);
xmin=llc(1); ymin=llc(2); zmin= llc(3);

[intpnts sortidx]= sortrows(intpnts);
intersection_status=intersection_status(sortidx);
intpnts = [intpnts;p0+[abs(urc(1))*2 0 0]];
intersection_status(length(intersection_status)+1)=0;

toggle_worthy=[0 3];
if st == 0
    curr_state = 0;
else
    curr_state = NA;
end

starti = 1;
for J = jj:ncol
%     pp = mapijk2xyz(ii,J,kk,[dx dy dz],[nrow ncol npln],[xmin ymin zmin]);
    pp = [(J-0.5)*dx+xmin (nrow-ii+0.5)*dy+ymin (npln-kk+0.5)*dz+zmin];
    if pp(1) > intpnts(starti,1) % && ismember(intersection_status(cc),toggle_worthy); 
        % Pixel ii,J,k is within the prism, it should be tagged as a buffer
        % pixel
        crossed = 0;
        endi=starti;
        for i=starti:size(intpnts,1)
            if pp(1)<intpnts(i,1)
                endi = i;
                break
            end
        end
        for i=starti:(endi-1)
            % count the number of valid crossings
            if ismember(intersection_status(i),toggle_worthy)
                crossed = crossed + 1;
            end
        end
        if mod(crossed,2)==1 
            if curr_state == 0, curr_state=NA;
            else curr_state=0; end
        end
        if P(ii,J,kk)==0
            P(ii,J,kk)=curr_state;
        end
        starti = endi;
    else % if we haven't crossed over an intersection point yet we should keep assigning the same tag code
        if P(ii,J,kk)==0
            P(ii,J,kk)=curr_state;
        end
    end
end

