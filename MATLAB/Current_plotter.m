dbdown=20
map='default';
cut_pec=-50000%-0.1011%-0.42;
eks_pec=2;
below_pec=0;

ith=5;
fid1 = fopen (['EigVec_',num2str(ith),'_freq_10600000000.000000.out'],'r');
fid2=fopen('pillbox_100000_node.inp','r');
fid3=fopen('pillbox_100000_elem.inp','r');

nnod2=fscanf(fid2,'%d',1);
npat=fscanf(fid3,'%d',1);

t1=cputime;
xyz2 = zeros(nnod2,3);
for i=1:nnod2
    tmp=fscanf(fid2,'%d',1);
    temp(1:3)=(fscanf(fid2,'%f',3))';
    xyz2(i,1)=temp(1);
    xyz2(i,2)=temp(2);
    xyz2(i,3)=temp(3);
end

t2=cputime;
display('surface node read time:')
t2-t1
t1=cputime;
ii=1;
eks=eks_pec;
cut=cut_pec;
below=below_pec;
ct2=zeros(npat,3);
if (below==1)
    for i=1:npat
        tmp=(fscanf(fid3,'%d',1))';
        ct2(ii,1:3)=(fscanf(fid3,'%d',3))';
%         if (xyz2(ct2(ii,1),eks)<=cut & xyz2(ct2(ii,2),eks)<=cut & xyz2(ct2(ii,3),eks)<=cut)
            ii=ii+1;
            %if any one of the nodes are outside the cut dont do anything
%         end
    end
else
    for i=1:npat
        tmp=(fscanf(fid3,'%d',1))';
        ct2(ii,1:3)=(fscanf(fid3,'%d',3))';
        %if (xyz2(ct2(ii,1),eks)>=cut & xyz2(ct2(ii,2),eks)>=cut & xyz2(ct2(ii,3),eks)>=cut ...
        %        & xyz2(ct2(ii,1),eks)<=5.999 & xyz2(ct2(ii,2),eks)<=5.999 & xyz2(ct2(ii,3),eks)<=5.999)
%         if (xyz2(ct2(ii,1),eks)>=cut & xyz2(ct2(ii,2),eks)>=cut & xyz2(ct2(ii,3),eks)>=cut)
        
            ii=ii+1;
            %if any one of the nodes are outside the cut dont do anything
%         end
    end
end
count_pec=ii-1
t2=cputime;
display('surface element read and processing time:')
t2-t1


% PROCESSING CURRENT INFORMATION.......................................

X=xyz2(:,1);
Y=xyz2(:,2);
Z=xyz2(:,3);

CV=zeros(nnod2,1);
% color for each node
CV(1:nnod2)=fscanf(fid1,'%f',nnod2)';

figure(3);
colormap(map);
CV = CV./max((abs(CV)));
CV =20*log10(abs(CV));

clf
a = trisurf(ct2(1:count_pec,1:3),X,Y,Z,CV); %pause;

c = colorbar('location','EastOutside');
caxis([-dbdown 0]);
axis equal;
shading interp;
axis off;
grid off;
view(-55,40);

fname=sprintf('%s%d.%s','curr_mode_',num2str(ith),'jpg');
print ('-djpeg','-r300',fname)


fclose('all');

axis equal
shading interp
