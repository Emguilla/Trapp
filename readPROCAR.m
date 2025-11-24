filename='PROCAR';
path='./';
fid=fopen([path,filename],'r');
fgetl(fid);
Data=textscan(fid,'%s %s %s %d',3,'commentStyle','%');
nkpts=Data{4}(1);
nbands=Data{4}(2);
natoms=Data{4}(3);
fgetl(fid)
Data=textscan(fid,'%s %d %s %f %f %f %s %s %f',1,'commentStyle','%');
k(p,:)=[Data{4}(1) Data{5}(1) Data{6}(1)];
w(p)=Data{9}(1);
for q=1:nbands
    Data=textscan(fid,'%s %d %s %s %f %s %s %f',1,'commentStyle','%');
    E(p,q)=Data{5}(1);
    Occ(p,q)=Data{8}(1);
    
end
fclose(fid);