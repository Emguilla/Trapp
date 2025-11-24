%function readDOSCAR

filename='DOSCAR';
path='./';
fid=fopen([path,filename],'r');
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
Data=textscan(fid,'%f',5,'commentStyle','%');
Emin=Data{1}(2);
Emax=Data{1}(1);
nE=Data{1}(3);
for p=1:nE
    Data=textscan(fid,'%f',3,'commentStyle','%');
    E(p)=Data{1}(1);
    DOS(p)=Data{1}(2);
end
fclose(fid)
plot(E,DOS)