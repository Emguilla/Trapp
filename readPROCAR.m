clear
close all
clc

filename='PROCAR';
path='./';
fid=fopen([path,filename],'r');
a=fgetl(fid);
disp(a)
Data=textscan(fgetl(fid),'%s %s %s %d %s %s %s %d %s %s %s %d',1,'commentStyle','%');
nkpts=Data{4}(1);
nbands=Data{8}(1);
natoms=Data{12}(1);
for p=1:2
    a=fgetl(fid);
    Data=textscan(fgetl(fid),'%s %d %s %f %f %f %s %s %f',1,'commentStyle','%');
    k(p,:)=[Data{4}(1) Data{5}(1) Data{6}(1)];
    w(p)=Data{9}(1);
    a=fgetl(fid);
    for q=1:nbands
        disp(['band no ',num2str(q)])
        Data=textscan(fgetl(fid),'%s %d %s %s %f %s %s %f',1,'commentStyle','%');
        E(p,q)=Data{5}(1);
        Occ(p,q)=Data{8}(1);
        a=fgetl(fid);
        a=fgetl(fid);
        disp('Line 4')
        a=fgetl(fid);
        disp(a)
        disp('Line 5')
        a=fgetl(fid);
        disp(a)
        disp('Line 6')
        a=fgetl(fid);
        disp(a)
        a=fgetl(fid);
    end
end
fclose(fid);