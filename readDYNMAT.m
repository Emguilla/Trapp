function DYNMAT=readDYNMAT(filename)
%==================================================================================================================================%
% readPOSCAR.m: Read of a DYNMAT file from VASP and creation of a DYNMAT structure file in MatLab (v0.1.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (12/10/2026) - Creation
%       author: EYG
%   version 0.1.1 (12/10/2026) - Add a way to retrieve the path to the folder based on the filename
%       contrib: EYG
%   version 0.1.2 (21/10/2026) - Add a default case for the path variable
%       contrib: EYG
%==================================================================================================================================%
% args:
%   filename:   path + name of the file to be read as a DYNMAT
%==================================================================================================================================%
path='./';
for p=1:length(filename)
    if strcmpi(filename(p),'/')||strcmpi(filename(p),'\')
        path=filename(1:p);
    end
end
if exist([path,'POSCAR'],'file')
    DYNMAT.POSCAR=readPOSCAR([path,'POSCAR']);
else
    DYNMAT.POSCAR=NaN;
end

% by default, we assume the DYNMAT to cover all DOF
partial=false;

% file opening
fid=fopen(filename,'r');

% Reading of the header, containing the number of atom types, the total amount of atoms, the number of DOF, and the masses of each
% atom species
Data=textscan(fid,'%d',3,'commentStyle','%');
DYNMAT.Natoms=Data{1}(2);
DYNMAT.Ntype=Data{1}(1);
DYNMAT.NDOF=Data{1}(3);
Data=textscan(fid,'%f',DYNMAT.Ntype,'commentStyle','%');
DYNMAT.masses=Data{1};
Data=textscan(fid,'%d',2,'commentStyle','%');

% start reading finite-difference results
p=0;

while ~isempty(Data{1})
    p=p+1;
    idx(p,:)=Data{1};
    Data=textscan(fid,'%f',3,'commentStyle','%');
    DYNMAT.dh{idx(p,1),idx(p,2)}=Data{1};
    Data=textscan(fid,'%f',3*DYNMAT.Natoms,'commentStyle','%');
    DYNMAT.dF{idx(p,1),idx(p,2)}=reshape(Data{1},3,DYNMAT.Natoms)';
    Data=textscan(fid,'%d',2,'commentStyle','%');
end
fclose(fid);