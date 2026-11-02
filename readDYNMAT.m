function DYNMAT=readDYNMAT(filename)
%==================================================================================================================================%
% readPOSCAR.m: Read of a DYNMAT file from VASP and creation of a DYNMAT structure file in MatLab (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (12/10/2026) - Creation
%       author: EYG
%   version 0.1.1 (12/10/2026) - Add a way to retrieve the path to the folder based on the filename
%       contrib: EYG
%   version 0.1.2 (21/10/2026) - Add a default case for the path variable
%       contrib: EYG
%   version 0.2 (02/11/2026) - Add a "mask" field to the DYNMAT structure to easily find which DOF were processed. In addition, all
%       author: EYG             variables are put into the DYNMAT structure at the very end of the script.
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
    POSCAR=readPOSCAR([path,'POSCAR']);
else
    POSCAR=NaN;
end

% file opening
fid=fopen(filename,'r');

% Reading of the header, containing the number of atom types, the total amount of atoms, the number of DOF, and the masses of each
% atom species
Data=textscan(fid,'%d',3,'commentStyle','%');
Natoms=Data{1}(2);
Ntype=Data{1}(1);
NDOF=Data{1}(3);
Data=textscan(fid,'%f',Ntype,'commentStyle','%');
masses=Data{1};
Data=textscan(fid,'%d',2,'commentStyle','%');

% Initialisation of the mask
mask=false(DYNMAT.Natoms,3);

% start reading finite-difference results
p=0;

while ~isempty(Data{1})
    p=p+1;
    idx(p,:)=Data{1};
    mask(idx(p,1),idx(p,2))
    Data=textscan(fid,'%f',3,'commentStyle','%');
    dh{idx(p,1),idx(p,2)}=Data{1};
    Data=textscan(fid,'%f',3*DYNMAT.Natoms,'commentStyle','%');
    dF{idx(p,1),idx(p,2)}=reshape(Data{1},3,DYNMAT.Natoms)';
    Data=textscan(fid,'%d',2,'commentStyle','%');
end
fclose(fid);

for p=1:DYNMAT.Natoms
    for q=1:3
        if ~mask(p,q)
            dh{p,q}=NaN;
            dF{p,q}=NaN;
        end
    end
end

DYNMAT.POSCAR=POSCAR;
DYNMAT.Natoms=Natoms;
DYNMAT.Ntype=Ntype;
DYNMAT.NDOF=NDOF;
DYNMAT.masses=masses;
DYNMAT.mask=mask;
DYNMAT.dh=dh;
DYNMAT.dF=dF;
