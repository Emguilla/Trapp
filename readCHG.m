function [POSCAR,CHG]=readCHG(filename)
%==================================================================================================================================%
% readPOSCAR.m: Read of a CHG file from VASP and creation of a CHG structure file in MatLab (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (26/04/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   filename:   path + name of the file to be read as a CHG file (works also for CHGCAR)
%==================================================================================================================================%
% Load periodic table to add masses to POSCAR structures
load('ptable.mat','ptable');
% file opening
fid=fopen(filename,'r');

% Reading of the header
POSCAR.Title=fgets(fid);
POSCAR.Title=strtrim(POSCAR.Title(1:length(POSCAR.Title)-1));

% Reading of lattice parameter
Data=textscan(fid,'%f',1,'commentStyle','%');
POSCAR.acell=Data{1};

% Reading of the vector delimiting the periodical cell
Data=textscan(fid,'%f',3,'commentStyle','%');
POSCAR.vec(1,1:3)=Data{1:3};
Data=textscan(fid,'%f',3,'commentStyle','%');
POSCAR.vec(2,1:3)=Data{1:3};
Data=textscan(fid,'%f',3,'commentStyle','%');
POSCAR.vec(3,1:3)=Data{1:3};

% Rearrange the size of vectors so that the lattice parameter is always 1
POSCAR.vec=POSCAR.vec*POSCAR.acell;
POSCAR.acell=1;

% Reading of the different chemical species present in the structure
GO=true;
k=0;
while GO==true
    Data=textscan(fid,'%s',1,'commentStyle','%');
    if isempty(str2num(Data{1}{1}))==true
        k=k+1;
        POSCAR.chemicals{k}=Data{1}{1};
    else
        GO=false;
    end
end

% Reading of the number of each element previously read
POSCAR.n_chemicals(1)=str2num(Data{1}{1});
if k>=2
    for p=2:k
        Data=textscan(fid,'%f',1,'commentStyle','%');
        POSCAR.n_chemicals(p)=Data{1};
    end
end
k_XDATCAR=1;
Data=fgetl(fid);
% Reading of the constraints or type of coordinates (direct or cartesian)
Data=fgetl(fid);
if contains(lower(Data),'selective')
    POSCAR.Selective_dynamics=true;
    Data=fgetl(fid);
else
    POSCAR.Selective_dynamics=false;
end
if contains(lower(Data),'direct')
    POSCAR.coord.Direct=true;
    POSCAR.coord.Cartesian=false;
elseif contains(lower(Data),'cartesian')
    POSCAR.coord.Direct=false;
    POSCAR.coord.Cartesian=true;
end

% Reading of the positions of the atoms and the associated constraints
if k_XDATCAR>1
    fgetl(fid);
    Data=fgetl(fid);
    if isempty(Data)||~strcmpi(Data(1),'D')
        Data=fgetl(fid);
        readmultiple=false;
    end
end
if POSCAR.Selective_dynamics==true
    Data=textscan(fid,'%f %f %f %s %s %s',sum(POSCAR.n_chemicals),'commentStyle','%');
    POSCAR.positions=[Data{1} Data{2} Data{3}];
    POSCAR.constraint=[str2bool(cell2mat(Data{4})) str2bool(cell2mat(Data{5})) str2bool(cell2mat(Data{6}))];
% Reading of the positions of the atoms (no constraints found)
else
    Data=textscan(fid,'%f %f %f',sum(POSCAR.n_chemicals),'commentStyle','%');
    POSCAR.positions=[Data{1} Data{2} Data{3}];
    POSCAR.constraint=NaN;
end
% Translation of direct to cartesian and vice-versa
if POSCAR.coord.Direct==true
    POSCAR.xred=POSCAR.positions;
    for p=1:sum(POSCAR.n_chemicals)
        POSCAR.positions(p,:)=POSCAR.acell*POSCAR.vec'*POSCAR.positions(p,:)';
    end
elseif POSCAR.coord.Cartesian==true
    for p=1:sum(POSCAR.n_chemicals)
        POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
    end
end

% Shift of stray atoms back into the cell (periodic boundary conditions)
POSCAR=POSCARfolding(POSCAR);
k=1;
        
% Based on the list of chemical element and their number, it is possible to deduce which positions correspond to which 
% element. However it is much more practical to explicitly associate to each position a chemical symbol
for p=1:length(POSCAR.chemicals)
    for q=1:POSCAR.n_chemicals(p)
        POSCAR.symbols{k}=POSCAR.chemicals{p};
        k=k+1;
    end
end

% Assignement of a mass and an atomic number to each possible chemical element (data from PTable.com as of 20/08/2025)
for p=1:sum(POSCAR.n_chemicals)
    for q=1:118
        if strcmpi(POSCAR.symbols{p},ptable.symbol{q})
            POSCAR.mass(p)=ptable.mass(q);
            POSCAR.Z(p)=q;
        end
    end
end

% Once the POSCAR is read, the CHG file is read
line=fgetl(fid);
line=fgetl(fid);

% Read volumetric data size
Data=textscan(fid,'%d',3,'commentStyle','%');
size_CHG=Data{1};

% Read volumetric data
Data=textscan(fid,'%f ',prod(size_CHG),'commentStyle','%');
vol=abs(dot(POSCAR.vec(1,:),cross(POSCAR.vec(2,:),POSCAR.vec(3,:))));
total_dens=reshape(Data{1},size_CHG(1),size_CHG(2),size_CHG(3))/vol;
CHG.n_tot=total_dens;
line=fgetl(fid);
Data=textscan(fid,'%d',3,'commentStyle','%');
size_CHG=Data{1};
% If the calculation is spin-polarised, the second set of volumetric data contains the magnetisation
if ~size(size_CHG,1)==0
    Data=textscan(fid,'%f ',prod(size_CHG),'commentStyle','%');
    mag=reshape(Data{1},size_CHG(1),size_CHG(2),size_CHG(3))/vol;
    % Computation of the spin up and spin down electronic densities
    CHG.n_u=(total_dens+mag)/2;
    CHG.n_d=(total_dens-mag)/2;
    CHG.n_mag=mag;
end

fclose(fid);