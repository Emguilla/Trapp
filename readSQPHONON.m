function XDATCAR=readSQPHONON(filename)
%==================================================================================================================================%
% readSQPHONON.m: Read of a SQ_PHONON file from HIVE and creation of a POSCAR array structure file in MatLab (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (06/09/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   filename:   path + name of the file to be read as a SQ_PHONON
%==================================================================================================================================%

% Find all files matching the pattern provided by "filename"
ldir=dir(filename);
if isempty(ldir)
    error(['Attempt at reading ',filename,' failed: no such a file'])
end

% Load periodic table to add masses to POSCAR structures
load('ptable.mat','ptable');
% Loop over all the files that matches the pattern of filename
for s=1:length(ldir)
    % file opening
    fid=fopen([ldir(s).folder,'/',ldir(s).name],'r');
    % Reading of the HIVE header
    fgetl(fid);
    fgetl(fid);
    Data=textscan(fid,'%s %s %s %s %d',1,'commentStyle','%');
    n_images=Data{5};
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    for r=1:n_images
        fgetl(fid);
        fgetl(fid);
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
        if POSCAR.Selective_dynamics==true
            Data=textscan(fid,'%f %f %f %s %s %s\n',sum(POSCAR.n_chemicals),'commentStyle','%');
            if isempty(Data{1})
                readmultiple=false;
            else
                POSCAR.positions=[Data{1} Data{2} Data{3}];
                POSCAR.constraint=[str2bool(cell2mat(Data{4})) str2bool(cell2mat(Data{5})) str2bool(cell2mat(Data{6}))];
            end
        % Reading of the positions of the atoms (no constraints found)
        else
            Data=textscan(fid,'%f %f %f\n',sum(POSCAR.n_chemicals),'commentStyle','%');
            if isempty(Data{1})
                readmultiple=false;
            else
                POSCAR.positions=[Data{1} Data{2} Data{3}];
            end
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
        XDATCAR(s,r)=POSCAR;
    end
    fclose(fid);
end
end