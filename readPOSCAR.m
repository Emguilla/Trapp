function POSCAR=readPOSCAR(filename)
%==================================================================================================================================%
% readPOSCAR.m: Read of a POSCAR file from VASP and creation of a POSCAR structure file in MatLab (v0.6)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.2 (20/80/2025) - Add masses (in uma) to the POSCAR structure
%       author: EYG
%   version 0.2.0.1 (21/08/2025) - Correction to a typo in a variable name
%       contrib: EYG
%   version 0.3 (25/08/2025) - Add the atomic number to the POSCAR structure under the field Z
%       author: EYG
%   version 0.4 (29/08/2025) - Add the possibility to read XDATCAR files
%       author: EYG
%   version 0.4.1 (31/08/2025) - For non-selective dynamics POSCAR structures, the field "constraint" is set to NaN.
%       contrib: EYG
%   version 0.5 (10/09/2025) - Modification of the reading loop for XDATCAR to automatically accomodate for POSCAR or XDATCAR. The
%       author: EYG             optional argument 'XDATCAR' has thus been removed.
%   version 0.5.1 (05/12/2025) - Removal of the vestigial bit that handle varargin and removal of varargin as an input
%       contrib: EYG
%   version 0.5.2 (01/01/2026) - The folding of stray atoms now occurs in a separate function "POSCARfolding.m"
%       contrib: EYG
%   version 0.6 (27/01/2026) - All the files that matches the pattern of "filename" are now read and stored in a POSCAR array. This
%       author: EYG             allows for reading of multiple POSCAR at once.
%==================================================================================================================================%
% args:
%   filename:   path + name of the file to be read as a POSCAR (works for CONTCAR and XDATCAR as well)
%==================================================================================================================================%

% Find all files matching the pattern provided by "filename"
ldir=dir(filename);

% Load periodic table to add masses to POSCAR structures
load('constant_fund.mat','ptable');

% Loop over all the files that matches the pattern of filename
for r=1:length(ldir)
    % Initialisation of the default parameters
    readmultiple=true;
    % file opening
    disp([ldir(r).folder,'/',ldir(r).name])
    fid=fopen([ldir(r).folder,'/',ldir(r).name],'r');
    
    % Reading of the header
    POSCAR.Title=fgets(fid);
    disp(POSCAR.Title)
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
    
    % This only happens if the reading of the positions was successfully done. Otherwise, it means the reading is finished.
    while readmultiple
        % Reading of the positions of the atoms and the associated constraints
        if k_XDATCAR>1
            fgetl(fid);
            Data=fgetl(fid);
            if isempty(Data)||~strcmpi(Data(1),'D')
                Data=fgetl(fid);
                readmultiple=false;
            end
        end
        if readmultiple
            if POSCAR.Selective_dynamics==true
                Data=textscan(fid,'%f %f %f %s %s %s',sum(POSCAR.n_chemicals),'commentStyle','%');
                if isempty(Data{1})
                    readmultiple=false;
                else
                    POSCAR.positions=[Data{1} Data{2} Data{3}];
                    POSCAR.constraint=[str2bool(cell2mat(Data{4})) str2bool(cell2mat(Data{5})) str2bool(cell2mat(Data{6}))];
                end
            % Reading of the positions of the atoms (no constraints found)
            else
                Data=textscan(fid,'%f %f %f',sum(POSCAR.n_chemicals),'commentStyle','%');
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
            
            if k_XDATCAR==1
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
            end
            XDATCAR(k_XDATCAR)=POSCAR;
            k_XDATCAR=k_XDATCAR+1;
        end
    end
    fclose(fid);
    CONTCAR(r)=XDATCAR;
    disp(CONTCAR(r).Title)
end
POSCAR=CONTCAR;