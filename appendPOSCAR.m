function POSCAR=appendPOSCAR(POSCAR1,POSCAR2)
%==================================================================================================================================%
% appendPOSCAR.m:   Append a POSCAR structure to another (v0.1.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (25/08/2025) - Creation
%       author: EYG
%   version 0.1.1 (29/11/2025) - Fix of the issue of the field constraint of POSCAR not being included in what was merged.
%       contrib: EYG
%==================================================================================================================================%
% args:
%   POSCAR1, POSCAR2:   POSCAR structure
%==================================================================================================================================%
load('constant_fund.mat','ptable');

% The resulting POSCAR will inherit all the generic info from the POSCAR1 structure (header, cell, constraints, coordinates, ...)
POSCAR.Title=POSCAR1.Title;
POSCAR.acell=POSCAR1.acell;
POSCAR.vec=POSCAR1.vec;
POSCAR.Selective_dynamics=POSCAR1.Selective_dynamics;
POSCAR.coord=POSCAR1.coord;

% The different chemical elements and their respective number of atoms is extracted from the POSCAR1 structure
count=0;
for p=1:length(POSCAR1.chemicals)
    chemicals{p}=POSCAR1.chemicals{p};
    n_chemicals(p)=POSCAR1.n_chemicals(p);
    positions{p}=POSCAR1.positions(count+1:count+POSCAR1.n_chemicals(p),:);
    if ~isnan(POSCAR1.constraint)
        constraint{p}=POSCAR1.constraint(count+1:count+POSCAR1.n_chemicals(p),:);
    end
    count=count+POSCAR1.n_chemicals(p);
end

% Here we try to match the type of chemical element from one POSCAR structure to the next. If there is a match, the n_chemicals
% field increased by the amount of atom in the POSCAR2 structure.
count=0;
for q=1:length(POSCAR2.chemicals)
    isdone=false;
    for p=1:length(POSCAR1.chemicals)
        if strcmpi(POSCAR1.chemicals{p},POSCAR2.chemicals{q})
            n_chemicals(p)=n_chemicals(p)+POSCAR2.n_chemicals(q);
            positions{p}(end+1:end+POSCAR2.n_chemicals(q),:)=POSCAR2.positions(count+1:count+POSCAR2.n_chemicals(q),:);
            if ~isnan(POSCAR1.constraint)
                constraint{p}(end+1:end+POSCAR2.n_chemicals(q),:)=POSCAR2.constraint(count+1:count+POSCAR2.n_chemicals(q),:);
            end
            count=count+POSCAR2.n_chemicals(q);
            isdone=true;
        end
    end
    % if there are additional chemicals in the POSCAR2 structure, they are simply appended to the list of chemicals
    if ~isdone
        chemicals{end+1}=POSCAR2.chemicals{q};
        positions{end+1}=POSCAR2.positions(count+1:count+POSCAR2.n_chemicals(q),:);
        if ~isnan(POSCAR1.constraint)
            constraint{end+1}=POSCAR2.constraint(count+1:count+POSCAR2.n_chemicals(q),:);
        end
        count=count+POSCAR2.n_chemicals(q);
    end
end
POSCAR.chemicals=chemicals;
POSCAR.n_chemicals=n_chemicals;

% The positions of the second file are appended to the ones of the previous file. Note that the new positions occurs "in batches",
% to match the sequence of chemical elements.
count=0;
for p=1:length(positions)
    pos_tmp=positions{p}(:,:);
    POSCAR.positions(count+1:count+POSCAR.n_chemicals(p),:)=positions{p}(:,:);
    if POSCAR.Selective_dynamics
        POSCAR.constraint(count+1:count+POSCAR.n_chemicals(p),:)=constraint{p}(:,:);
    else
        POSCAR.constraint=NaN;
    end
    count=count+POSCAR.n_chemicals(p);
end

% By default, the positions considered are expressed in cartesian coordinates. Once the full list of position is obtained, the
% direct coordinates are computed using the lattice vector of POSCAR1.
for p=1:sum(POSCAR.n_chemicals)
    POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
end

% The symbol field is completed using the sequence of the chemical elements
k=1;
for p=1:length(POSCAR.chemicals)
    for q=1:POSCAR.n_chemicals(p)
        POSCAR.symbols{k}=POSCAR.chemicals{p};
        k=k+1;
    end
end

% Assignement of a mass to each possible chemical element (data from PTable.com as of 20/08/2025)
for p=1:sum(POSCAR.n_chemicals)
    for q=1:118
        if strcmpi(POSCAR.symbols{p},ptable.symbol{q})
            POSCAR.mass(p)=ptable.mass(q);
            POSCAR.Z(p)=q;
        end
    end
end
