function writePOSCAR(POSCAR,filename,varargin)
%==================================================================================================================================%
% readPOSCAR.m: write a POSCAR MatLab structure in a file following VASP format (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (02/12/2025) - Creation
%       author: EYG
%   version 0.2 (05/12/2025) - Add optional argument to decide to sort or not the atoms in the POSCAR.
%==================================================================================================================================%
% args:
%   POSCAR:     POSCAR structure to be written
%   filename:   path + name of the file in which the POSCAR structure must be stored
%   opt. args:  'sort', followed by true or false to decide whether to sort the atoms.
%                   (default: atoms are grouped and sorted according to their atomic number)
%==================================================================================================================================%
% Initialisation of the default parameters
sort_opt=true;
% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'sort'
                sort_opt=varargin{p+1};
        end
    end
end
fid=fopen(filename,'w');
% Writing of the header
fprintf(fid,[POSCAR.Title,'\n']);
% Writing of lattice parameter
fprintf(fid,'%10.6f\n',POSCAR.acell);
% Writing of the vector delimiting the periodical cell
fprintf(fid,'    %10.6f    %10.6f    %10.6f\n',POSCAR.vec(1,:));
fprintf(fid,'    %10.6f    %10.6f    %10.6f\n',POSCAR.vec(2,:));
fprintf(fid,'    %10.6f    %10.6f    %10.6f\n',POSCAR.vec(3,:));
% sort atoms if not disabled
if sort_opt
    [lZ,ntype,ltype]=unique(POSCAR.Z);
    for p=1:length(lZ)
        fseq(p)=sum(ltype==ntype(p));
        species_seq{p}=POSCAR.symbols{ntype(p)};
    end
    [~,idx]=sort(ltype);
    POSCAR.positions=POSCAR.positions(idx,:);
    POSCAR.xred=POSCAR.xred(idx,:);
    POSCAR.symbols=POSCAR.symbols(idx);
    POSCAR.mass=POSCAR.mass(idx);
    POSCAR.Z=POSCAR.Z(idx);
    if ~isnan(POSCAR.constraint)
        POSCAR.constraint=POSCAR.constraint(idx,:);
    end
    POSCAR.n_chemicals=fseq;
    POSCAR.chemicals=species_seq;
end
% Writing of the different chemical species present in the structure and their respective amount
blank_space='     ';
for p=1:length(POSCAR.chemicals)
	fprintf(fid,[blank_space(1:end-length(POSCAR.chemicals{p})),POSCAR.chemicals{p}]);
end
for p=1:length(POSCAR.chemicals)
    fprintf(fid,'\n  %3d',POSCAR.n_chemicals(1));
end
% Writing of the constraints or type of coordinates (direct or cartesian) and the positions of the atoms
if POSCAR.Selective_dynamics==true
    fprintf(fid,'\nSelective Dynamics\n');
    fprintf(fid,'Cartesian');
    for p=1:sum(POSCAR.n_chemicals)
        fprintf(fid,'\n %10.6f %10.6f %10.6f',POSCAR.positions(p,:));
        for q=1:3
            if POSCAR.constraint(p,q)==true
                fprintf(fid,' T');
            else
                fprintf(fid,' F');
            end
        end
    end
else
    fprintf(fid,'\nCartesian\n');
    for p=1:sum(POSCAR.n_chemicals)
        fprintf(fid,' %10.6f %10.6f %10.6f\n',POSCAR.positions(p,:));
    end
end
fclose(fid);