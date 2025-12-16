function F=readForces(path,varargin)
%==================================================================================================================================%
% readForces.m: Extraction of the forces acting on a system from the raw output of VASP (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (28/08/2025) - Creation
%       author: EYG
%   version 0.2 (11/09/2025) - The recording of the grep search is now optional, and default is no save
%       author: EYG
%==================================================================================================================================%
% args:
%   path:   Location of the directory where the OUTCAR file is stored
%==================================================================================================================================%
save=false;
% read optional arguments
if exist('varargin')
    for p=1:2:length(varargin)
        switch varargin{p}
            case 'save'
                save=varargin{p+1};
        end
    end
end
POSCAR=readPOSCAR([path,'/POSCAR']);
n_atoms=sum(POSCAR.n_chemicals);
% check if the OUTCAR file has already been processed. If not, use the home-made grep function for MatLab
if ~exist([path,'Forces.dat'],'file')
    grep([path,'/OUTCAR'],'TOTAL-FORCE','fwd',n_atoms+3,'prt',[path,'/Forces.dat']);
end
fid=fopen([path,'/Forces.dat']);
GO=true;
k=1;
while GO==true
    Data=textscan(fid,'%f %f %f %f %f %f',(n_atoms+5));
    szd=size(Data{1});
    if szd(1)==0
        disp(Data{1})
        if fgetl(fid)==-1
            GO=false;
        end
    else
        F{k}=[Data{4} Data{5} Data{6}].*POSCAR.constraint;
        k=k+1;
    end
end
fclose(fid);
if ~save
    delete([path,'Forces.dat'])
end
end