function F=readMaxForces(path,varargin)
%==================================================================================================================================%
% readEnergy.m: Extraction of the maximum force exerced on the atoms in the system from the raw output of VASP (v0.2)
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
% check if the OUTCAR file has already been processed. If not, use the home-made grep function for MatLab
if ~exist([path,'MaxForces.dat'],'file')
    grep([path,'OUTCAR'],'RMS','fwd',0,'prt',[path,'MaxForces.dat']);
end
fid=fopen([path,'MaxForces.dat']);
Data=fgetl(fid);
k=1;
while ischar(Data)
    % Since each line looks like:
    % "  FORCES: max atom, RMS     0.003924    0.000473"
    % Only the character from 29 to 36 are converted into a double
    F(k)=str2double(Data(29:36));
    Data=fgetl(fid);
    k=k+1;
end
fclose(fid);
if ~save
    delete([path,'MaxForces.dat'])
end
end