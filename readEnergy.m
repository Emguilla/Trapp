function E=readEnergy(path,varargin)
%==================================================================================================================================%
% readEnergy.m: Extraction of the energy of a system from the raw output of VASP (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.1.0.1 (20/08/2025) - Add a few comments
%       contrib: EYG
%   version 0.2 (28/08/2025) - Previously only the first energy of the system was read. Now all of them are read.
%       author: EYG
%   version 0.3 (11/09/2025) - The recording of the grep search is now optional, and default is no save
%       author: EYG
%==================================================================================================================================%
% args:
%   path:       Location of the directory where the OUTCAR file is stored
%   opt arg:    'save', followed by true or false to specify whether or not the grep search is to be saved.
%                   (default: false)
%==================================================================================================================================%
if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end
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
if ~exist([path,'Energy.dat'],'file')
    grep([path,'OUTCAR'],'free  energy','fwd',0,'prt',[path,'Energy.dat']);
end
fid=fopen([path,'Energy.dat']);
Data=fgetl(fid);
k=1;
while ischar(Data)
    % Since each line looks like:
    % "  free  energy   TOTEN  =     -1760.36803168 eV"
    % Only the character from end-20 to end-3 are converted into a double
    E(k)=str2double(Data(end-20:end-3));
    Data=fgetl(fid);
    k=k+1;
end
fclose(fid);
if ~save
    delete([path,'Energy.dat'])
end