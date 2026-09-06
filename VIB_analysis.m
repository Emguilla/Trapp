function VIB=VIB_analysis(path)
%==================================================================================================================================%
% VIB_analysis.m: Data extraction of the (multiple) phonon calculations (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (03/07/2026) - Creation
%       author: EYG
%   version 0.2 (06/09/2026) - The output is now a structure array instead of a structure containing arrays.
%       author: EYG
%==================================================================================================================================%
% args:
%   path:   Location of the directory where the phonon calculation are located
%==================================================================================================================================%

% Check that "path" ends with a slash or backslash
if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end

% If the directory contains a DYNMAT file, it means it is a tight transition state or a groundstate calculation
if exist([path,'DYNMAT'],'file')
    ldir.name='.';
% Otherwise, it is a loose transition state and contains multiple directories, one per image.
else
    ldir=subdir_ordering(path);
end
n_TS=length(ldir);
for p=1:n_TS
    VIB(p).POSCAR=readPOSCAR([path,ldir(p).name,'/POSCAR']);
    if n_TS>1
        VIB(p).POSCAR.Title=ldir(p).name;
    end
    if exist([path,ldir(p).name,'/VibRotor.dat'],'file')
        Data=load([path,ldir(p).name,'/VibRotor.dat']);
        if Data(1)==1
            VIB(p).Rotor.mode='free';
        elseif Data(1)==2
            VIB(p).Rotor.mode='hindered';
        end
        VIB(p).Rotor.Vib_idx=Data(2);
        VIB(p).Rotor.sigma=Data(3);
        VIB(p).Rotor.I=Data(4);
        VIB(p).Rotor.RotEnergyBarrier=Data(5);
    else
        VIB(p).Rotor=NaN;
    end
    Energy=readEnergy([path,ldir(p).name,'/'],'save',true);
    VIB(p).E=Energy(1);
    VIB(p).nu=readFrequencies([path,ldir(p).name,'/']);
    if sum(VIB(p).nu==0)~=3
        error('Severe problem found: there are fewer or more than three zero-frequency vibrational mode!')
    end
    VIB(p).date=datetime('today');
end

% The Vibration structure is stored in a ".mat" file for convenience
save([path(1:end-1),'.mat'],'VIB')
