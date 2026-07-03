function VIB=VIB_analysis(path)
%==================================================================================================================================%
% VIB_analysis.m: Data extraction of the (multiple) phonon calculations (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (03/07/2026) - Creation
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
    POSCAR(p)=readPOSCAR([path,ldir(p).name,'/POSCAR']);
    if exist([path,ldir(p).name,'/VibRotor.dat'],'file')
        Data=load([path,ldir(p).name,'/VibRotor.dat']);
        if Data(1)==1
            Rotor(p).mode='free';
        elseif Data(1)==2
            Rotor(p).mode='hindered';
        end
        Rotor(p).Vib_idx=Data(2);
        Rotor(p).sigma=Data(3);
        Rotor(p).I=Data(4);
        Rotor(p).RotEnergyBarrier=Data(5);
    else
        Rotor(p)=NaN;
    end
    E=readEnergy([path,ldir(p).name,'/'],'save',true);
    Energy(p)=E(1);
    nu{p}=readFrequencies([path,ldir(p).name,'/']);
    if sum(nu{p}==0)~=3
        error('Severe problem found: there are fewer or more than three zero-frequency vibrational mode!')
    end
end

% Construction of the VIB structure
VIB.POSCAR=POSCAR;
VIB.E=Energy;
VIB.Freq=nu;
VIB.Rotor=Rotor;

% The data is stored in a ".mat" file for convenience
save([path(1:end-1),'_VIB.mat'],'VIB')
