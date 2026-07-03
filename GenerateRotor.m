function GenerateRotor(path,idxRot,idxCenter,alpha_deg)
%==================================================================================================================================%
% GenerateRotor.m:  Computation and writing of rotor rotated images (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (03/07/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   'path':         path to the vibration directory
%   'idxRot':       indices of the atoms involved in the rotation
%   'idxCenter':    index of the atom in the center of the rotor
%   'alpha_deg':    Array of rotation angles (in degrees)
%==================================================================================================================================%

% Check if the rotation has already been computed to avoid 
if exist('Rotation','dir')
    error('Rotation directory already exists!')
end
if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end
mkdir([path,'Rotation'])

% Get POSCAR
POSCAR=readPOSCAR([path,'POSCAR']);
% Get origin of the rotor
origin=POSCAR.positions(idxCenter,:);
% Compute rotation axis
n=origin-sum(POSCAR.positions(idxRot,:))/3;
% Extraction of the rotor as a separate POSCAR
CH3=subPOSCAR(POSCAR,idxRot);

% Rotation of the Rotor and replacement of the initial positions with the positions of the rotated rotor for each angle
for p=2:length(alpha_deg)
    nCH3=rotn(origin,n,alpha_deg(p),CH3);
    POSCAR(p)=POSCAR(1);
    POSCAR(p).positions(idxRot,:)=nCH3.positions;
end

% Rendering of the rotation
ReactionRendering(POSCAR,[0 0.5],[0.2 0.8],[0.15 0.55],'save_data','infinite','rotate_view_angle',[0 90]);

% Writing of the rotated angle to POSCAR files
for p=1:length(POSCAR)
    mkdir([path,'Rotation/',num2str(p,'%03i')]);
    writePOSCAR(POSCAR(p),[path,'Rotation/',num2str(p,'%03i'),'/POSCAR'])
end




