function RotorPreProcessing(path,POSCAR,vibidx,idx_origin,idx_n,idx_rotor,sigma)
%==================================================================================================================================%
% RotorPreProcessing.m: Creation of the rotor file to be read by "hr_freq.m" or "fr_freq.m" (v0.3)
%==================================================================================================================================%
% Version history:
%   version 0.1 (21/08/2025) - Creation
%       author: EYG
%   version 0.2 (25/08/2025) - Add the possibility to specify a vector instead of the index of the second atom on the rotation axis
%       author: EYG
%   version 0.3 (19/02/2026) - Removal of the mode argument. This argument is now an input from the user based on its analysis of
%       author: EYG             the relative height of the rotational energy barrier. The vibrational mode to be replaced by the 
%                               VibRot mode is provided via the new input vibidx. For sanity purposes, HIVE_analysis is called to
%                               generate a rendering of that particular mode.
%==================================================================================================================================%
% args:
%   mode:       'hindered' -> create a "hindered_rotor.dat" file.
%               'free' -> create a "free_rotor.dat" file.
%   path:       path to the POSCAR file and where the file are to be written
%   vibidx:     Index of the vibrational mode of interest, starting from the lowest frequencies.
%   idx_origin: index (in the POSCAR structure) of the atom used as the anchor of the rotor adsorbate on the surface
%   idx_n:      integer scalar - indices (in the POSCAR structure) of the second atom that defines, along with the anchor atom at 
%                   the origin, the rotation axis of the adsorbate. 
%               double array - cartesian coordinate of the vector to be used as the rotational axis
%   idx_rotor:  indices (in the POSCAR structure) of all the atoms that are part of the rotor, including anchor and origin.
%   sigma:      number of equivalent positions around a 360 rotation of the rotor (e.g., sigma of CH3 is equal to 3)
%==================================================================================================================================%
% Default parameters
% path='C:\Users\emgui\OneDrive - Université de Namur\Systems\Diamond\100\Gam-only\02-passivated\LibReac\A0A0+CH3_A1A0+CH4\tTS';
% POSCAR='POSCAR';
% vibidx=2;
% idx_origin=6;
% idx_n=36;
% idx_rotor=[6 33:36];
% sigma=3;
%==================================================================================================================================%
if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end
% If the POSCAR input is a character string, the corresponding file is read into a POSCAR structure
if ischar(POSCAR)
    POSCAR=readPOSCAR(POSCAR);
end
idx=unique(idx_rotor);
Rotor_0=subPOSCAR(POSCAR,idx_rotor);
aLim=[min(Rotor_0.xred(:,1)) max(Rotor_0.xred(:,1))];
bLim=[min(Rotor_0.xred(:,2)) max(Rotor_0.xred(:,2))];
cLim=[min(Rotor_0.xred(:,3)) max(Rotor_0.xred(:,3))];
POSCARs=HIVE_analysis(vibidx);
for p=1:length(POSCARs)
    CONTCARs(p)=subPOSCAR(POSCARs(p),idx_rotor);
    aLim=[min([aLim(1) min(CONTCARs(p).xred(:,1))]) max([aLim(2) max(CONTCARs(p).xred(:,1))])];
    bLim=[min([bLim(1) min(CONTCARs(p).xred(:,2))]) max([bLim(2) max(CONTCARs(p).xred(:,2))])];
    cLim=[min([cLim(1) min(CONTCARs(p).xred(:,3))]) max([cLim(2) max(CONTCARs(p).xred(:,3))])];
end
ReactionRendering(CONTCARs,aLim,bLim,cLim);
% extract the origin and rotation axis vector from POSCAR structure. if the "idx_n" variable is a vector, this vector is taken as
% the rotation axis
origin=POSCAR.positions(idx_origin,:);
if isscalar(idx_n)
    n=POSCAR.positions(idx_n,:)-POSCAR.positions(idx_origin,:);
elseif length(idx_n)==3
    n=idx_n;
end
% Creation of a sub-POSCAR structure containing only the rotor
idx=unique(idx_rotor);
Rotor_0=subPOSCAR(POSCAR,idx);
% Rotation of the rotor to match the rotational axis and the z-axis, and translation of the anchor to the origin
[Rotor_z,~,~,~]=RotorProjectionZ(origin,n,Rotor_0);
% Computation of the inertia moment
[~,~,Iz]=InertiaMoment(Rotor_z);
% If the mode is "hindered", a rotation of 180/sigma degree is performed around the n-axis and the resulting structure is saved in the
% file "POSCAR_rotn_XXX"
alpha_deg=180/sigma;
Rotor_60=rotn(origin,n,alpha_deg,Rotor_0);
if exist([path,'Rotation/0000/OUTCAR'],'file')&&exist([path,'Rotation/',num2str(alpha_deg,'%04i'),'/OUTCAR'],'file')
    NEB=NEB_analysis('path',[path,'Rotation/']);
    W=max(NEB.energies)-min(NEB.energies);
    fprintf(['The rotational energy barrier seems to be about ',num2str(W*1000),' meV. Do you consider that to be a free or hindered rotation?\n 1) Free\n2) Hindered\n'])
    OK=false;
    while ~OK
        mode=input('==> ');
        switch mode
            case {1,2}
                OK=true;
            otherwise
                fprintf('Please enter 1 or 2.')
        end
    end
    % For either mode, the value of the index of the vibrational mode, the symmetry number and inertia moment are written.
    fid=fopen([path,'VibRotor.dat'],'w');
    fprintf(fid,'%i',mode);
    fprintf(fid,'   %04i',vibidx);
    fprintf(fid,'    %03i',sigma);
    fprintf(fid,'    %10.6e',Iz);
    fprintf(fid,'    %10.6e',W);
    fclose(fid);
else
    if ~exist([path,'Rotation'],'dir')
        mkdir([path,'Rotation'])
    end
    if ~exist([path,'Rotation/0000'],'dir')
        mkdir([path,'Rotation/0000'])
    end
    if ~exist([path,'Rotation/',num2str(alpha_deg,'%04i')],'dir')
        mkdir([path,'Rotation/',num2str(alpha_deg,'%04i')])
    end
    writePOSCAR(POSCAR,[path,'Rotation/0000/POSCAR'])
    for p=1:length(idx)
        POSCAR.positions(idx(p),:)=Rotor_60.positions(p,:);
        POSCAR.xred(idx(p),:)=Rotor_60.xred(p,:);
    end
    writePOSCAR(POSCAR,[path,'Rotation/',num2str(alpha_deg,'%04i'),'/POSCAR'])
end
end
