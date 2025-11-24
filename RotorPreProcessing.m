function RotorPreProcessing(mode,path,idx_origin,idx_n,idx_rotor,sigma)
%==================================================================================================================================%
% RotorPreProcessing.m: Creation of the rotor file to be read by "hr_freq.m" or "fr_freq.m" (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (21/08/2025) - Creation
%       author: EYG
%   version 0.2 (25/08/2025) - Add the possibility to specify a vector instead of the index of the second atom on the rotation axis
%       author: EYG
%==================================================================================================================================%
% args:
%   mode:       'hindered' -> create a "hindered_rotor.dat" file.
%               'free' -> create a "free_rotor.dat" file.
%   path:       path to the POSCAR file and where the file are to be written
%   idx_origin: index (in the POSCAR structure) of the atom used as the anchor of the rotor adsorbate on the surface
%   idx_n:      integer scalar - indices (in the POSCAR structure) of the second atom that defines, along with the anchor atom at 
%                   the origin, the rotation axis of the adsorbate. 
%               double array - cartesian coordinate of the vector to be used as the rotational axis
%   idx_rotor:  indices (in the POSCAR structure) of all the atoms that are part of the rotor, including 
%   sigma:      number of equivalent positions around a 360 rotation of the rotor (e.g., sigma of CH3 is equal to 3)
%==================================================================================================================================%
% Default parameters
%   mode='hindered';
%   path='E:\Reactions\034_100_0A+CH3_to_1A+CH4\TS\';
%   filename=[path,'POSCAR'];
%   idx_origin=6;
%   idx_n=[6 36];
%   idx_rotor=[6 33:36];
%   sigma=3;
%==================================================================================================================================%
POSCAR=readPOSCAR(filename);
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
if strcmpi(mode,'hindered')
    alpha_deg=180/sigma;
    Rotor_60=rotn(origin,n,alpha_deg,Rotor_0);
    for p=1:length(idx)
        POSCAR.positions(idx(p),:)=Rotor_60.positions(p,:);
        POSCAR.xred(idx(p),:)=Rotor_60.xred(p,:);
    end
    writePOSCAR(POSCAR,[path,'POSCAR_rotn_',num2str(alpha_deg,'%03i')])
    fid=fopen([path,'hindered_rotor.dat'],'w');
elseif strcmpi(mode,'free')
    fid=fopen([path,'free_rotor.dat'],'w');
end
% For either mode, the value of 
fprintf(fid,'%03i',sigma);
fprintf(fid,'    %10.6e',Iz);
fclose(fid);
end
