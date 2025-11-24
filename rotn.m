function POSCAR=rotn(origin,n,alpha_deg,POSCAR)
%==================================================================================================================================%
% rotn.m: rotation of a structure around an arbitrary axis n (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (21/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   origin:     position of the atom which anchor the rotor to the surface
%   n:          direction of the rotation axis of the rotor
%   alpha_deg:  rotation angle (in degree)
%   POSCAR:     POSCAR structure
%==================================================================================================================================%
[POSCAR,~,theta,phi]=RotorProjectionZ(origin,n,POSCAR);
for p=1:sum(POSCAR.n_chemicals)
    POSCAR.positions(p,:)=rotz(alpha_deg)*POSCAR.positions(p,:)';
    Ry=roty(theta*180/pi);
    POSCAR.positions(p,:)=Ry*POSCAR.positions(p,:)';
    Rz=rotz(phi*180/pi);
    POSCAR.positions(p,:)=Rz*POSCAR.positions(p,:)';
    POSCAR.positions(p,:)=POSCAR.positions(p,:)+origin;
end
end