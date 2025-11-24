function [POSCAR,r,theta,phi]=RotorProjectionZ(origin,n,POSCAR)
%==================================================================================================================================%
% RotorProjectionZ.m: translation of a rotor and rotation of its axis to match the z-axis (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (21/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   origin: position of the atom which anchor the rotor to the surface
%   n:      direction of the rotation axis of the rotor
%   POSCAR: POSCAR structure
%==================================================================================================================================%
% calculation of the spherical coordinate of the axis vector
r=norm(n);
theta=acos(n(3)/r);
phi=acos(n(1)/sqrt(n(1)^2+n(2)^2));
% Computation of the rotation matrix to move the rotor to the origin and the z-axis
Rz=rotz(-phi);
Ry=roty(-theta);
% The positions (in cartesian coordinates) are first shifted to the origin, then the structure is rotated around the z-axis, and 
% finally around the y-axis. The direct coordinate are modified accordingly.
for p=1:sum(POSCAR.n_chemicals)
    POSCAR.positions(p,:)=POSCAR.positions(p,:)-origin;
    POSCAR.positions(p,:)=Rz*POSCAR.positions(p,:)';
    POSCAR.positions(p,:)=Ry*POSCAR.positions(p,:)';
    POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
end
end