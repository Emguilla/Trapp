function POSCAR=rotn(origin,n,alpha_deg,POSCAR)
%==================================================================================================================================%
% rotn.m: rotation of a structure around an arbitrary axis n (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (21/08/2025) - Creation
%       author: EYG
%   version 0.2 (21/12/2025) - Inclusion of the "RotorProjectionZ.m" code
%       author: EYG
%==================================================================================================================================%
% args:
%   origin:     position of the atom which anchor the rotor to the surface
%   n:          direction of the rotation axis of the rotor
%   alpha_deg:  rotation angle (in degree)
%   POSCAR:     POSCAR structure
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
% The arbitrary axis now matches the z-axis, accordingly, a rotation of alpha_deg is performed around the z-axis
Rn=rotz(-alpha_deg*pi/180);
for p=1:sum(POSCAR.n_chemicals)
    POSCAR.positions(p,:)=Rn*POSCAR.positions(p,:)';
    POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
end
% The structure is rotated around the z- and y-axis to orient the n-axis in its initial direction, and then positions (in cartesian 
% coordinates) are are shifted back to the anchor)
Ry=roty(theta);
Rz=rotz(phi);
for p=1:sum(POSCAR.n_chemicals)
    POSCAR.positions(p,:)=Ry*POSCAR.positions(p,:)';
    POSCAR.positions(p,:)=Rz*POSCAR.positions(p,:)';
    POSCAR.positions(p,:)=POSCAR.positions(p,:)+origin;
    POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
end
end