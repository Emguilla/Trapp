function [Ix,Iy,Iz]=InertiaMoment(POSCAR)
%==================================================================================================================================%
% InertiaMoment.m:  Calculation of the moment of inertia of a POSCAR structure (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (25/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   POSCAR: POSCAR structure
%==================================================================================================================================%
load('constant_fund.mat','uma')
Ix=0;Iy=0;Iz=0;
for p=1:sum(POSCAR.n_chemicals)
    Ix=Ix+POSCAR.mass(p)*uma*(norm(POSCAR.positions(p,[2 3]))*1e-10)^2;
    Iy=Iy+POSCAR.mass(p)*uma*(norm(POSCAR.positions(p,[1 3]))*1e-10)^2;
    Iz=Iz+POSCAR.mass(p)*uma*(norm(POSCAR.positions(p,[1 2]))*1e-10)^2;
end
end