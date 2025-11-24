function [ds,ds_vec,ds_xred]=ds_POSCAR(P1,P2)
%==================================================================================================================================%
% ds_POSCAR.m:  Calculation of the distances between two configurations of the same system (v0.3)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.2 (31/08/2025) - Simplification of the boundary-crossing check using direct coordinates
%       author: EYG
%   version 0.2.1 (01/09/2025) - Add "ds_xred" as an output and removed ds_norm
%       contrib: EYG
%   version 0.3 (02/09/2025) - change the update of ds_vec to account for non-orthorombic cells
%       author: EYG
%==================================================================================================================================%
% args:
%   P1, P2: POSCAR structures
%==================================================================================================================================%
% extraction of the cell geometry and cartesian coordinate of the atoms
a=P1.vec(1,:);
b=P1.vec(2,:);
c=P1.vec(3,:);
ds_vec=P2.positions-P1.positions;
ds_xred=P2.xred-P1.xred;
% since some atoms can "hop" through the boundary, the distance between the two p-ith atoms of P1 and P2 is determined as the
% shortest distance possible accounting for the periodic copies of the atom from P2. /!\ The side-effect is that if the reaction
% investigated involves an atom travelling across the whole cell, it could appear as not moving between P1 and P2 /!\
% (one explanation might be that your supercell is not large enough, if you even use a supercell at all ...)
boundary_hopping=(ones(size(ds_xred)).*(ds_xred<-0.5)-ones(size(ds_xred)).*(ds_xred>0.5));
ds_xred=ds_xred+boundary_hopping;
ds_vec=ds_vec+boundary_hopping(:,1)*a+boundary_hopping(:,2)*b+boundary_hopping(:,3)*c;
% The norm of the vector between the position of atom p in the P1 and P2 structures is also stored
% The scalar distance between the two geometries is defined using a L2-norm
ds=sqrt(sum(vecnorm(ds_vec').^2));
end