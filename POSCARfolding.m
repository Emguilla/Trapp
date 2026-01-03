function POSCAR=POSCARfolding(POSCAR)
%==================================================================================================================================%
% POSCARfolding.m: Shift of stray atoms back into the cell according to periodic boundary conditions (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (03/01/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   POSCAR:         POSCAR structure, or path+filename of a POSCAR file
%==================================================================================================================================%
for p=1:sum(POSCAR.n_chemicals)
    update=false;
    for q=1:3
        % if the direct coordinate is outside the interval [0;1], the atoms is shifted by one lattice cell until it is in the cell.
        while POSCAR.xred(p,q)<0
            update=true;
            POSCAR.xred(p,q)=POSCAR.xred(p,q)+1;
        end
        while POSCAR.xred(p,q)>=1
            update=true;
            POSCAR.xred(p,q)=POSCAR.xred(p,q)-1;
        end
    end
    % if any of the direct coordinates of atom p have been modified, the cartesian coordinates are updated as well. 
    if update
        POSCAR.positions(p,:)=POSCAR.acell*POSCAR.vec'*POSCAR.xred(p,:)';
    end
end
