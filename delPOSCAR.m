function POSCAR=delPOSCAR(POSCAR,k)
%==================================================================================================================================%
% delPOSCAR.m:  Deletion of atoms #k in a POSCAR structure (v0.3)
%==================================================================================================================================%
% Version history:
%   version 0.1 (20/08/2025) - Creation
%       author: EYG
%   version 0.2 (20/08/2025) - Add the removal of the values of the "mass" field of the POSCAR structure
%       author: EYG
%   version 0.3 (25/09/2025) - Modification of the n_chemicals removal to include the possibility to have an empty POSCAR structure
%       author: EYG             + Add the removal of the constraint field of the POSCAR structure
%==================================================================================================================================%
% args:
%   POSCAR: POSCAR structure
%   k:      Atom numbers in POSCAR structure to be removed (can be an array)
%==================================================================================================================================%
% loop over the atoms to be removed to change the number of each species in the "n_chemicals" field of the POSCAR structure and
% check if this number goes to zero after the removal. If so, the corresponding string in the "chemical" field of the POSCAR
% structure is removed.
for pk=1:length(k)
    tmp_symbol=POSCAR.symbols(k(pk));
    for p=1:length(POSCAR.chemicals)
        if strcmpi(tmp_symbol,POSCAR.chemicals(p))
            POSCAR.n_chemicals(p)=POSCAR.n_chemicals(p)-1;
        end
    end
end
remove_chem=[];
for p=1:length(POSCAR.chemicals)
    if POSCAR.n_chemicals(p)==0
        remove_chem=[remove_chem p];
    end
end
if length(remove_chem)>=1
    POSCAR.chemicals(remove_chem)=[];
    POSCAR.n_chemicals(remove_chem)=[];
end
% Actual deletion of the atoms in the POSCAR structure
POSCAR.positions(k,:)=[];
POSCAR.symbols(k)=[];
POSCAR.xred(k,:)=[];
POSCAR.mass(k)=[];
POSCAR.Z(k)=[];
POSCAR.constraint(k,:)=[];
end