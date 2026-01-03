function CONTCAR=sortPOSCAR(POSCAR,varargin)
%==================================================================================================================================%
% sortPOSCAR.m: sort atoms in POSCAR according to their atomic number or in a specified order (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (03/01/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   POSCAR:     POSCAR structure
%   opt. args:  'idx', followed by an array of integer providing the order in which the atoms must be re-arranged.
%                   (default: atoms are ordered according to their atomic numbers)
%==================================================================================================================================%
idx_sort=false;
% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'idx'
                idx_sort=true;
                idx=varargin{p+1};
        end
    end
end
if ~idx_sort
    [~,idx]=sort(POSCAR.Z);
end
% update of all the fields of the POSCAR structure
CONTCAR=POSCAR;
CONTCAR.positions=POSCAR.positions(idx,:);
CONTCAR.constraint=POSCAR.constraint(idx,:);
CONTCAR.xred=POSCAR.xred(idx,:);
CONTCAR.symbols=POSCAR.symbols(idx);
CONTCAR.mass=POSCAR.mass(idx);
CONTCAR.Z=POSCAR.Z(idx);

