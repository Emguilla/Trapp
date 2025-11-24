function POSCAR=subPOSCAR(POSCAR,idx)
%==================================================================================================================================%
% subPOSCAR.m:  Creation of a subPOSCAR structure containing only the atoms #idx from a POSCAR structure (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (20/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   POSCAR: POSCAR structure
%   idx:    Atom numbers in POSCAR structure to be removed (can be an array)
%==================================================================================================================================%
bin_idx=1:sum(POSCAR.n_chemicals);
bin_idx(unique(idx))=[];
POSCAR=delPOSCAR(POSCAR,bin_idx);
end