function [POSCAR,bonds]=Proximity_check(POSCAR,varargin)
%==================================================================================================================================%
% Proximity_check.m:    Check distances between all possible pairs of atoms to identify and optionally modify their positions 
%                       (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (09/03/2026) - Creation
%       author: EYG
%   version 0.1.1 (12/03/2026) - Fix bug happening with the activation of the optional argument 'min'.
%       author: EYG
%   version 0.2 (30/04/2026) - Removal of the recursive call for efficiency. The distances are only checked for different atoms, 
%       author: EYG             this may change in future updates. Added an 'absmin' optional argument.
%==================================================================================================================================%
% args:
%   POSCAR:             POSCAR structure, or path+filename of a POSCAR file
%   opt. args:          'min', followed by the ratio wrt the maximum bond distance allowed between each pair of atoms. 
%                           (default: 0.82)
%                       'absmin', followed by the minimum distance between atoms. Does not apply to same-type atoms.
%                           (default: 1)
%                       'pushback', followed by true or false to specify whether pairs of too close atoms should be repelled.
%                           (default: false)
%                       'verbose', followed by true or false to display information about the rendering (timing and number of atoms)
%                           (default: false)
%==================================================================================================================================%
load('ptable.mat')
% Initialisation of the default parameters and handling of the optional arguments
pushback=false;
verbose=false;
DBOND_MIN=ptable.bond_length;
DBOND_MIN_abs=1;

% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'min'
                DBOND_MIN_ratio=varargin{p+1};
                DBOND_MIN=DBOND_MIN_ratio*DBOND_MIN;
            case 'absmin'
                DBOND_MIN_abs=varargin{p+1};
                DBOND_MIN=DBOND_MIN_abs;
            case 'pushback'
                pushback=varargin{p+1};
            case 'verbose'
                verbose=varargin{p+1};
        end
    end
end

Natoms=sum(POSCAR.n_chemicals);
% get all combinations of atoms
pairs = nchoosek(1:Natoms,2);
% get all distance between pairs, then restrict the list to the distance below a threshold (here the maximum distance in the DBOND
% matrix). This is only a first clean-up of all the possible pairs.
ds=sqrt(sum((POSCAR.positions(pairs(:,1),:)-POSCAR.positions(pairs(:,2),:)).^2,2));
ks=find(ds<max(DBOND_MIN(:)))';
Z=POSCAR.Z;
% Browsing through the list, check the symbols and keep the pair only if the distance is lower than the one defined in the DBOND
% matrix.
k=1;
ks_n=[];
for p_ks=1:length(ks)
    if ds(ks(p_ks))<DBOND_MIN(Z(pairs(ks(p_ks),1)),Z(pairs(ks(p_ks),2)))
        if verbose
            disp(['Atoms ',num2str(pairs(ks(p_ks),1)),' and ',num2str(pairs(ks(p_ks),2)),' are too close!'])
        end
        if pushback
            if Z(pairs(ks(p_ks),1))<Z(pairs(ks(p_ks),2))
                idx_A=1;
                idx_B=2;
            elseif Z(pairs(ks(p_ks),1))>Z(pairs(ks(p_ks),2))
                idx_A=2;
                idx_B=1;
            end
            pos_A=POSCAR.positions(pairs(ks(p_ks),idx_A));
            pos_B=POSCAR.positions(pairs(ks(p_ks),idx_B));
            d=(pos_A-pos_B)/norm(pos_A-pos_B);
            min_d=DBOND_MIN(Z(pairs(ks(p_ks),1)),Z(pairs(ks(p_ks),2)));
            POSCAR.positions(pairs(ks(p_ks),idx_A))=POSCAR.positions(pairs(ks(p_ks),idx_B))+min_d*d;
        end
    end
end
for p=1:Natoms
    POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
end

bonds.idx=pairs(ks_n,:);
bonds.ds=ds(ks_n);

% end