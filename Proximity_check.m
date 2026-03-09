function [POSCAR,bonds]=Proximity_check(POSCAR,varargin)
%==================================================================================================================================%
% Proximity_check.m:    Check distances between all possible pairs of atoms to identify and optionally modify their positions (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (09/03/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   POSCAR:             POSCAR structure, or path+filename of a POSCAR file
%   opt. args:          'min', followed by the ratio wrt the maximum bond distance allowed between each pair of atoms. 
%                           (default: 0.82)
%                       'pushback', followed by true or false to specify whether pairs of too close atoms should be repelled.
%                           (default: false)
%                       'verbose', followed by true or false to display information about the rendering (timing and number of atoms)
%                           (default: false)
%==================================================================================================================================%
load('ptable.mat')
% Initialisation of the default parameters and handling of the optional arguments
pushback=false;
verbose=false;
DBOND_MIN_ratio=ptable*bond_length;

% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'min'
                DBOND_MIN_ratio=varargin{p+1};
                DBOND_MIN=DBOND_MIN_ratio*DBOND_MIN;
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
    for pZ=unique(Z)
        for qZ=unique(Z)
            if Z(pairs(ks(p_ks),1))==pZ&&Z(pairs(ks(p_ks),2))==qZ
                if ds(ks(p_ks))<DBOND_MIN(pZ,qZ)
                    ks_n(k)=ks(p_ks);
                    DBOND_update(k)=DBOND_MIN(pZ,qZ);
                    k=k+1;
                end
            end
        end
    end
end

% Optional displacement of atoms. The function is called recursively to ensure that displaced atoms are not displaced too close to
% other atoms.
if pushback
    is_moved=false(Natoms,1);
    for p=1:length(ks_n)
        [~,Zmin]=min([Z(pairs(ks_n(p),:))]);
        if ~is_moved(pairs(ks_n(p),Zmin))
            r=POSCAR.positions(pairs(ks_n(p),:),:);
            r1=r(Zmin,:);
            r(Zmin,:)=[];
            r2=r;
            dr=(r1-r2)/norm((r1-r2))*DBOND_update(p);
            if norm(r1-r2-dr)>1e-3
                POSCAR.positions(pairs(ks_n(p),Zmin),:)=r2+dr;
                if verbose
                    disp(['changing distance between atoms ',num2str(pairs(ks_n(p),1)),' and ',num2str(pairs(ks_n(p),2)),' by moving atom ',num2str(pairs(ks_n(p),Zmin))])
                end
                is_moved(pairs(ks_n(p),Zmin))=true;
            end
        end
    end
    if any(is_moved)
        for p=1:Natoms
            POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
        end
        [POSCAR,~]=Proximity_check(POSCAR,'min',DBOND_MIN_ratio,'pushback',true);
    end
end

bonds.idx=pairs(ks_n,:);
bonds.ds=ds(ks_n);

% end