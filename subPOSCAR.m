function POSCAR=subPOSCAR(POSCAR,idx,varargin)
%==================================================================================================================================%
% subPOSCAR.m:  Creation of a subPOSCAR structure containing only the atoms #idx from a POSCAR structure (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (20/08/2025) - Creation
%       author: EYG
%   version 0.2 (19/02/2026) - Add an optional argument to define the subPOSCAR with respect to the positions of the atoms in terms
%       author: EYG             of the unit cell vectors. This new optional argument, when activated, superseeds the idx argument.
%==================================================================================================================================%
% args:
%   POSCAR:     POSCAR structure
%   idx:        Atom numbers in POSCAR structure to be removed (can be an array)
%   opt. args:  'Lim', followed by the limits (in direct coordinates) of the sub/supercell in a 3x2 array.
%                   (default: not defined)
%==================================================================================================================================%
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'lim'
                idx=0; % the use of the "Lim" optional argument superseeds the idx argument
                rendering=true;
                aLim=varargin{p+1}(1,:);
                bLim=varargin{p+1}(2,:);
                cLim=varargin{p+1}(3,:);
            otherwise
                warning(['Unknown argument: "',varargin{p},'"'])
        end
    end
end
if idx~=0
    bin_idx=1:sum(POSCAR.n_chemicals);
    bin_idx(unique(idx))=[];
    POSCAR=delPOSCAR(POSCAR,bin_idx);
else
    POSCAR_init=POSCAR;
    % Get the limit of the supercell
    for p=floor(aLim(1)):ceil(aLim(2))-1
        for q=floor(bLim(1)):ceil(bLim(2))-1
            for r=floor(cLim(1)):ceil(cLim(2))-1
                if ~(p==0&&q==0&&r==0)
                    POSCAR_tmp=POSCAR_init;
                    POSCAR_tmp.positions=POSCAR_tmp.positions+p*POSCAR.vec(1,:)+q*POSCAR.vec(2,:)+r*POSCAR.vec(3,:);
                    POSCAR_tmp.xred=POSCAR.xred+p+q+r;
                    POSCAR=appendPOSCAR(POSCAR,POSCAR_tmp);
                end
            end
        end
    end
    % Deletion of all the atoms that are beyond the boundaries set by XLim, YLim and ZLim
    count=0;
    max_count=sum(POSCAR.n_chemicals);
    k=1;
    while count<max_count
        if POSCAR.xred(k,1)<aLim(1)||POSCAR.xred(k,1)>aLim(2)
            POSCAR=delPOSCAR(POSCAR,k);
            k=k-1;
        elseif POSCAR.xred(k,2)<bLim(1)||POSCAR.xred(k,2)>bLim(2)
            POSCAR=delPOSCAR(POSCAR,k);
            k=k-1;
        elseif POSCAR.xred(k,3)<cLim(1)||POSCAR.xred(k,3)>cLim(2)
            POSCAR=delPOSCAR(POSCAR,k);
            k=k-1;
        end
        k=k+1;
        count=count+1;
    end
end
end