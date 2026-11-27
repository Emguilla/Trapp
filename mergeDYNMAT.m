function DYNMAT=mergeDYNMAT(D1,D2,varargin)
%==================================================================================================================================%
% mergeDYNMAT.m:    Merge DYNMAT structures (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (27/11/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   D1, D2: DYNMAT structures
% opt. args:    'force_merge', followed by true or false to determine whether to force the merge of the DYNMAT in case duplicate 
%                   DOF are found, and use first argument as reference.
%                   (default: false)
%==================================================================================================================================%

% Initialisation of default parameters
force_merge=false;

if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'force_merge'
                force_merge=varargin{p+1};
            otherwise
                error(['Unknown argument: ',varargin{p}])
        end
    end
end

% check compatibility between POSCAR structures
[SIGERR,msg_err]=compatPOSCAR(D1.POSCAR,D2.POSCAR);
if SIGERR~=0
    error(msg_err);
end

% check completeness of the merged DYNMAT structure
if any((D1.mask(:)+D2.mask(:))>1)
    if force_merge
        warning('Duplicate DOF were found. The finite-difference from the first DYNMAT argument will be considered as reference.')
    else
        error('Duplicate DOF were found. The finite-difference from the first DYNMAT argument will be considered as reference.')
    end
end
if any((D1.mask(:)+D2.mask(:))<1)
    warning('Some DOF are missing. Resulting merged DYNMAT will be incomplete.')
end

% Use D2 as a starting point, and fill out/replace finite-difference values everywhere 'D1.mask' is true
DYNMAT=D2;
for p=1:Natoms
    for q=1:3
        if D1.mask(p,q)
            DYNMAT.mask(p,q)=true;
            DYNMAT.dh{p,q}=D1.dh{p,q};
            DYNMAT.dF{p,q}=D1.dF{p,q};
        end
    end
end


