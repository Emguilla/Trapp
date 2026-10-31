function VIB=VIB_analysis(path)
%==================================================================================================================================%
% VIB_analysis.m: Data extraction and diagonalisation of the (multiple) phonon calculations (v0.4)
%==================================================================================================================================%
% Version history:
%   version 0.1 (03/07/2026) - Creation
%       author: EYG
%   version 0.2 (06/09/2026) - The output is now a structure array instead of a structure containing arrays.
%       author: EYG
%   version 0.3 (12/10/2026) - Overhaul of the code. The diagonalisation is now performed within the toolbox instead of relying on 
%       author: EYG             HIVE.
%   version 0.3.1 (12/10/2026) - Modification to accomodate the changes to "ASR_diagonal" and "ASR_Lagrange" scripts.
%       contrib: EYG
%   version 0.4 (31/10/2026) - Rotor data is now loaded from the "Rotor" field of the VIB structure, instead of the "VibRotor.dat"
%       author: EYG             file. The displacement associated with each phonon is then projected onto the rotational motion to
%                               determine which frequency should be substituted by the free or hindered rotation frequency.
%==================================================================================================================================%
% args:
%   path:   Location of the directory where the phonon calculation are located
%==================================================================================================================================%
load('constant_fund.mat')

% Check that "path" ends with a slash or backslash
if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end

% If the directory contains a DYNMAT file, it means it is a tight transition state or a groundstate calculation
if exist([path,'DYNMAT'],'file')
    ldir.name='.';
% Otherwise, it is a loose transition state and contains multiple directories, one per image.
else
    ldir=subdir_ordering(path);
end

n_TS=length(ldir);
if n_TS>1
    nbytes = fprintf('Processing image 0 out of %d', n_TS);
end
for p=1:n_TS
    if n_TS>1
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('Processing image %d out of %d\n', p, n_TS);
    end
    % Reading of the POSCAR to which are applied finite differences
    VIB(p).POSCAR=readPOSCAR([path,ldir(p).name,'/POSCAR']);
    if n_TS>1
        VIB(p).POSCAR.Title=ldir(p).name;
    end
    % Check existence of a rotor substitution
    if exist([path,ldir(p).name,'/Rotation/ROTOR.mat'],'file')
        load([path,ldir(p).name,'/Rotation/ROTOR.mat'])
        VIB(p).Rotor=ROTOR;
    elseif exist([path,ldir(p).name,'/Rotation'],'dir')
        warning('Rotation data not processed yet, ignoring rotor substitution')
        VIB(p).Rotor=NaN;
    else
        VIB(p).Rotor=NaN;
    end
    % Only the energy of the groundstate is relevant here
    Energy=readEnergy([path,ldir(p).name,'/'],'save',true);
    VIB(p).E=Energy(1);
    %Reading of the DYNMAT
    DYNMAT=readDYNMAT([path,ldir(p).name,'/DYNMAT']);
    Natoms=DYNMAT.Natoms;
    masses=DYNMAT.POSCAR.mass;
    Atomic=NaN(Natoms,Natoms);
    n=1:Natoms;
    % Reconstruction of the Hessian matrix (and Atomic Hessian)
    for u=n
        for v=n
            tmp_d=0;
            for alpha=1:3
                for beta=1:3
                    H(3*(u-1)+alpha,3*(v-1)+beta)=-DYNMAT.dF{u,alpha}(v,beta)/norm(DYNMAT.dh{u,alpha});
                end
                tmp_d=tmp_d+norm(DYNMAT.dF{u,alpha}(v,:))^2;
            end
            % Atomic Hessian matrix is the taken as the norm of the 3x3 interaction matrix between atoms u and v
            Atomic(u,v)=sqrt(tmp_d);
        end
    end
    % Symmetrisation of the Hessian matrix
    H=(H'+H)/2;
    
    % Enforcing the Acoustic Sum Rule (ASR) wrt the diagonal term
    H=ASR_diagonal(H);
    % Enforcing the ASR using lagrange multiplier (seems to reduce leaking from the TS mode, idk why yet ... tbd)
    H=ASR_Lagrange(H);

    for u=1:Natoms
        s(u)=0;
        for v=1:Natoms
            for alpha=1:3
                for beta=1:3
                    s(u)=s(u)+H(3*(u-1)+alpha,3*(v-1)+beta);
                end
            end
        end
    end
    if max(s)>1e-8
        disp('Acoustic sum')
        disp(s)
        error('Oops, something went wrong with the ASR ...')
    end
    
    % Conversion from eV/A^2 to J/m^2
    H=H*e/1e-20;
    % Mass matrix
    M=repelem(masses,3)*uma;
    % Calculation of the mass-weighted, Dynamical matrix
    D=H./sqrt(M'*M);
    % Diagonalisation
    [V,L]=eig(D);
    % Extraction of eigenvalues
    [lambda,idx]=sort(diag(L));
    freq_Hz = sqrt(lambda)/(2*pi);
    VIB(p).DYNMAT=DYNMAT;
    % Removal of negligible frequencies
    VIB(p).nu=freq_Hz.*(abs(freq_Hz)>1e8);
    % Hessian Matrix
    VIB(p).H=H;
    % Dynamical Matrix
    VIB(p).D=D;
    % Un-mass-weighted eigenvector matrix
    V=V(:,idx)./sqrt(repelem(DYNMAT.POSCAR.mass(:),3));
    VIB(p).V=V;
    
    % Extraction of the displacement associated with the rotation and with each phonon. Each phonon displacement is then projected
    % to find the best match for the rotation and determine which frequency to substitute with the rotation-like frequency of the 
    % hindered or free rotor.
    if isstruct(VIB(p).Rotor)
        Rotor_0=rotn(VIB(p).Rotor.POSCAR.positions(VIB(p).Rotor.INPUT.idx_origin,:),VIB(p).Rotor.INPUT.n,5,subPOSCAR(VIB(p).Rotor.POSCAR,VIB(p).Rotor.INPUT.idx_rotor));
        CONTCAR=VIB(p).Rotor.POSCAR;
        CONTCAR.positions(VIB(p).Rotor.INPUT.idx_rotor,:)=Rotor_0.positions;
        ds_rot=CONTCAR.positions-VIB(p).Rotor.POSCAR.positions;
        ds_rot=ds_rot/norm(ds_rot);
        for vibidx=1:Natoms*3
            [ds_vib,~]=PhononMotion(V,vibidx);
            VIB(p).Rotor.proj(vibidx)=abs(dot(ds_rot(:),ds_vib(:)))/norm(ds_rot(:))/norm(ds_vib(:));
        end
        [~,VIB(p).Rotor.vibidx]=max(VIB(p).Rotor.proj);
    end
    % Check the presence of only 3 translation, zero-frequency modes
    if sum(VIB(p).nu==0)~=3
        error('Severe problem found: there are fewer or more than three zero-frequency vibrational mode!')
    end
    % Add a timestamp for archiving
    VIB(p).date=datetime('today');
end

% The vibration structure is stored in a ".mat" file for convenience, without the H, D, V and DYNMAT fields (they are too massive 
% to be stored systematically). 
tmp=VIB;
VIB=rmfield(VIB,{'DYNMAT','H','D','V'});
save([path(1:end-1),'.mat'],'VIB')
VIB=tmp;