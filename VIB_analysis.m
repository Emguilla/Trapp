function VIB=VIB_analysis(path)
%==================================================================================================================================%
% VIB_analysis.m: Data extraction and diagonalisation of the (multiple) phonon calculations (v0.3.1)
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
for p=1:n_TS
    % Reading of the POSCAR to which are applied finite differences
    VIB(p).POSCAR=readPOSCAR([path,ldir(p).name,'/POSCAR']);
    if n_TS>1
        VIB(p).POSCAR.Title=ldir(p).name;
    end
    % Check existence of a rotor substitution
    if exist([path,ldir(p).name,'/VibRotor.dat'],'file')
        Data=load([path,ldir(p).name,'/VibRotor.dat']);
        if Data(1)==1
            VIB(p).Rotor.mode='free';
        elseif Data(1)==2
            VIB(p).Rotor.mode='hindered';
        end
        VIB(p).Rotor.Vib_idx=Data(2);
        VIB(p).Rotor.sigma=Data(3);
        VIB(p).Rotor.I=Data(4);
        VIB(p).Rotor.RotEnergyBarrier=Data(5);
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
    VIB.L=L;
    [lambda,idx]=sort(diag(L));
    VIB.idx=idx;
    freq_Hz = sqrt(lambda)/(2*pi);
    VIB(p).DYNMAT=DYNMAT;
    % Removal of negligible frequencies
    VIB(p).nu=freq_Hz.*(abs(freq_Hz)>1e8);
%========================================================%
% Work in progress: extraction of displacement vectors
%========================================================%
    V=V(:,idx)./sqrt(repelem(DYNMAT.POSCAR.mass(:),3));
    POSCAR=readPOSCAR([ldir(p).name,'/Rotation/*/POSCAR']);
    nframes=36;
    amplitude=1;
    mode_idx=2;
    for vibidx=1:Natoms*3%mode_idx
        disp(vibidx)
        Motion_Vector=V(:,vibidx)/norm(V(:,vibidx));
        Motion_Vector=reshape(Motion_Vector,3,Natoms)';
        for q=1:nframes
            phase=2*pi*(q-1)/nframes;
            VIB(p).modes(vibidx,q)=DYNMAT.POSCAR;
            VIB(p).modes(vibidx,q).positions=DYNMAT.POSCAR.positions+sin(phase)*Motion_Vector;
        end
        VIB(p).Motion_Vector=Motion_Vector;

        VIB(p).proj(vibidx)=max(mode_projection(POSCAR,VIB(p).modes(vibidx,:)),[],'all');
    end
%========================================================%
%========================================================%
%========================================================%
    % Storage
    VIB(p).H=H; % Hessian Matrix
    VIB(p).D=D; % Dynamical Matrix
    VIB(p).V=V; % Un-mass-weighted eigenvector matrix
    % Check the presence of only 3 translation, zero-frequency modes
    if sum(VIB(p).nu==0)~=3
        error('Severe problem found: there are fewer or more than three zero-frequency vibrational mode!')
    end
    % Add a timestamp for archiving
    VIB(p).date=datetime('today');
end

% The vibration structure is stored in a ".mat" file for convenience
%save([path(1:end-1),'.mat'],'VIB')




%========================================================%
% plot data
%========================================================%
% figure;
% imagesc(3*n,3*n,log10(abs(H)))
% pbaspect([1 1 1])
% xticks(0:3:Natoms*3)
% yticks(0:3:Natoms*3)
% xticklabels('')
% xtl=xticklabels();
% xtl{2}='1';
% xtl{end}=num2str(Natoms);
% xticklabels(xtl);
% yticklabels('')
% ytl=yticklabels();
% ytl{2}='1';
% ytl{end}=num2str(Natoms);
% yticklabels(ytl);
% grid on
% set(gca,'YDir','normal');
% set(gca,'XTickLabelRotation',0);
% set(gca,'fontsize',12,'fontname','cambria math')
% 
% figure;
% imagesc(n,n,log10(Atomic))
% hold on
% pbaspect([1 1 1])
% if isstruct(DYNMAT.POSCAR)
%     for u=1:length(DYNMAT.POSCAR.n_chemicals)
%         nticks(u)=sum(DYNMAT.POSCAR.n_chemicals(1:u));
%         if u>1
%             plot([nticks(u-1)+0.5 nticks(u-1)+0.5],[0 Natoms],'w--')
%             plot([0 Natoms],[nticks(u-1)+0.5 nticks(u-1)+0.5],'w--')
%         end
%     end
% end
% set(gca,'YDir','normal');
% set(gca,'XTickLabelRotation',0);
% set(gca,'fontsize',12,'fontname','cambria math')

