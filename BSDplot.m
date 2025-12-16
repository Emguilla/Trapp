function BSDplot(EIGENVAL,DOSCAR,varargin)
if ~isstruct(EIGENVAL)
    EIGENVAL=readEIGENVAL(EIGENVAL);
end
if ~isstruct(DOSCAR)
    DOSCAR=readDOSCAR(DOSCAR);
end
pltBS=true;
pltDOS=true;
filename='BSD';
symmetry=NaN;
issymmetry=false;
POSCAR=NaN;
KPOINTS=NaN;
if exist('varargin')
    for p=1:2:length(varargin)
        switch varargin{p}
            case 'symmetry'
                symmetry=varargin{p+1};
            case 'lm-projection'
                lm_decomp=varargin{p+1};
            case 'KPOINTS'
                KPOINTS=varargin{p+1};
                if ~isstruct(KPOINTS)
                    KPOINTS=readKPOINTS(KPOINTS);
                end
            case 'POSCAR'
                POSCAR=varargin{p+1};
                if ~isstruct(POSCAR)
                    POSCAR=readPOSCAR(POSCAR);
                end
            case 'save'
                filename=varargin{p+1};
            case 'BS'
                pltBS=varargin{p+1};
            case 'DOS'
                pltDOS=varargin{p+1};
        end
    end
end
while ~issymmetry
    issymmetry=true;
    switch lower(symmetry)
        case {1,'cub','cubic'}
            symmetry=1;% Cubic
            % ! Gamma-point %! M-point %! R-point %! X-point
            HS=[0 0 0; 1/2 1/2 0; 1/2 1/2 1/2; 0 1/2 0];
            HS_name={'\Gamma','M','R','X'};
        case {2,'fcc','face-centered cubic'}
            symmetry=2;% Face-centered cubic
            % ! Gamma-point %! X-point %! W-point %! L-point %! K-point %! U-point
            HS=[0 0 0; 3/8 3/8 3/4; 1/2 1/2 1/2; 5/8  1/4 5/8; 1/2 1/4 3/4; 1/2 0 1/2];
            HS_name={'\Gamma','K','L','U','W','X'};
        case {3,'bcc','body-centered cubic'}
            symmetry=3;% Body-centered cubic
            %! Gamma-point %! H-point %! P-point %! N-point
            HS=[0 0 0; 1/2 -1/2 1/2; 1/4 1/4 1/4; 0 0 1/2];
            HS_name={'\Gamma','H','P','N'};
        case {4,'tet','tetragonal'}
            symmetry=4;% Tetragonal
            %! Gamma-point %! A-point %! M-point %! R-point %! X-point %! Z-point
            HS=[0 0 0; 1/2 1/2 1/2; 1/2 1/2 0; 0 1/2 1/2; 0 1/2 0; 0 0 1/2];
            HS_name={'\Gamma','A','M','R','X','Z'};
        case {5,'bct','body-centered tetragonal'}
            symmetry=5;% Body-centered tetragonal
            BCT1=false;
            BCT2=false;
            if ~isnan(POSCAR)
                lattice=POSCAR.vec*POSCAR.acell;
                if lattice(1,1)==-lattice(2,1)&&lattice(1,1)==-lattice(1,3)...
                        &&lattice(2,2)==-lattice(1,2)&&lattice(2,2)==-lattice(3,2)...
                        &&lattice(3,3)==-lattice(1,3)&&lattice(3,3)==-lattice(2,3)...
                        &&lattice(1,2)==lattice(2,1) % primitive cell case
                    a=abs(lattice(1,1)*2);
                    c=abs(lattice(3,3)*2);
                elseif all((lattice-diag(lattice).*eye(3,3))==0) % conventional cell case
                    a=lattice(1,1);
                    c=lattice(3,3);
                else
                    error('Unknown lattice vectors. Check POSCAR.')
                end
                if c<a
                    BCT1=true;
                    eta=(1+(c^2)/(a^2))/4;
                    %! Gamma-point %! M-point %! N-point %! P-point %! X-point %! Z-point %! Z1-point
                    HS=[0 0 0; -1/2 1/2 1/2; 0 1/2 0; 1/4 1/4 1/4; 0 0 1/2; eta eta -eta; -eta 1-eta eta];
                    HS_name={'\Gamma','M','N','P','X','Z','Z_1'};
                elseif a>c
                    BCT2=true;
                    eta=(1+(a^2)/(c^2))/4;
                    zeta=(a^2)/(2*c^2);
                    %! Gamma-point %! N-point %! P-point %! Sigma-point %! Sigma1-point %! X-point %! Y-point %! Y1-point %! Z-point
                    HS=[0 0 0; 0 1/2 0; 1/4 1/4 1/4; -eta eta eta; eta 1-eta -eta; 0 0 1/2; -zeta zeta 1/2; 1/2 1/2 -zeta; 1/2 1/2 -1/2];
                    HS_name={'\Gamma','N','P','\Sigma','\Sigma_1','X','Y','Y_1','Z'};
                else
                    warning('You choose a BCT symmetry, but the POSCAR indicates a cubic symmetry. Defaulting to "cub" option.')
                    is_symmetry=false;
                    symmetry=1;
                end
            else
                error('Body-centered tetragonal bandstructure plot requires to set the optional argument POSCAR!')
            end
        case {6,'orc','orthorhombic'}
            symmetry=6;% Orthorhombic
            %! Gamma-point %! M-point %! N-point %! P-point %! X-point %! Z-point %! Z1-point
            HS=[0 0 0; 1/2 1/2 1/2; 1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2; 1/2 0 0; 0 1/2 0; 0 0 1/2];
            HS_name={'\Gamma','N','P','\Sigma','\Sigma_1','X','Y','Y_1','Z'};
        case {7,'orcf','face-centered orthorhombic'}
            symmetry=7;% Face-centered orthorhombic
            ORCF1=false;
            ORCF2=false;
            ORCF3=false;
            if ~isnan(POSCAR)
                lattice=POSCAR.vec*POSCAR.acell;
                if all(diag(lattice)==0) % primitive cell case
                    a=lattice(2,1)*2;
                    b=lattice(1,2)*2;
                    c=lattice(1,3)*2;
                elseif all((lattice-diag(lattice).*eye(3,3))==0) % conventional cell case
                    a=lattice(1,1);
                    b=lattice(2,2);
                    c=lattice(3,3);
                else
                    error('Unknown lattice vectors. Check POSCAR.')
                end
                if ~(c>b&&b>a)
                    error('There seems to be a mistake in the ordering of the lattice vectors. Check POSCAR.')
                end
                if (1/a^2)>(1/b^2+1/c^2)
                    ORCF1=true;
                elseif (1/a^2)<(1/b^2+1/c^2)
                    ORCF2=true;
                elseif (1/a^2)==(1/b^2+1/c^2)
                    ORCF3=true;
                end
                if ORCF1||ORCF3
                    zeta=(1+a^2/b^2-a^2/c^2)/4;
                    eta=(1+a^2/b^2+a^2/c^2)/4;
                    %! Gamma-point %! A-point %! A1-point %! L-point %! T-point %! X-point %! X1-point %! Y-point %! Z-point  
                    HS=[0 0 0; 1/2 1/2+zeta zeta; 1/2 1/2-zeta 1-zeta; 1/2 1/2 1/2;1 1/2 1/2; 0 eta eta; 1 1-eta 1-eta; 1/2 0 1/2; 1/2 1/2 0];
                    HS_name={'\Gamma','A','A_1','L','T','X','X_1','Y','Z'};
                elseif ORCF2
                    eta=(1+a^2/b^2-a^2/c^2)/4;
                    delta=(1+b^2/a^2-b^2/c^2)/4;
                    phi=(1+c^2/b^2-c^2/a^2)/4;
                    %! Gamma-point %! C-point %! C1-point %! D-point %! D1-point %! L-point %! H-point %! H1-point %! X-point %! Y-point %! Z-point  
                    HS=[0 0 0; 1/2 1/2-eta 1-eta; 1/2 1/2+eta eta; 1/2-delta 1/2 1-delta; 1/2+delta 1/2 delta; 1/2 1/2 1/2; 1-phi 1/2-phi 1/2; phi 1/2+phi 1/2; 0 1/2 1/2; 1/2 0 1/2; 1/2 1/2 0];
                    HS_name={'\Gamma','C','C_1','D','D_1','L','H','H_1','X','Y','Z'};
                end
            else
                error('Face-centered orthorhombic bandstructure plot requires to set the optional argument POSCAR!')
            end
        case {8,'orci','body-centered orthorhombic'}
            symmetry=8;% Body-centered orthorhombic
            if ~isnan(POSCAR)
                lattice=POSCAR.vec*POSCAR.acell;
                if lattice(1,1)==-lattice(2,1)&&lattice(1,1)==-lattice(1,3)...
                        &&lattice(2,2)==-lattice(1,2)&&lattice(2,2)==-lattice(3,2)...
                        &&lattice(3,3)==-lattice(1,3)&&lattice(3,3)==-lattice(2,3)...
                        &&lattice(1,2)~=lattice(2,1) % primitive cell case
                    a=abs(lattice(1,1)*2);
                    b=abs(lattice(2,2)*2);
                    c=abs(lattice(3,3)*2);
                elseif all((lattice-diag(lattice).*eye(3,3))==0) % conventional cell case
                    a=lattice(1,1);
                    b=lattice(2,2);
                    c=lattice(3,3);
                else
                    error('Unknown lattice vectors. Check POSCAR.')
                end
                if ~(c>b&&b>a)
                    error('There seems to be a mistake in the ordering of the lattice vectors. Check POSCAR.')
                end
                zeta=(1+a^2/c^2)/4;
                eta=(1+b^2/c^2)/4;
                delta=(b^2-a^2)/(4*c^2);
                mu=(a^2+b^2)/(4*c^2);
                %! Gamma-point %! L-point %! L1-point %! L2-point %! R-point %! S-point %! T-point %! W-point %! X-point %! X1-point %! Y-point %! Y1-point %! Z-point  
                HS=[0 0 0; -mu mu 1/2-delta; mu -mu 1/2+delta; 1/2-delta 1/2+delta -mu; 0 1/2 0; 1/2 0 0; 0 0 1/2; 1/4 1/4 1/4;-zeta zeta zeta; zeta 1-zeta -zeta;eta -eta eta;1-eta eta -eta;1/2 1/2 -1/2];
                HS_name={'\Gamma','L','L_1','L_2','R','S','T','W','X','X_1','Y','Y_1','Z'};
            else
                error('Body-centered orthorhombic bandstructure plot requires to set the optional argument POSCAR!')
            end
        case {9,'orcc','c-centered orthorhombic'}
            symmetry=9;% C-centered orthorhombic
            if ~isnan(POSCAR)
                lattice=POSCAR.vec*POSCAR.acell;
                if all(lattice([3 6:8])==0)&&lattice(1,1)==lattice(2,1)&&lattice(1,2)==-lattice(2,2) % primitive cell case
                    a=lattice(1,1)*2;
                    b=lattice(2,2)*2;
                    c=lattice(3,3);
                elseif all((lattice-diag(lattice).*eye(3,3))==0) % conventional cell case
                    a=lattice(1,1);
                    b=lattice(2,2);
                    c=lattice(3,3);
                else
                    error('Unknown lattice vectors. Check POSCAR.')
                end
                if ~(c>b&&b>a)
                    error('There seems to be a mistake in the ordering of the lattice vectors. Check POSCAR.')
                end
                zeta=(1+a^2/c^2)/4;
                eta=(1+b^2/c^2)/4;
                delta=(b^2-a^2)/(4*c^2);
                mu=(a^2+b^2)/(4*c^2);
                %! Gamma-point %! L-point %! L1-point %! L2-point %! R-point %! S-point %! T-point %! W-point %! X-point %! X1-point %! Y-point %! Y1-point %! Z-point  
                HS=[0 0 0; -mu mu 1/2-delta; mu -mu 1/2+delta; 1/2-delta 1/2+delta -mu; 0 1/2 0; 1/2 0 0; 0 0 1/2; 1/4 1/4 1/4;-zeta zeta zeta; zeta 1-zeta -zeta;eta -eta eta;1-eta eta -eta;1/2 1/2 -1/2];
                HS_name={'\Gamma','L','L_1','L_2','R','S','T','W','X','X_1','Y','Y_1','Z'};
            else
                error('C-centered orthorhombic bandstructure plot requires to set the optional argument POSCAR!')
            end
        case {10,'hex','hexagonal'}
            symmetry=10;% Hexagonal
            %! Gamma-point %! A-point %! H-point %! K-point %! L-point %! M-point 
            HS=[0 0 0;0 0 1/2;1/3 1/3 1/2;1/3 1/3 0;1/2 0 1/2;1/2 0 0];
            HS_name={'\Gamma','A','H','K','L','M'};
        case {11,'rhl','rhombohedral'}
            symmetry=11;% Rhombohedral
            RHL1=false;
            RHL2=false;
            if ~isnan(POSCAR)
                lattice=POSCAR.vec*POSCAR.acell;
                if all(lattice(6:8)==0)
                    alpha=pi-acos((norm(lattice(1,:))^2+norm(lattice(2,:))^2-norm(lattice(1,:)+lattice(2,:))^2)/(2*norm(lattice(1,:))*norm(lattice(2,:))));
                    a=lattice(1,1)/cos(alpha/2);
                    if ~all((lattice(:)-[a*cos(alpha/2) -a*sin(alpha/2) 0;a*cos(alpha/2) a*sin(alpha/2) 0;a*cos(alpha)/cos(alpha/2) 0 a*sqrt(1-cos(alpha)^2/cos(alpha/2)^2)])==0)
                        error('Unknown lattice vectors. Check POSCAR.')
                    end
                else
                    error('Unknown lattice vectors. Check POSCAR.')
                end
                if alpha<pi/2
                    RHL1=true;
                    eta=(1+4*cos(alpha))/(2+4*os(alpha));
                    nu=3/4+eta/2;
                    %! Gamma-point %! B-point %! B1-point %! F-point %! L-point %! L1-point %! P-point %! P1-point %! P2-point %! Q-point %! X-point %! Z-point
                    HS=[0 0 0; eta 1/2 1-eta; 1/2 1-eta eta-1; 1/2 1/2 0; 1/2 0 0;0 0 -1/2;eta nu nu; 1-nu 1-nu 1-eta; nu nu eta-1;1-nu nu 0;nu 0 -nu;1/2 1/2 1/2];
                    HS_name={'\Gamma','B','B_1','F','L','L_1','P','P_1','P_2','Q','X','Z'};
                elseif alpha>pi/2
                    RHL2=true;
                    eta=1/(2*tan(alpha/2)^2);
                    nu=3/4-eta/2;
                    %! Gamma-point %! F-point %! L-point %! P-point %! P1-point %! Q-point %! Q1-point %! Z-point
                    HS=[0 0 0;1/2 -1/2 0; 1/2 0 0; 1-nu -nu 1-nu;nu nu-1 nu-1; eta eta eta;1-eta -eta -eta; 1/2 -1/2 1/2];
                    HS_name={'\Gamma','F','L','P','P1','Q','Q_1','Z'};
                else
                    warning('You choose a RHL symmetry, but the POSCAR indicates a cubic symmetry (alpha=pi/2). Defaulting to "cub" option.')
                    is_symmetry=false;
                    symmetry=1;
                end
            else
                error('Rhombohedral bandstructure plot requires to set the optional argument POSCAR!')
            end
        case {12,'mcl','monoclinic'}
            symmetry=12;% Monoclinic
            if ~isnan(POSCAR)
                lattice=POSCAR.vec*POSCAR.acell;
                if all(lattice([2:4 7:8])==0)
                    alpha=pi-acos((norm(lattice(2,:))^2+norm(lattice(3,:))^2-norm(lattice(2,:)+lattice(3,:))^2)/(2*norm(lattice(2,:))*norm(lattice(3,:))));
                    a=lattice(1,1);
                    b=lattice(2,2);
                    c=lattice(3,3)/sin(alpha);
                    if ~all((lattice(:)-[a 0 0;0 b 0;0 c*cos(alpha) c*sin(alpha)])==0)
                        error('Unknown lattice vectors. Check POSCAR.')
                    end
                else
                    error('Unknown lattice vectors. Check POSCAR.')
                end
                eta=(1-b*cos(alpha)/c)/(2*sin(alpha)^2);
                nu=1/2-eta*c*cos(alpha)/b;
                %! Gamma-point %! A-point %! C-point %! D-point %! D1-point %! E-point %! H-point %! H1-point %! H2-point %! M-point %! M1-point %! M2-point %! X-point %! Y-point %! Y1-point %! Z-point
                HS=[0 0 0; 1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2; 1/2 0 -1/2; 1/2 1/2 1/2; 0 eta 1-nu; 0 1-eta nu; 0 eta -nu; 1/2 eta 1-nu; 1/2 1-eta nu; 1/2 eta -nu; 0 1/2 0; 0 0 1/2; 0 0 -1/2; 1/2 0 0];
                HS_name={'\Gamma','A','C','D','D1','E','H','H_1','H_2','M','M_1','M_2','X','Y','Y_1','Z'};
            else
                error('Monoclinic bandstructure plot requires to set the optional argument POSCAR!')
            end
        case {13,'mclc','c-centered monoclinic'}
            symmetry=13;% C-centered monoclinic
            MCLC1=false;
            MCLC2=false;
            MCLC3=false;
            MCLC4=false;
            MCLC5=false;
            if ~isnan(POSCAR)
                lattice=POSCAR.vec*POSCAR.acell;
                if all(lattice([4 7:8])==0) % Conventional cell case
                    alpha=pi-acos((norm(lattice(2,:))^2+norm(lattice(3,:))^2-norm(lattice(2,:)+lattice(3,:))^2)/(2*norm(lattice(2,:))*norm(lattice(3,:))));
                    a=lattice(1,1);
                    b=lattice(2,2);
                    c=lattice(3,3)/sin(alpha);
                    if ~all((lattice(:)-[a 0 0;0 b 0;0 c*cos(alpha) c*sin(alpha)])==0)
                        error('Unknown lattice vectors. Check POSCAR.')
                    end
                else
                    error('Unknown lattice vectors. Check POSCAR.')
                end
                a1=[a/2 b/2 0];
                a2=[-a/2 b/2 0];
                a3=[0 c*cos(alpha) c*sin(alpha)];
                V=dot(a1,cross(a2,a3));
                b1=2*pi*cross(a2,a3)/V;
                b2=2*pi*cross(a3,a1)/V;
                b3=2*pi*cross(a1,a2)/V;
                kgam=pi-acos((norm(b1)^2+norm(b2)^2-norm(b1+b2)^2)/(2*norm(b1)*norm(b2)));
                if kgam>pi/2
                    MCLC1=true;
                elseif kgam==pi/2
                    MCLC2=true;
                elseif (b*cos(alpha)/c+b^2*sin(alpha)^2/a^2)<1
                    MCLC3=true;
                elseif (b*cos(alpha)/c+b^2*sin(alpha)^2/a^2)==1
                    MCLC4=true;
                elseif (b*cos(alpha)/c+b^2*sin(alpha)^2/a^2)>1
                    MCLC5=true;
                end
                if MCLC1||MCLC2
                    zeta=(2-b*cos(alpha)/c)/(4*sin(alpha)^2);
                    eta=1/2+2*zeta*c*cos(alpha)/b;
                    psi=3/4+a^2/(4*b^2*sin(alpha)^2);
                    phi=psi+(3/4-psi)*cos(alpha)/c;
                    %! Gamma-point %! N-point %! N1-point %! F-point %! F1-point %! F2-point %! F3-point %! I-point %! I1-point %! L-point %! M-point %! X-point %! X1-point %! X2-point %! Y-point %! Y1-point %! Z-point
                    HS=[0 0 0; 1/2 0 0; 0 -1/2 0; 1-zeta 1-zeta 1-eta; zeta zeta eta; -zeta -zeta 1-eta; 1-zeta -zeta 1-eta; phi 1-phi 1/2; 1-phi phi-1 1/2; 1/2 1/2 1/2; 1/2 0 1/2; 1-psi psi-1 0; psi 1-psi 0; psi-1 -psi 0; 1/2 1/2 0; -1/2 -1/2 0; 0 0 1/2];
                    HS_name={'\Gamma','N','N_1','F','F_1','F_2','F_3','I','I_1','L','M','X','X1','X2','Y','Y1','Z'};
                elseif MCLC3||MCLC4
                    mu=(1+b^2/a^2)/4;
                    delta=b*c*cos(alpha)/(2*a^2);
                    zeta=mu-1/4+(1-b*cos(alpha)/c)/(4*sin(alpha)^2);
                    eta=1/2+2*zeta*c*cos(alpha)/b;
                    phi=1+zeta-2*mu;
                    psi=eta-2*delta;
                    %! Gamma-point %! F-point %! F1-point %! F2-point %! H-point %! H1-point %! H2-point %! I-point %! M-point %! N-point %! N1-point %! X-point %! Y-point %! Y1-point %! Y2-point %! Y3-point %! Z-point
                    HS=[0 0 0; 1-phi 1-phi 1-psi; phi phi-1 psi; 1-phi -phi 1-psi; zeta zeta eta; 1-zeta -zeta 1-eta; -zeta -zeta 1-eta; 1/2 -1/2 1/2; 1/2 0 1/2; 1/2 0 0; 0 -1/2 0; 1/2 -1/2 0; mu mu delta; 1-mu -mu -delta; -mu -mu -delta; mu mu-1 delta; 0 0 1/2];
                    HS_name={'\Gamma','F','F_1','F_2','H','H_1','H_2','I','M','N','N_1','X','Y','Y1','Y2','Y3','Z'};
                elseif MCLC5
                    zeta=(b^2/a^2+(1-b*cos(alpha)/c)/sin(alpha)^2)/4;
                    eta=1/2+2*zeta*c*cos(alpha)/b;
                    mu=eta/2+b^2/(4*a^2)-b*c*cos(alpha)/(2*a^2);
                    nu=2*mu-zeta;
                    omega=(4*nu-1-b^2*sin(alpha)^2/a^2)*c/(2*b*cos(alpha));
                    delta=zeta*c*cos(alpha)/b+omega/2-1/4;
                    rho=1-zeta*a^2/b^2;
                    %! Gamma-point %! F-point %! F1-point %! F2-point %! H-point %! H1-point %! H2-point %! I-point %! I1-point %! L-point %! M-point %! N-point %! N1-point %! X-point %! Y-point %! Y1-point %! Y2-point %! Y3-point %! Z-point
                    HS=[0 0 0; nu nu omega; 1-nu 1-nu 1-omega; nu nu-1 omega; zeta zeta eta; 1-zeta -zeta 1-eta; -zeta -zeta 1-eta; rho 1-rho 1/2; 1-rho rho-1 1/2; 1/2 1/2 1/2; 1/2 0 1/2; 1/2 0 0; 0 -1/2 0; 1/2 -1/2 0; mu mu delta; 1-mu -mu -delta; -mu -mu -delta; mu mu-1 delta; 0 0 1/2];
                    HS_name={'\Gamma','F','F_1','F_2','H','H_1','H_2','I','I_1','L','M','N','N_1','X','Y','Y1','Y2','Y3','Z'};
                end
            else
                error('C-centered monoclinic bandstructure plot requires to set the optional argument POSCAR!')
            end
        case {14,'tri','triclinic'}
            symmetry=14;% Triclinic
            TRI1A=false;
            TRI1B=false;
            TRI2A=false;
            TRI2B=false;
            if ~isnan(POSCAR)
                lattice=POSCAR.vec*POSCAR.acell;
                if all(lattice([3 7:8])==0)&&lattice(1,1)==-lattice(2,1)&&lattice(1,2)==lattice(2,2)
                    alpha=pi-acos((norm(lattice(2,:))^2+norm(lattice(3,:))^2-norm(lattice(2,:)+lattice(3,:))^2)/(2*norm(lattice(2,:))*norm(lattice(3,:))));
                    beta=pi-acos((norm(lattice(3,:))^2+norm(lattice(1,:))^2-norm(lattice(3,:)+lattice(1,:))^2)/(2*norm(lattice(3,:))*norm(lattice(1,:))));
                    gamma=pi-acos((norm(lattice(1,:))^2+norm(lattice(2,:))^2-norm(lattice(1,:)+lattice(2,:))^2)/(2*norm(lattice(1,:))*norm(lattice(2,:))));
                    a=lattice(1,1);
                    b=lattice(2,2)/cos(gamma);
                    c=lattice(3,3)/cos(beta);
                    a1=[a 0 0];
                    a2=[b*cos(gamma) b*sin(gamma) 0];
                    a3=[c*cos(beta) c*(cos(alpha)-cos(beta)*cos(gamma))/(sin(gamma))];
                    if ~all((lattice(:)-[a1;a2;a3])==0)
                        error('Unknown lattice vectors. Check POSCAR.')
                    end
                else
                    error('Unknown lattice vectors. Check POSCAR.')
                end
                V=dot(a1,cross(a2,a3));
                b1=2*pi*cross(a2,a3)/V;
                b2=2*pi*cross(a3,a1)/V;
                b3=2*pi*cross(a1,a2)/V;
                kalp=pi-acos((norm(b2)^2+norm(b3)^2-norm(b2+b3)^2)/(2*norm(b2)*norm(b3)));
                kbet=pi-acos((norm(b1)^2+norm(b3)^2-norm(b1+b3)^2)/(2*norm(b1)*norm(b3)));
                kgam=pi-acos((norm(b1)^2+norm(b2)^2-norm(b1+b2)^2)/(2*norm(b1)*norm(b2)));
                if kgam>pi/2&&kbet>pi/2&&kalp>pi/2
                    TRI1A=true;
                elseif kgam<pi/2&&kbet<pi/2&&kalp<pi/2
                    TRI1B=true;
                elseif kgam==pi/2&&kbet>pi/2&&kalp>pi/2
                    TRI2A=true;
                elseif kgam==pi/2&&kbet<pi/2&&kalp<pi/2
                    TRI2B=true;
                end
                if TRI1A||TRI2A
                    %! Gamma-point %! L-point %! M-point %! N-point %! R-point %! X-point %! Y-point %! Z-point
                    HS=[0 0 0;1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2; 1/2 1/2 1/2; 1/2 0 0; 0 1/2 0; 0 0 1/2];
                    HS_name={'\Gamma','L','M','N','R','X','Y','Z'};
                elseif TRI1B||TRI2B
                    %! Gamma-point %! L-point %! M-point %! N-point %! R-point %! X-point %! Y-point %! Z-point
                    HS=[0 0 0; 1/2 -1/2 0;  0 0 1/2;  -1/2 -1/2 1/2; 0 -1/2 1/2; 0 -1/2 0; 1/2 0 0; -1/2 0 1/2];
                    HS_name={'\Gamma','L','M','N','R','X','Y','Z'};
                end
            else
                error('C-centered monoclinic bandstructure plot requires to set the optional argument POSCAR!')
            end
        % case {15,'custom'}
        %     % Custom path
        %     fid=fopen(custom_filename);
        %     %! Gamma-point %! L-point %! M-point
        %     HS=[0 0 0;1/3 1/3 0;1/2 0 0];
        %     HS_name={'\Gamma','K','M'};
        %     fclose(fid);
        otherwise
            issymmetry=false;
            disp('You did not specify the symmetry points, or inserted something unrecognized')
            disp('Please choose one of the option below')
            disp('1) CUB/Cubic')
            disp('2) FCC/Face-centered cubic')
            disp('3) BCC/Body-centered cubic')
            disp('4) TET/Tetragonal')
            disp('5) BCT/Body-centered tetragonal')
            disp('6) ORC/Orthorhombic')
            disp('7) ORCF/Face-centered orthorhombic')
            disp('8) ORCI/Body-centered orthorhombic')
            disp('9) ORCC/C-centered orthorhombic')
            disp('10) HEX/Hexagonal')
            disp('11) RHL/Rhombohedral')
            disp('12) MCL/Monoclinic')
            disp('13) MCLC/C-centered monoclinic')
            disp('14) TRI/Triclinic')
            disp('15) Custom (not yet fully functional)')
            disp('16) None (not yet fully functional)')
            symmetry=input('Insert corresponding number ==> ');
    end
end

k=EIGENVAL.k;
EFermi=DOSCAR.EFermi;
if EIGENVAL.Ispin
    E_u=EIGENVAL.E_u;
    Occ_u=EIGENVAL.Occ_u;
    E_d=EIGENVAL.E_d;
    Occ_d=EIGENVAL.Occ_d;
    nkpts=size(E_u,1);
    nbands=size(E_u,2);
else
    E=EIGENVAL.E;
    Occ=EIGENVAL.Occ;
    nkpts=size(E,1);
    nbands=size(E,2);
end

idx_rm=[];
HS_str={};
HS_idx=[];
is_HS=zeros(nkpts,1);
is_HS([1 end])=true;
angle=zeros(nkpts-1,1);
angle([1 end])=0;
for p=2:nkpts-2
    d1=k(p,:)-k(p-1,:);
    d2=k(p+2,:)-k(p+1,:);
    angle(p)=atan2(norm(cross(d1,d2)),dot(d1,d2))*180/pi;
end
angle=angle.*(angle>1);
idx_chdir=find(angle>0);
for p=1:length(idx_chdir)
    if (angle(idx_chdir(p)-1)==0&&angle(idx_chdir(p)+1)==0)||(angle(idx_chdir(p)-1)>0&&angle(idx_chdir(p)+1)>0)
        is_HS(idx_chdir(p))=true;
        is_HS(idx_chdir(p)+1)=true;
    end
end
idxHS=find(is_HS);
for p=1:length(idxHS)
    HSdef=false;
    for q=1:length(HS(:,1))
        if norm(k(idxHS(p),:)-HS(q,:))<1e-3
            HS_str=[HS_str,HS_name{q}];
            HSdef=true;
        end
    end
    if ~HSdef
        HS_str=[HS_str,'?'];
        warning(['k-point #',num2str(idxHS(p)),' is a high symmetry point, but not specified in the symmetry you selected.'])
    end
end

HS_discontinuity=zeros(length(HS_str)/2+1);
HS_label{1}=HS_str{1};
for p=2:2:length(HS_str)-1
    if strcmpi(HS_str{p+1},HS_str{p})
        HS_label{p/2+1}=HS_str{p};
    else
        HS_label{p/2+1}=[HS_str{p},'/',HS_str{p+1}];
        HS_discontinuity(p/2+1)=true;
    end
end
HS_label{end+1}=HS_str{end};

BZ_int=reshape(find(is_HS),2,length(find(is_HS))/2)';
n_int=length(BZ_int(:,1));
for p=1:n_int
    for q=BZ_int(p,1):BZ_int(p,2)
        dk_int{p}(q-BZ_int(p,1)+1)=norm(k(q,:)-k(BZ_int(p,1),:));
    end
    deltak(p)=dk_int{p}(end);
    dk_int{p}=dk_int{p};
end
kx=[0 dk_int{1}(2:end)];
for p=2:n_int
    kx=[kx sum(deltak(1:p-1))+dk_int{p}(2:end)];
end
idx_rm=find(is_HS(2:end).*is_HS(1:end-1));
k(idx_rm,:)=[];
is_HS(idx_rm)=[];
if EIGENVAL.Ispin
    E_u(idx_rm,:)=[];
    E_d(idx_rm,:)=[];
    Occ_u(idx_rm,:)=[];
    Occ_d(idx_rm,:)=[];
    plot(kx,E_u-EFermi,'r')
    hold on
    plot(kx,E_d-EFermi,'b')
else
    E(idx_rm,:)=[];
    Occ(idx_rm,:)=[];
    plot(kx,E-EFermi,'k')
    hold on
end

xlim([min(kx) max(kx)])
yl=ylim;
kx_HS=kx(find(is_HS));
for p=2:length(kx_HS)-1
    if HS_discontinuity(p)
        h=plot([kx_HS(p) kx_HS(p)],[yl(1) yl(2)],'--','color',[1 1 1]*0.75);
        uistack(h,'bottom');
    else
        h=plot([kx_HS(p) kx_HS(p)],[yl(1) yl(2)],'color',[1 1 1]*0.75);
        uistack(h,'bottom');
    end
end
xticks(kx_HS)
xticklabels(HS_label);
set(gca,'fontsize',12,'fontname','cambria math')

E=DOSCAR.E;
DOS=DOSCAR.DOS;
lmDOS=DOSCAR.PDOS;
Ispin=DOSCAR.spin;
EFermi=DOSCAR.EFermi;
fig=figure('Units','pixels','Position',[200 120 647.2 425]); % [X X 560 425] by default
ax=axes('Parent',fig);
hold(ax,'on');
box on
set(gca,'Units','normalized')
set(gca,'Position',[0.1125 0.1087 0.6708 0.8054])
xlim([min(E) max(E)])
% Plot curves and store handles
if ~Ispin
    displayname='Total DOS';
    DOS_tot=plot(ax,E-EFermi,DOS(:,1),'LineWidth',1,'DisplayName',displayname);
    uicontrol(fig,'Style','checkbox','String',displayname,'Value',1,'Units','pixels','Position',[520 364.5 100 25],...
        'Callback',@(src,~)toggleVisibility(src,DOS_tot));
else
    displayname='Total DOS up';
    DOS_tot_u=plot(ax,E-EFermi,DOS(:,1),'LineWidth',1,'DisplayName',displayname);
    uicontrol(fig,'Style', 'checkbox','String',displayname','Value', 1,'Units','pixels','Position',[520 364.5 75 25],...
        'Callback',@(src,~)toggleVisibility(src,DOS_tot_u));
    displayname='Total DOS down';
    DOS_tot_d=plot(ax,E-EFermi,-DOS(:,2),'LineWidth',1,'DisplayName',displayname);
    uicontrol(fig,'Style','checkbox','String',displayname,'Value',1,'Units','pixels','Position',[520 339.5-25 75 25],...
        'Callback',@(src,~)toggleVisibility(src,DOS_tot_d));
end
if ~isnan(lmDOS)
    lmDOS_label=DOSCAR.lm_labels;
    for p=1:length(lmDOS(1,:))
        if Ispin
            lmDOS_plt(p)=plot(ax,E-EFermi,((-1)^p-1)*lmDOS(:,p),'LineWidth',1,'DisplayName',lmDOS_label{p},'visible','off');
        else
            lmDOS_plt(p)=plot(ax,E-EFermi,lmDOS(:,p),'LineWidth',1,'DisplayName',lmDOS_label{p},'visible','off');
        end
        uicontrol(fig,'Style','checkbox','String',lmDOS_label{p},'Value',0,'Units','pixels',...
            'Position',[520 364.5-(p+Ispin)*25 75 25],'Callback',@(src,~)toggleVisibility(src,lmDOS_plt(p)));
    end
    grid on
end
end