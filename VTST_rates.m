function R_out=VTST_rates(path,press,V,T,molecules,concentrations,gas_Fwd,gas_Bwd,varargin)
%==================================================================================================================================%
% VTST_rates.m:  transition state theory calculation of reaction rate coefficient (v0.5.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation (based on former bits and pieces from my thesis)
%       author: EYG
%   version 0.2 (25/08/2025) - Modification of the OUTPUT structure which the function returns
%       author: EYG
%   version 0.2.1 (28/08/2025) - Correction to the code to accomodate the modifications in 'readEnergy.m'
%       contrib: EYG
%   version 0.3 (12/09/2025) - Modification of the code to include the new formatting of the EnergyPathway structure type. This
%       author: EYG             removes the POSCAR and energy reading, along with the reaction coordinates computation, which is 
%                               now performed inside NEB_analysis.m
%   version 0.4 (02/03/2026) - Change of the input to specify the gas species via a string array and an array specifying the 
%       author: EYG             concentrations. Both normal and variational calculations can be handled with this program. Multiple
%                               sanity checks are performed: on the number of degree of freedom, on the fact that all DOF are free,
%                               on the difference between the final structure of the NEB calculations and the inputs of the
%                               vibrational calculations, etc...
%   version 0.5 (06/09/2026) - Add an optional argument to remove images from the calculation. In addition, information related to
%       author: EYG             the vibrational spectrum are now contained in the VIB structure, saved in a ".mat" file. The 
%                               EnergyPathway structure is no longer generated on-the-fly but read from a ".mat" file stored in
%                               path/images. ZPE is now computed for each image and the endpoints within the T loop.
%   version 0.5.1 (12/10/2026) - The warning to the user about the excess imaginary frequency now discriminate between one and more
%       contrib: EYG                than one extra imaginary frequency
%==================================================================================================================================%
% args:
%   path:           Location of reaction directory
%   press:          Pressure in Pascals
%   V:              Volume in cubic meter
%   T:              Temperature in Kelvin (array)
%   molecules:      String array of the formula of the gas species
%                   -> Only H, H2, CH3 and CH4 have been studied so far, carefully check results for other molecules
%   concentrations: Array of concentration of the corresponding chemicals in molecules
%   gas_Fwd:        index of the molecule to be considered as part of the reactants (default: NaN)
%   gas_Bwd:        index of the molecule to be considered as part of the products (default: NaN)
% opt. args:            'rm_images', followed by the indices of the images to be removed from the calculation.
%                           (default: NaN)
%                       'fh_crit', 
%==================================================================================================================================%
load('constant_fund.mat','Na','kB','e','h')
rm_image_idx=NaN;
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'rm_image'
                rm_image_idx=varargin{p+1};
        end
    end
end

if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end
E_tol=1e-3; % Tolerance for the energy comparisons
ds_tol=1e-3; % Tolerance for the POSCAR/CONTCAR comparisons
disp(['a tolerance of ',num2str(E_tol*1e3),' meV is applied to the energy comparisons.'])
disp(['a tolerance of ',num2str(ds_tol),' Angstroms is applied to the POSCAR/CONTCAR comparisons.'])

% Auto-check if there is a need for a variational calculation.
TS_dir=dir([path,'*.mat']);
if length(TS_dir)~=1
    error('There are either zero or more than one VIB structure in the directory!')
end
load([path,TS_dir.name]);
n_lTS=length(VIB);

% Read images to extract the geometry and energy of the system along the minimum energy pathway
if ~exist([path,'Images/NEB.mat'])
    error('No NEB.mat file in the image directory')
else
    load([path,'Images/NEB.mat'])
end

% Reading of the POSCARs
for p=1:n_lTS
    if n_lTS>1
        ds=ds_POSCAR(NEB.CONTCAR(p),VIB(p).POSCAR);
        if ds>ds_tol
            warning(sprintf(['The CONTCAR of the NEB do not match with the POSCAR of vibrational calculation for image ', ...
                VIB(p).POSCAR.Title,'.\nThere is a difference of ',num2str(ds,'%6.2e'),' Angstrom.']))
        end
    else
        for q=1:length(NEB.CONTCAR)
            ds(q)=ds_POSCAR(NEB.CONTCAR(q),VIB.POSCAR);
        end
        [minds,mindsidx]=min(ds);
        if minds>ds_tol
            warning(sprintf(['The CONTCARs of the NEB do not match with the POSCAR of vibrational calculation.' ...
                '\nClosest match is image ',NEB.POSCAR(mindsidx-1).Title,' with a difference of ',num2str(minds,'%6.2e'),' Angstrom.']))
        end
    end
    nDOF(p)=3*sum(VIB(p).POSCAR.n_chemicals);
end
if any(nDOF~=nDOF(1))
    error('Severe problem found: not all calculations were computed with the same number of atoms!')
end
nDOF=nDOF(1);

% Reading of the initial energy of the vibrational calculation and sanity check.
for p=1:n_lTS
    if n_lTS>1
        dE=abs(VIB(p).E-NEB.energies(p,end));
        if dE>E_tol
            warning(sprintf(['There is a mismatch (of ',num2str(dE,'%8.3E'),' eV) between the energy values of the NEB ' ...
                'calculation\n(for image ',num2str(p),') and that of the corresponding vibrational calculation']))
        end
        TS.E_Fwd(p)=NEB.energies(p,end)-NEB.energies(1,end);
        TS.E_Bwd(p)=NEB.energies(p,end)-NEB.energies(end,end);
    else
        [minE,minEidx]=min(abs(VIB.E-NEB.energies(:,end)));
        if minE>E_tol
            warning(sprintf(['There is a mismatch (of ',num2str(minE,'%8.3E'),' eV) between the energy values of the NEB ' ...
                'calculation\n(for image ',num2str(minEidx-1),') and that of the corresponding vibrational calculation']))
        end
        TS.E_Fwd=max(NEB.energies(:,end)-NEB.energies(1,end));
        TS.E_Bwd=max(NEB.energies(:,end)-NEB.energies(end,end));
    end
end

% Reading of the vibrational frequencies and calculation of the zero-point energy of the slabs.
for p=1:n_lTS
    if sum(VIB(p).nu==0)~=3
        error('Severe problem found: there are fewer or more than three zero-frequency vibrational mode!')
    end
end

% Only show once if an image is not taken into account due to the absence of an imaginary frequency
dispwarning=false(1,n_lTS);

% Extraction of the geometries of endpoints and TS
s0_dir=dir([path,'Reactants/*.mat']);
Data=load([path,'Reactants/',s0_dir.name],'VIB');
s0=Data.VIB;
sf_dir=dir([path,'Products/*.mat']);
Data=load([path,'Products/',sf_dir.name]);
sf=Data.VIB;
if ~all(s0.POSCAR.constraint)
    error('Not all degrees of freedom were computed for the reactants.')
end
if ~all(sf.POSCAR.constraint)
    error('Not all degrees of freedom were computed for the products.')
end

% Extraction of the raw frequencies provided by HIVE
if sum(s0.nu==0)~=3||sum(sf.nu==0)~=3
    error('Severe error found: there seems to be less than 3 zero frequency vibrational modes for the reactants or products.')
elseif any(s0.nu<0)||any(sf.nu<0)
    error('Severe error found: there is an imaginary frequency in the reactant or products')
end

% Past this point everything is T-dependent (FOR LOOP)
for p=1:length(T)
    % Calculation of the gas-related properties (total partition function, concentration, ...) given the input parameters
    gas=gas_modelling(press,V,T(p),molecules,concentrations);
    % Calculation of the vibrational partition function of the slabs
    [~,nu_vib_I(p),nu_rot_free_I(p)]=VibRot_Freq(s0,T(p),'mode',1);
    [~,~,nu_rot_hind_I(p)]=VibRot_Freq(s0,T(p),'mode',2);
    Qv_vib_I=Qv(nu_vib_I(p),T(p));
    Qv_rot_free_I=Qv(nu_rot_free_I(p),T(p));
    Qv_rot_hind_I=Qv(nu_rot_hind_I(p),T(p));
    if Qv_rot_free_I<Qv_rot_hind_I
        Rotor_mode=1;
    else
        Rotor_mode=2;
    end
    [nu_I,nu_vib_I(p),nu_rot_I(p)]=VibRot_Freq(s0,T(p),'mode',Rotor_mode);
    [~,nu_vib_F(p),nu_rot_free_F(p)]=VibRot_Freq(sf,T(p),'mode',1);
    [~,~,nu_rot_hind_F(p)]=VibRot_Freq(sf,T(p),'mode',2);
    Qv_vib_F=Qv(nu_vib_F(p),T(p));
    Qv_rot_free_F=Qv(nu_rot_free_F(p),T(p));
    Qv_rot_hind_F=Qv(nu_rot_hind_F(p),T(p));
    if Qv_rot_free_F<Qv_rot_hind_F
        Rotor_mode=1;
    else
        Rotor_mode=2;
    end
    [nu_F,nu_vib_F(p),nu_rot_F(p)]=VibRot_Freq(sf,T(p),'mode',Rotor_mode);
    % Removal of the zero-frequency modes
    nu_I=nu_I(nu_I~=0);
    nu_F=nu_F(nu_F~=0);
    % Computation of the zero-point energy
    ZPE_I=ZPE(nu_I);
    ZPE_F=ZPE(nu_F);
    % Computation of the vibrational partition functions of reactants and products
    Qv_I=Qv(nu_I,T(p));
    Qv_F=Qv(nu_F,T(p));
    for q=1:n_lTS
        if any(rm_image_idx==q)
            TS.nu(q,:)=NaN(1,nDOF-4);
            TS.nu_d(q)=NaN;
            kappa_Fwd_var(p,q)=NaN;
            kappa_Bwd_var(p,q)=NaN;
            A_Fwd_var(p,q)=NaN;
            A_Bwd_var(p,q)=NaN;
            k_Fwd_var(p,q)=NaN;
            k_Bwd_var(p,q)=NaN;
        else
            [~,~,nu_rot_free(p,q)]=VibRot_Freq(VIB(q),T(p),'mode',1);
            [~,~,nu_rot_hind(p,q)]=VibRot_Freq(VIB(q),T(p),'mode',2);
            if Qv(nu_rot_free(p,q),T(p))<Qv(nu_rot_hind(p,q),T(p))
                Rotor_mode=1;
            else
                Rotor_mode=2;
            end
            [nu{q},nu_vib(p,q),nu_rot(p,q)]=VibRot_Freq(VIB(q),T(p));
            if sum(nu{q}==0)~=3
                error('Severe problem found: there are fewer or more than three zero-frequency vibrational mode!')
            elseif sum(nu{q}<0)<1
                if ~dispwarning(q)
                    dispwarning(q)=true;
                    warning(['Vibration calculation of image ',num2str(q),' (',VIB(q).POSCAR.Title,') has no imaginary frequency, and will be ignored.'])
                end
                TS.nu(q,:)=NaN(1,nDOF-4);
                TS.nu_d(q)=NaN;
            elseif sum(nu{q}<0)>1
                if ~dispwarning(q)
                    dispwarning(q)=true;
                    if sum(nu{q}<0)==2
                        warning(['Vibration calculation of image ',num2str(q),' (',VIB(q).POSCAR.Title,') has one extra imaginary frequency, and will be ignored.'])
                    else
                        warning(['Vibration calculation of image ',num2str(q),' (',VIB(q).POSCAR.Title,') has more than one extra imaginary frequency, and will be ignored.'])
                    end
                end
                TS.nu(q,:)=NaN(1,nDOF-4);
                TS.nu_d(q)=NaN;
            else
                TS.nu(q,:)=nu{q}(nu{q}>0);
                TS.nu_d(q)=-nu{q}(nu{q}<0);
            end
            TS.ZPE(q)=ZPE([TS.nu(q,:)]);
            TS.delta_ZPE_Fwd(q)=(TS.ZPE(q))-(ZPE_I);
            TS.delta_ZPE_Bwd(q)=(TS.ZPE(q))-(ZPE_F);
            % Computation of the Skodje-Thrular coefficient to account for quantum tunneling through the energy barrier
            if NEB.energies(q,end)>NEB.energies(1,end)&&NEB.energies(q,end)>NEB.energies(end,end)
                kappa_Fwd_var(p,q)=SkodjeTruhlar(T(p),TS.nu_d(q),TS.E_Fwd(q)*e);
                kappa_Bwd_var(p,q)=SkodjeTruhlar(T(p),TS.nu_d(q),TS.E_Bwd(q)*e);
            else
                kappa_Fwd_var(p,q)=1;
                kappa_Bwd_var(p,q)=1;
            end
            TS.Qv(q,:)=Qv(TS.nu(q,:),T(p));
            % Each type of gas-phase reactant or product has a specific effect on the calculation. In addition, when no gas-phase 
            % species is involved in either direction, there is no N_A in the prefactor A (order 0 reaction).
            A_gas_Fwd=1;
            A_gas_Bwd=1;
            if ~isnan(gas_Fwd)
                A_gas_Fwd=Na*exp(gas(gas_Fwd).ZPE/(kB*T(p)))/gas(gas_Fwd).Q_tot;
            end
            if ~isnan(gas_Bwd)
                A_gas_Bwd=Na*exp(gas(gas_Bwd).ZPE/(kB*T(p)))/gas(gas_Bwd).Q_tot;
            end
            % Eyring-Polanyi-Evans equation
            A_Fwd_var(p,q)=A_gas_Fwd*kappa_Fwd_var(p,q)*(kB*T(p)/h)*Qv_div(TS.Qv(q,:),Qv_I)*exp(-TS.delta_ZPE_Fwd(q)/(kB*T(p)));
            A_Bwd_var(p,q)=A_gas_Bwd*kappa_Bwd_var(p,q)*(kB*T(p)/h)*Qv_div(TS.Qv(q,:),Qv_F)*exp(-TS.delta_ZPE_Bwd(q)/(kB*T(p)));
            % The reaction rate coefficient (at the temperature T) is given by the product of the prefactor and the exponential of the 
            % energy barrier
            k_Fwd_var(p,q)=A_Fwd_var(p,q)*exp(-max(VIB(q).E-NEB.energies(1,end))*e/(kB*T(p)));
            k_Bwd_var(p,q)=A_Bwd_var(p,q)*exp(-max(VIB(q).E-NEB.energies(end,end))*e/(kB*T(p)));
        end
    end
    [k_Fwd(p),idx_Fwd(p)]=min(k_Fwd_var(p,:));
    [k_Bwd(p),idx_Bwd(p)]=min(k_Bwd_var(p,:));
    if idx_Fwd~=idx_Bwd
        disp('different indices between forward and backward')
    end
    kappa_Fwd(p)=kappa_Fwd_var(p,idx_Fwd(p));
    kappa_Bwd(p)=kappa_Bwd_var(p,idx_Bwd(p));
    A_Fwd(p)=A_Fwd_var(p,idx_Fwd(p));
    A_Bwd(p)=A_Bwd_var(p,idx_Bwd(p));
    E_Fwd(p)=TS.E_Fwd(idx_Fwd(p));
    E_Bwd(p)=TS.E_Bwd(idx_Bwd(p));
    % The effective reaction rate coefficient is obtained by multiplicating the concentration of the gas phase species by the
    % reaction rate coefficient
    r_Fwd(p)=k_Fwd(p);
    r_Bwd(p)=k_Bwd(p);
    n_Fwd(p)=1;
    n_Bwd(p)=1;
    if ~isnan(gas_Fwd)
        r_Fwd(p)=k_Fwd(p)*gas(gas_Fwd).n;
        n_Fwd(p)=gas(gas_Fwd).n;
    end
    if ~isnan(gas_Bwd)
        r_Bwd(p)=k_Bwd(p)*gas(gas_Bwd).n;
        n_Bwd(p)=gas(gas_Bwd).n;
    end
end

% The function outputs a structure, containing almost every variable in an orderly fashion
% First, all the input parameters are stored in the input field
slashtrim=find(path(1:end-1)=='/');
if isempty(slashtrim)
    R_out.Input.Reaction_ID=path(1:end-1);
else
    R_out.Input.Reaction_ID=path(slashtrim(end)+1:end-1);
end
R_out.Input.T=T;
R_out.Input.P=press;
R_out.Input.V=V;
R_out.Input.gas_Fwd=gas_Fwd;
R_out.Input.gas_Bwd=gas_Bwd;
R_out.Input.molecules=molecules;
R_out.Input.concentrations=concentrations;

% Second, the gas properties are stored using the already-existing structure format for each gas
R_out.gas=gas;

% Third, the "static" properties of reactants, products, images and TS are saved. Each subfield of TST is itself a structure, see
% above for their precise definition
R_out.TST.Reactants=s0;
R_out.TST.Products=sf;
R_out.TST.Products=sf;
R_out.TST.nu_vib_I=nu_vib_I;
R_out.TST.nu_rot_I=nu_rot_I;
R_out.TST.nu_rot_F=nu_rot_F;
R_out.TST.nu_vib_F=nu_vib_F;
R_out.TST.NEB=NEB;
R_out.TST.TS=TS;
R_out.TST.idx_Fwd=idx_Fwd;
R_out.TST.idx_Bwd=idx_Bwd;
R_out.TST.VIB=VIB;
R_out.TST.nu_rot=nu_rot;
R_out.TST.nu_rot_hind=nu_rot_hind;
R_out.TST.nu_rot_free=nu_rot_free;
R_out.TST.nu_vib=nu_vib;

% Finally, the reaction rate coefficients and related parameters are stored. The energies are expressed in eV, the effective
% reaction rate coefficient in s^-1, and the units of A and k depends on the presence or absence of a gas phase particle.
R_out.rates.Fwd.dE=E_Fwd;
R_out.rates.Fwd.kappa=kappa_Fwd;
R_out.rates.Fwd.A=A_Fwd;
R_out.rates.Fwd.A_var=A_Fwd_var;
R_out.rates.Fwd.k_var=k_Fwd_var;
R_out.rates.Fwd.k=k_Fwd;
R_out.rates.Fwd.n=n_Fwd;
R_out.rates.Fwd.r=r_Fwd;
R_out.rates.Bwd.dE=E_Bwd;
R_out.rates.Bwd.kappa=kappa_Bwd;
R_out.rates.Bwd.A=A_Bwd;
R_out.rates.Bwd.A_var=A_Bwd_var;
R_out.rates.Bwd.k_var=k_Bwd_var;
R_out.rates.Bwd.k=k_Bwd;
R_out.rates.Bwd.n=n_Bwd;
R_out.rates.Bwd.r=r_Bwd;
% The fit is saved as well, but only if the temperature array is large enough to perform the fit: if there are more than 5 
% temperature value, a fit to the equation A0*T^n*exp(Ea/(kB*T)) for both the forward and backward reaction (works best with a 
% significant amount of datapoints, minimum is usually around 5 but I use ~ 2000 to be accurate)
if length(T)>5
    [A0_Fwd,nT_Fwd,Ea_Fwd]=fit_arrh(T,k_Fwd);
    [A0_Bwd,nT_Bwd,Ea_Bwd]=fit_arrh(T,k_Bwd);
    k_Fwd_fit=A0_Fwd*T.^nT_Fwd.*exp(-Ea_Fwd./(kB*T));
    k_Bwd_fit=A0_Bwd*T.^nT_Bwd.*exp(-Ea_Bwd./(kB*T));
    disp(size(k_Fwd_fit-k_Fwd))
    diff_Fwd=100*abs(k_Fwd_fit-k_Fwd)./k_Fwd;
    diff_Bwd=100*abs(k_Bwd_fit-k_Bwd)./k_Bwd;
    mu_Fwd=mean(diff_Fwd);
    mu_Bwd=mean(diff_Bwd);
    sigma_Fwd=std(diff_Fwd);
    sigma_Bwd=std(diff_Bwd);
    R_out.rates.fit.Fwd.A0=A0_Fwd;
    R_out.rates.fit.Fwd.n=nT_Fwd;
    R_out.rates.fit.Fwd.Ea=Ea_Fwd/e;
    R_out.rates.fit.Fwd.mu=mu_Fwd;
    R_out.rates.fit.Fwd.sigma=sigma_Fwd;
    R_out.rates.fit.Bwd.A0=A0_Bwd;
    R_out.rates.fit.Bwd.n=nT_Bwd;
    R_out.rates.fit.Bwd.Ea=Ea_Bwd/e;
    R_out.rates.fit.Bwd.mu=mu_Bwd;
    R_out.rates.fit.Bwd.sigma=sigma_Bwd;
else
    R_out.rates.fit='not applicable';
end

% clean-up of all the variable but the OUTPUT structure, then this structure is save in a .mat file for later use
% clearvars -except R_out
save(['Reaction_',R_out.Input.Reaction_ID,'.mat'],'R_out')
