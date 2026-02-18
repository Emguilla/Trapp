function R_out=VTST_rates(path,press,V,T,molecules,concentrations,gas_Fwd,gas_Bwd)
%==================================================================================================================================%
% VTST_rates.m:  transition state theory calculation of reaction rate coefficient (v0.3)
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
%   version 0.4 (19/02/2026) - Change of the input to specify the gas species via a string array and an array specifying the 
%                               concentrations. Both normal and variational calculations can be handled with this program. Multiple
%                               sanity checks are performed: on the number of degree of freedom, on the fact that all DOF are free,
%                               on the difference between the final structure of the NEB calculations and the inputs of the
%                               vibrational calculations, etc...
% Upcoming modifs: allow for difference between the NEB images and the actual images which phonon spectrum was computed. On way to
% do so could be by checking which NEB image is closer to the phonon image (with a tolerance of 0.01 Angstrom maybe?). In case there
% are important image missing (e.g. endpoints), use the energy of the NEB calculation instead of the first energy of the phonon
% calculation?
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
%==================================================================================================================================%
% Default parameters
press=25e3;
V=1e-6;
T=1200;
molecules={'H','H2','CH3','CH4'};
concentrations=[0.01 1 0 0];
path='A0A0+CH3_A1A0+CH4';
gas_Fwd=3;
gas_Bwd=4;
path='A0A0+H_A1A0+H2';
gas_Fwd=1;
gas_Bwd=2;
path='A0A0_A1A0+H';
gas_Fwd=NaN;
gas_Bwd=1;
%==================================================================================================================================%
load('constant_fund.mat','Na','kB','e','h')
if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end

E_tol=1e-3; % Tolerance for the energy comparisons
ds_tol=1e-3; % Tolerance for the POSCAR/CONTCAR comparisons
disp(['a tolerance of ',num2str(E_tol*1e3),' meV is applied to the energy comparisons.'])
disp(['a tolerance of ',num2str(ds_tol),' Angstroms is applied to the POSCAR/CONTCAR comparisons.'])

% Auto-check if there is a need for a variational calculation.
if exist([path,'lTS'],'dir')
    ldir=dir([path,'lTS/']);
    ldir=ldir(3:end);
    variational='lTS/';
    n_lTS=length(ldir);
elseif exist([path,'tTS'],'dir')
    ldir=dir([path,'tTS*']);
    variational='';
    n_lTS=1;
end

% Read images to extract the geometry and energy of the system along the minimum energy pathway
NEB=NEB_analysis('path',[path,'/Images']);

% Extraction of the geometries of endpoints and TS
s0.POSCAR=readPOSCAR([path,'Reactants/POSCAR']);
if ~all(s0.POSCAR.constraint)
    error('Not all degrees of freedom were computed for the reactants.')
end
sf.POSCAR=readPOSCAR([path,'Products/POSCAR']);
if ~all(sf.POSCAR.constraint)
    error('Not all degrees of freedom were computed for the products.')
end

% Extraction of the raw frequencies provided by HIVE
s0.nu=readFrequencies([path,'Reactants/']);
sf.nu=readFrequencies([path,'Products/']);
if sum(s0.nu==0)~=3||sum(sf.nu==0)~=3
    error('Severe error found: there seems to be less than 3 zero frequency vibrational modes for the reactants or products.')
elseif any(s0.nu<0)||any(sf.nu<0)
    error('Severe error found: there is an imaginary frequency in the reactant or products')
end
% Removal of the zero-frequency modes
s0.nu=s0.nu(s0.nu~=0);
sf.nu=sf.nu(sf.nu~=0);
% Computation of the zero-point energy
s0.ZPE=ZPE(s0.nu);
sf.ZPE=ZPE(sf.nu);

% Reading of the POSCARs
for p=1:n_lTS
    TS.POSCAR(p)=readPOSCAR([path,variational,ldir(p).name,'/POSCAR']);
    if ~isempty(variational)
        ds=ds_POSCAR(NEB.CONTCAR(p),TS.POSCAR(p));
        if ds>ds_tol
            warning(sprintf(['The CONTCAR of the NEB do not match with the POSCAR of vibrational calculation in ', ...
                path,variational,ldir(p).name,'.\nThere is a difference of ',num2str(ds,'%6.2e'),' Angstrom.']))
        end
    else
        for q=1:length(NEB.CONTCAR)
            ds(q)=ds_POSCAR(NEB.CONTCAR(q),TS.POSCAR);
        end
        [minds,mindsidx]=min(ds);
        if minds>ds_tol
            warning(sprintf(['The CONTCARs of the NEB do not match with the POSCAR of vibrational calculation.' ...
                '\nClosest match is image ',num2str(mindsidx-1),' with a difference of ',num2str(minds,'%6.2e'),' Angstrom.']))
        end
    end
    nDOF(p)=3*sum(TS.POSCAR(p).n_chemicals);
end
if any(nDOF~=nDOF(1))
    error('Severe problem found: not all calculations were computed with the same number of atoms!')
end
nDOF=nDOF(1);

% Reading of the initial energy of the vibrational calculation and sanity check.
for p=1:n_lTS % Initial and final images have already been read.
    E=readEnergy([path,variational,ldir(p).name,'/'],'save',true);
    TS.E(p)=E(1);
    if ~isempty(variational)
        dE=abs(TS.E(p)-NEB.energies(p,end));
        if dE>E_tol
            warning(sprintf(['There is a mismatch (of ',num2str(dE,'%8.3E'),' eV) between the energy values of the NEB ' ...
                'calculation\n (for image ',num2str(p),') and that of the corresponding vibrational calculation']))
        end
        TS.E_Fwd(p)=NEB.energies(p,end)-NEB.energies(1,end);
        TS.E_Bwd(p)=NEB.energies(p,end)-NEB.energies(end,end);
    else
        [minE,minEidx]=min(abs(TS.E-NEB.energies(:,end)));
        if minE>E_tol
            warning(sprintf(['There is a mismatch (of ',num2str(minE,'%8.3E'),' eV) between the energy values of the NEB ' ...
                'calculation\n (for image ',num2str(minEidx-1),') and that of the corresponding vibrational calculation']))
        end
        TS.E_Fwd=max(NEB.energies(:,end)-NEB.energies(1,end));
        TS.E_Bwd=max(NEB.energies(:,end)-NEB.energies(end,end));
    end
end

% Reading of the vibrational frequencies and calculation of the zero-point energy of the slabs.
nonTS=[];
for p=1:n_lTS
    nu=readFrequencies([path,variational,ldir(p).name,'/']);
    if sum(nu==0)~=3
        error('Severe problem found: there are fewer or more than three zero-frequency vibrational mode!')
    elseif sum(nu<0)~=1
        warning(['Vibration calculation of image ',num2str(p),' has no imaginary frequency, and will be ignored.'])
        nonTS=[nonTS p];
        TS.nu(p,:)=NaN*zeros(1,nDOF-4);
        TS.nu_d(p)=NaN;
        TS.ZPE(p)=NaN;
        TS.delta_ZPE_Fwd(p)=NaN;
        TS.delta_ZPE_Bwd(p)=NaN;
    else
        TS.nu(p,:)=nu(nu>0);
        TS.nu_d=-nu(nu<0);
        TS.ZPE(p)=ZPE([TS.nu(p,:)]);
        TS.delta_ZPE_Fwd(p)=(TS.ZPE(p))-(s0.ZPE);
        TS.delta_ZPE_Bwd(p)=(TS.ZPE(p))-(sf.ZPE);
    end
end

% Past this point everything is T-dependent (FOR LOOP)
for p=1:length(T)
    % Calculation of the gas-related properties (total partition function, concentration, ...) given the input parameters
    gas=gas_modelling(press,V,T(p),molecules,concentrations);
    % Calculation of the vibrational partition function of the slabs
    s0.Qv=Qv(s0.nu,T(p));
    sf.Qv=Qv(sf.nu,T(p));
    for q=1:n_lTS
        % Computation of the Skodje-Thrular coefficient to account for quantum tunneling through the energy barrier
        if any((TS.E>NEB.energies(1,end)).*(TS.E>NEB.energies(end,end)))
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
        A_Fwd_var(p,q)=A_gas_Fwd*kappa_Fwd_var(p,q)*(kB*T(p)/h)*Qv_div(TS.Qv(q,:),s0.Qv)*exp(-TS.delta_ZPE_Fwd(q)/(kB*T(p)));
        A_Bwd_var(p,q)=A_gas_Bwd*kappa_Bwd_var(p,q)*(kB*T(p)/h)*Qv_div(TS.Qv(q,:),sf.Qv)*exp(-TS.delta_ZPE_Bwd(q)/(kB*T(p)));
        % The reaction rate coefficient (at the temperature T) is given by the product of the prefactor and the exponential of the 
        % energy barrier
        k_Fwd_var(p,q)=A_Fwd_var(p,q)*exp(-max(NEB.energies(:,end)-NEB.energies(1,end))*e/(kB*T(p)));
        k_Bwd_var(p,q)=A_Bwd_var(p,q)*exp(-max(NEB.energies(:,end)-NEB.energies(end,end))*e/(kB*T(p)));
    end
    [k_Fwd(p),idx_Fwd]=min(k_Fwd_var(p,:));
    [k_Bwd(p),idx_Bwd]=min(k_Bwd_var(p,:));
    if idx_Fwd~=idx_Bwd
        disp('different indices between forward and backward')
    end
    kappa_Fwd(p)=kappa_Fwd_var(p,idx_Fwd);
    kappa_Bwd(p)=kappa_Bwd_var(p,idx_Bwd);
    A_Fwd(p)=A_Fwd_var(p,idx_Fwd);
    A_Bwd(p)=A_Bwd_var(p,idx_Bwd);
    E_Fwd(p)=TS.E_Fwd(idx_Fwd);
    E_Bwd(p)=TS.E_Bwd(idx_Fwd);
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
disp(slashtrim)
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
R_out.TST.EnergyPathway=NEB;
R_out.TST.TS=TS;

% Finally, the reaction rate coefficients and related parameters are stored. The energies are expressed in eV, the effective
% reaction rate coefficient in s^-1, and the units of A and k depends on the presence or absence of a gas phase particle.
R_out.rates.Fwd.dE=E_Fwd;
R_out.rates.Fwd.kappa=kappa_Fwd;
R_out.rates.Fwd.A=A_Fwd;
R_out.rates.Fwd.k=k_Fwd;
R_out.rates.Fwd.n=n_Fwd;
R_out.rates.Fwd.r=r_Fwd;
R_out.rates.Bwd.dE=E_Bwd;
R_out.rates.Bwd.kappa=kappa_Bwd;
R_out.rates.Bwd.A=A_Bwd;
R_out.rates.Bwd.k=k_Bwd;
R_out.rates.Bwd.n=n_Bwd;
R_out.rates.Bwd.r=r_Bwd;
% The fit is saved as well, but only if the temperature array is large enough to perform the fit: if there are more than 5 
% temperature value, a fit to the equation A0*T^n*exp(Ea/(kB*T)) for both the forward and backward reaction (works best with a 
% significant amount of datapoints, minimum is usually around 5 but I use ~ 2000 to be accurate)
if length(T)>5
    [A0_Fwd,nT_Fwd,Ea_Fwd]=fit_arrh(T,k_Fwd);
    [A0_Bwd,nT_Bwd,Ea_Bwd]=fit_arrh(T,k_Bwd);
    R_out.rates.fit.Fwd.A0=A0_Fwd;
    R_out.rates.fit.Fwd.n=nT_Fwd;
    R_out.rates.fit.Fwd.Ea=Ea_Fwd/e;
    R_out.rates.fit.Bwd.A0=A0_Bwd;
    R_out.rates.fit.Bwd.n=nT_Bwd;
    R_out.rates.fit.Bwd.Ea=Ea_Bwd/e;
else
    R_out.rates.fit='not applicable';
end

% clean-up of all the variable but the OUTPUT structure, then this structure is save in a .mat file for later use
% clearvars -except R_out
save([R_out.Input.Reaction_ID,'/Reaction_',R_out.Input.Reaction_ID,'.mat'],'R_out')
