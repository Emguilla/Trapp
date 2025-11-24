function R_out=TST_rates(path,press,V,T,gas_Fwd,gas_Bwd,ratio_H_H2,ratio_CH4_H2,ratio_CH3_CH4)
%==================================================================================================================================%
% TST_rates.m:  transition state theory calculation of reaction rate coefficient (v0.3)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation (based on former bits and pieces from my thesis)
%       author: EYG
%   version 0.2 (25/08/2025) - Modification of the OUTPUT structure which the function returns
%       author: EYG
%   version 0.2.1 (28/08/2025) - Correction to the code to accomodate the modifications in 'readEnergy.m'
%       contrib: EYG
%   version 0.3 (12/09/2025) - Modification of the code to include the new formatting of the EnergyPathway structure type. This
%       author: EYG           removes the POSCAR and energy reading, along with the reaction coordinates computation, which is now
%                               performed inside NEB_analysis.m
%==================================================================================================================================%
% args:
%   path:                                       Location of reaction directory
%   press:                                      Pressure in Pascals
%   V:                                          Volume in cubic meter
%   T:                                          Temperature in Kelvin (array)
%   gas_Fwd, gas_Bwd:                           Radicals/molecules involved in the forward and backward direction of the reaction
%                                               -> Only H, H2, CH3 and CH4 are currently implemented
%   ratio_H_H2, ratio_CH4_H2, ratio_CH3_CH4:    Ratio between the different species.
%==================================================================================================================================%
% Default parameters
%     path='026_100_0A+H_to_1A+H2/';
%     T=1200;
%     V=1e-6;
%     press=25e3;
%     gas_Fwd='H';
%     gas_Bwd='H2';
%     ratio_H_H2=0.01;
%     ratio_CH4_H2=0.00;
%     ratio_CH3_CH4=0.1;
%==================================================================================================================================%
load('constant_fund.mat','Na','kB','e','h')
if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end
% Read images to extract the geometry and energy of the system along the minimum energy pathway
EnergyPathway=NEB_analysis('path',[path,'/Images']);
E=EnergyPathway.energies(:,end)*e;
E_Fwd=max(E-E(1));
E_Bwd=max(E-E(end));
% Extraction of the geometries of endpoints and TS
s0.POSCAR=readPOSCAR([path,'Endpoints/Reactants/POSCAR']);
sf.POSCAR=readPOSCAR([path,'Endpoints/Products/POSCAR']);
TS.POSCAR=readPOSCAR([path,'TS/POSCAR']);
% Extraction of the raw frequencies provided by HIVE
s0.nu=readFrequencies([path,'Endpoints/Reactants/']);
sf.nu=readFrequencies([path,'Endpoints/Products/']);
TS.nu=readFrequencies([path,'TS/']);
% Identification of the frequency of the vibrational mode at the saddle point (imaginary frequency)
TS.nu_d=-TS.nu(end);

% Removal of the 3-zeros frequencies related to the translational DOF of the slab
s0.nu=s0.nu(1:end-3);
sf.nu=sf.nu(1:end-3);
TS.nu=TS.nu(1:end-4); % Since we already stored nu_d, we can remove it as well

% If relevant, the last frequency in the list provided by HIVE is replaced with that of a hindered rotor
[s0.nu,hr_s0]=hr_freq(s0.nu,[path,'Endpoints/Reactants/']);
if hr_s0; s0.rot.mode='hindered';s0.rot.nu=s0.nu(end);end
[sf.nu,hr_sf]=hr_freq(sf.nu,[path,'Endpoints/Products/']);
if hr_sf; sf.rot.mode='hindered';sf.rot.nu=sf.nu(end);end
[TS.nu,hr_TS]=hr_freq(TS.nu,[path,'TS/']);
if hr_TS; TS.rot.mode='hindered';TS.rot.nu=TS.nu(end);end
% Past this point everything is T-dependent (FOR LOOP)
for p=1:length(T)
    % Calculation of the gas-related properties (total partition function, concentration, ...) given the input parameters
    [H(p),H2(p),CH3(p),CH4(p)]=gas_modelling(press,V,T(p),ratio_H_H2,ratio_CH4_H2,ratio_CH3_CH4);
    % If relevant, the last frequency in the list provided by HIVE is replaced with that of a free rotor
    [s0.nu,fr_s0]=fr_freq(s0.nu,[path,'Endpoints/Reactants/'],T(p));
    if fr_s0; s0.rot.mode='free';s0.rot.nu(p)=s0.nu(end);end
    [sf.nu,fr_sf]=fr_freq(sf.nu,[path,'Endpoints/Products/'],T(p));
    if fr_sf; sf.rot.mode='free';sf.rot.nu(p)=sf.nu(end);end
    [TS.nu,fr_TS]=fr_freq(TS.nu,[path,'TS/'],T(p));
    if fr_TS; TS.rot.mode='free';TS.rot.nu(p)=TS.nu(end);end
    % Calculation of the zero-point energy of the slabs. The ZPE of the molecules are included later
    s0.ZPE=ZPE(s0.nu);
    sf.ZPE=ZPE(sf.nu);
    TS.ZPE=ZPE(TS.nu);
    % Calculation of the vibrational partition function of the slabs
    s0.Qv=Qv(s0.nu,T(p));
    sf.Qv=Qv(sf.nu,T(p));
    TS.Qv=Qv(TS.nu,T(p));
    % Computation of the Skodje-Thrular coefficient to account for quantum tunneling through the energy barrier
    kappa_Fwd(p)=SkodjeTruhlar(T(p),TS.nu_d,E_Fwd);
    kappa_Bwd(p)=SkodjeTruhlar(T(p),TS.nu_d,E_Bwd);
    % Each type of gas-phase reactant or product has a specific effect on the calculation. In addition, when no gas-phase species is
    % involved in either direction, there is no Na in the prefactor A (order 0 reaction).
    switch gas_Fwd
        case 'H'
            delta_ZPE=(TS.ZPE)-(s0.ZPE);
            A_Fwd(p)=kappa_Fwd(p)*Na*(kB*T(p)/h)*Qv_div(TS.Qv,s0.Qv)/H(p).Q_tot*exp(-delta_ZPE/(kB*T(p)));
            n_Fwd(p)=H(p).n;
        case 'H2'
            delta_ZPE=(TS.ZPE)-(s0.ZPE+H2(p).ZPE);
            A_Fwd(p)=kappa_Fwd(p)*Na*(kB*T(p)/h)*Qv_div(TS.Qv,s0.Qv)/H2(p).Q_tot*exp(-delta_ZPE/(kB*T(p)));
            n_Fwd(p)=H2(p).n;
        case 'CH3'
            delta_ZPE=(TS.ZPE)-(s0.ZPE+CH3(p).ZPE);
            A_Fwd(p)=kappa_Fwd(p)*Na*(kB*T(p)/h)*Qv_div(TS.Qv,s0.Qv)/CH3(p).Q_tot*exp(-delta_ZPE/(kB*T(p)));
            n_Fwd(p)=CH3(p).n;
        case 'CH4'
            delta_ZPE=(TS.ZPE)-(s0.ZPE+CH4(p).ZPE);
            A_Fwd(p)=kappa_Fwd(p)*Na*(kB*T(p)/h)*Qv_div(TS.Qv,s0.Qv)/CH4(p).Q_tot*exp(-delta_ZPE/(kB*T(p)));
            n_Fwd(p)=CH4(p).n;
        case ''
            delta_ZPE=(TS.ZPE)-(s0.ZPE);
            A_Fwd(p)=kappa_Fwd(p)*(kB*T(p)/h)*Qv_div(TS.Qv,s0.Qv)*exp(-delta_ZPE/(kB*T(p)));
            n_Fwd(p)=1;
        otherwise
            error('This molecule has not been implemented yet')
    end
    % Same goes for the backward reaction
    switch gas_Bwd
        case 'H'
            delta_ZPE=(TS.ZPE)-(sf.ZPE);
            A_Bwd(p)=kappa_Bwd(p)*Na*(kB*T(p)/h)*Qv_div(TS.Qv,sf.Qv)/H(p).Q_tot*exp(-delta_ZPE/(kB*T(p)));
            n_Bwd(p)=H(p).n;
        case 'H2'
            delta_ZPE=(TS.ZPE)-(sf.ZPE+H2(p).ZPE);
            A_Bwd(p)=kappa_Bwd(p)*Na*(kB*T(p)/h)*Qv_div(TS.Qv,sf.Qv)/H2(p).Q_tot*exp(-delta_ZPE/(kB*T(p)));
            n_Bwd(p)=H2(p).n;
        case 'CH3'
            delta_ZPE=(TS.ZPE)-(sf.ZPE+CH3(p).ZPE);
            A_Bwd(p)=kappa_Bwd(p)*Na*(kB*T(p)/h)*Qv_div(TS.Qv,sf.Qv)/CH3(p).Q_tot*exp(-delta_ZPE/(kB*T(p)));
            n_Bwd(p)=CH3(p).n;
        case 'CH4'
            delta_ZPE=(TS.ZPE)-(sf.ZPE+CH4(p).ZPE);
            A_Bwd(p)=kappa_Bwd(p)*Na*(kB*T(p)/h)*Qv_div(TS.Qv,sf.Qv)/CH4(p).Q_tot*exp(-delta_ZPE/(kB*T(p)));
            n_Bwd(p)=CH4(p).n;
        case ''
            delta_ZPE=(TS.ZPE)-(sf.ZPE);
            A_Bwd(p)=kappa_Bwd(p)*(kB*T(p)/h)*Qv_div(TS.Qv,sf.Qv)*exp(-delta_ZPE/(kB*T(p)));
            n_Bwd(p)=1;
        otherwise
            error('This molecule has not been implemented yet')
    end
    % The reaction rate coefficient (at the temperature T) is given by the product of the prefactor and the exponential of the 
    % energy barrier
    k_Fwd(p)=A_Fwd(p)*exp(-E_Fwd/(kB*T(p)));
    k_Bwd(p)=A_Bwd(p)*exp(-E_Bwd/(kB*T(p)));
    % The effective reaction rate coefficient is obtained by multiplicating the concentration of the gas phase species by the
    % reaction rate coefficient
    r_Fwd(p)=k_Fwd(p)*n_Fwd(p);
    r_Bwd(p)=k_Bwd(p)*n_Bwd(p);
end
% If there are more than one temperature value, a fit to the equation A0*T^n*exp(Ea/(kB*T)) for both the forward and backward
% reaction (works best with a significant amount of datapoints, minimum is usually around 5 but I use ~ 2000 to be accurate)
if length(T)>5
    [A0_Fwd,nT_Fwd,Ea_Fwd]=fit_arrh(T,k_Fwd);
    [A0_Bwd,nT_Bwd,Ea_Bwd]=fit_arrh(T,k_Bwd);
else
    A0_Fwd=NaN;
    nT_Fwd=NaN;
    Ea_Fwd=NaN;
    A0_Bwd=NaN;
    nT_Bwd=NaN;
    Ea_Bwd=NaN;
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
R_out.Input.ratio_H_H2=ratio_H_H2;
R_out.Input.ratio_CH4_H2=ratio_CH4_H2;
R_out.Input.ratio_CH3_CH4=ratio_CH3_CH4;

% Second, the gas properties are stored using the already-existing structure format for each gas
R_out.gas.H=H;
R_out.gas.H2=H2;
R_out.gas.CH3=CH3;
R_out.gas.CH4=CH4;

% Third, the "static" properties of reactants, products, images and TS are saved. Each subfield of TST is itself a structure, see
% above for their precise definition
R_out.TST.Reactants=s0;
R_out.TST.Products=sf;
R_out.TST.EnergyPathway=EnergyPathway;
R_out.TST.TS=TS;

% Finally, the reaction rate coefficients and related parameters are stored. The energies are expressed in eV, the effective
% reaction rate coefficient in s^-1, and the units of A and k depends on the presence or absence of a gas phase particle.
R_out.rates.Fwd.dE=E_Fwd/e;
R_out.rates.Fwd.kappa=kappa_Fwd;
R_out.rates.Fwd.A=A_Fwd;
R_out.rates.Fwd.k=k_Fwd;
R_out.rates.Fwd.n=n_Fwd;
R_out.rates.Fwd.r=r_Fwd;
R_out.rates.Bwd.dE=E_Bwd/e;
R_out.rates.Bwd.kappa=kappa_Bwd;
R_out.rates.Bwd.A=A_Bwd;
R_out.rates.Bwd.k=k_Bwd;
R_out.rates.Bwd.n=n_Bwd;
R_out.rates.Bwd.r=r_Bwd;
% The fit is saved as well, but only if the temperature array is large enough to perform the fit.
if length(T)>5
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
save([path,'Reaction_',R_out.Input.Reaction_ID,'.mat'],'R_out')
