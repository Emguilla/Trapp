function gas=gas_modelling(press,V,T,molecules,concentrations)
%==================================================================================================================================%
% gas_modelling.m:  Calculation of the gas-related properties (total partition function, concentration) of different gas-phase
%                   species (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.2 (18/02/2026) - Complete overhaul of the code. The output is now an array of gas structures which fields contains 
%       author: EYG             the chemical formula of the gas species, its molar concentration, translational, rotational and
%                               vibrational partition function, its total partition function and the contribution to the exponential
%                               prefactor. The inputs still include the pressure, the volume, the temperature, but whereas
%                               previously the ratios between the different gas species were required, the input demands the
%                               chemical formulas and relative concentrations.
%==================================================================================================================================%
% args:
%   press:          Pressure in Pascals
%   V:              Volume in cubic meter
%   T:              Temperature in Kelvin
%   molecules:      String array of the formula of the gas species
%                   -> Only H, H2, CH3 and CH4 have been studied so far, carefully check results for other molecules
%   concentrations: Array of concentration of the corresponding chemicals in molecules
%==================================================================================================================================%
load('constant_fund.mat','R','Na','kB')
% Ideal gas law to determine the total number of moles
ntot=press*V./(R*T);
n_gas=length(molecules);
for p=1:n_gas
    % Extraction of the chemical formula of the gas species
    gas(p).chemicals=molecules{p};
    % Calculation of the relative ratio between all species, renormalisation of the relative concentration, and multiplication of 
    % the total number of mole to get the concentration of each species
    gas(p).n=ntot*concentrations(p)/sum(concentrations);
    % Computation of the translational partition functions
    gas(p).Qt=Qt(gas(p).chemicals,V,T);
    % Computation of the rotational partition functions
    gas(p).Qr=Qr(gas(p).chemicals,T);
    % Frequency of the vibrational modes of the molecules (IIRC those are experimental values, but I can't find the reference I 
    % used years ago)
    switch gas(p).chemicals
        case 'H'
            gas(p).nu=[];
        case 'H2'
            gas(p).nu=129.644861*1e12;
        case 'CH3'
            gas(p).nu=[97.077214; 94.086405; 81.898746;
                        40.730409; 16.535369; 11.913584]*1e12;
        case 'CH4'
            gas(p).nu=[92.769776; 91.681553; 87.871511;
                        76.795138; 43.368385; 39.476422;
                        32.049292; 27.151476; 10.182439]*1e12;
    end
    % Computation of the vibrational partition functions
    gas(p).Qv=prod(Qv(gas(p).nu,T));
    % Calculation of the zero-point energy
    gas(p).ZPE=ZPE(gas(p).nu);
    % Calculation of the total partition function
    gas(p).Q_tot=gas(p).Qt*gas(p).Qr*gas(p).Qv;
    % Calculation of the contribution to the exponential prefactor
    gas(p).A=Na*Q_tot*exp(gas(p).ZPE/(kB*T));
end