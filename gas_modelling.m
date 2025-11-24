function [H,H2,CH3,CH4]=gas_modelling(press,V,T,ratio_H_H2,ratio_CH4_H2,ratio_CH3_CH4)
%==================================================================================================================================%
% gas_modelling.m:  Calculation of the gas-related properties (total partition function, concentration) of different gas-phase
%                   species (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   press:                                      Pressure in Pascals
%   V:                                          Volume in cubic meter
%   T:                                          Temperature in Kelvin
%   ratio_H_H2, ratio_CH4_H2, ratio_CH3_CH4:    Ratio between the different species.
%                                               -> Only H, H2, CH3 and CH4 are currently implemented
%==================================================================================================================================%
load('constant_fund.mat','R')
% Frequency of the vibrational modes of the H2, CH3 and CH4 molecules (IIRC those are experimental value, but I can't find the
% reference I used years ago
nu_H2=129.644861*1e12;
nu_CH3=[97.077214;94.086405;81.898746;40.730409;16.535369;11.913584]*1e12;
nu_CH4=[92.769776;91.681553;87.871511;76.795138;43.368385;39.476422;32.049292;27.151476;10.182439]*1e12;

% Ideal gas law to determine the total number of moles
ntot=press*V/(R*T);

% Calculation of the relative ratio between all species
cH2=1;
cH=cH2*ratio_H_H2;
cCH4=cH2*ratio_CH4_H2;
cCH3=min([cH cCH4])*ratio_CH3_CH4;

% Renormalisation of the relative concentration, and multiplication the total number of mole to get the concentration of each
% species
H.n=ntot*cH/(cH2+cH+cCH4+cCH4);
H2.n=ntot*cH2/(cH2+cH+cCH4+cCH4);
CH3.n=ntot*cCH3/(cH2+cH+cCH4+cCH4);
CH4.n=ntot*cCH4/(cH2+cH+cCH4+cCH4);

% Computation of the translational partition functions
H.Qt=Qt('H',V,T);
H2.Qt=Qt('H2',V,T);
CH3.Qt=Qt('CH3',V,T);
CH4.Qt=Qt('CH4',V,T);

% Computation of the rotational partition functions
H2.Qr=Qr('H2',T);
CH3.Qr=Qr('CH3',T);
CH4.Qr=Qr('CH4',T);

% Computation of the vibrational partition functions
H2.Qv=prod(Qv(nu_H2,T));
CH3.Qv=prod(Qv(nu_CH3,T));
CH4.Qv=prod(Qv(nu_CH4,T));

% Calculation of the zero-point energies
H2.ZPE=ZPE(nu_H2);
CH3.ZPE=ZPE(nu_CH3);
CH4.ZPE=ZPE(nu_CH4);

% Calculation of the total partition function
H.Q_tot=H.Qt;
H2.Q_tot=H2.Qt*H2.Qr*H2.Qv;
CH3.Q_tot=CH3.Qt*CH3.Qr*CH3.Qv;
CH4.Q_tot=CH4.Qt*CH4.Qr*CH4.Qv;
end