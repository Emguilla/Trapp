%==================================================================================================================================%
% constant_fund.m:  Creation of a matlab module containing fundamental constants (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (20/08/2025) - Creation
%       author: EYG
%   version 0.2 (09/03/2026) - Removal of the part related to the periodic table to a dedicated script.
%       author: EYG
%==================================================================================================================================%
% The point of this script being to save variable into a module, everything is cleaned first
clear all
% Fundamental constants
h=6.62607015e-34;
kB=1.380649e-23;
Na=6.0221408e+23;
e=1.602176565e-19;
c=299792458;
uma=1.66053906660e-27;
G=6.67430e-11;
mu_0=1.25663706212;
eps_0=8.8541878128e-12;
m_e=9.1093837015e-31;
m_p=1.67262192369e-27;
m_n=1.67492749804e-27;
% Derived fundamental constants
hbar=h/(2*pi);
R=Na*kB;

% recording the variable into a module "constant_fund.mat", to be placed in the matlab path for easy access
save('constant_fund.mat')