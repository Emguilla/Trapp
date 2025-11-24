function kappa=SkodjeTruhlar(T,nu,E)
%==================================================================================================================================%
% SkodjeTruhlar.m:  Computation of the Skodje-Thrular coefficient to account for quantum tunneling through the energy barrier (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   T:  Temperature in K
%   nu: Frequency of the transition state across the saddle point
%   E:  Energy barrier height
%==================================================================================================================================%
load('constant_fund.mat','kB','h')
alpha=2*pi/(h*nu);
beta=1/(kB*T);
if alpha>beta
    kappa=(beta*pi/alpha)/sin(beta*pi/alpha)-beta/(alpha-beta)*exp((beta-alpha)*E);
else
    kappa=beta/(alpha-beta)*(exp((beta-alpha)*E)-1);
end
