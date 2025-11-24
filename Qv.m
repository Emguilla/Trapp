function Q=Qv(F,T)
%==================================================================================================================================%
% Qv.m: Calculation of the vibration partition function as an array (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   F:  Frequency array
%   T:  Temperature in K
%==================================================================================================================================%
load('constant_fund.mat','kB','h')
% The vibrational partition function is stored as an array, each value corresponding to a specific vibrational mode (to avoid
% overflow errors)
for p=1:length(F)
    Q(p)=1./(1-exp(-h*F(p)/(kB*T)));
end
end