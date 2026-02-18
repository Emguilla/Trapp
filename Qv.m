function Q=Qv(F,T)
%==================================================================================================================================%
% Qv.m: Calculation of the vibration partition function as an array (v0.1.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.1.1 (18/02/2026) - Added an exception resulting in Q=1 when the frequency array is empty.
%       contrib: EYG
%==================================================================================================================================%
% args:
%   F:  Frequency array
%   T:  Temperature in K
%==================================================================================================================================%
load('constant_fund.mat','kB','h')
% The vibrational partition function is stored as an array, each value corresponding to a specific vibrational mode (to avoid
% overflow errors)
if isempty(F)
    Q=1;
else
    for p=1:length(F)
        Q(p)=1./(1-exp(-h*F(p)/(kB*T)));
    end
end