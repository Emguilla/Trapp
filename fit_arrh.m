function [A0,n,Ea]=fit_arrh(T,k)
%==================================================================================================================================%
% fit_arrh.m:   Fit of a temperature-dependent reaction rate coefficient to an equation A0*T^n*exp(Ea/(kB*T)) (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   T:  Temperatures in Kelvin
%   k:  reaction rate coefficient as a function of the temperature array T
%==================================================================================================================================%
load('constant_fund.mat')
% Linearisation of the data
[x,y]=prepareCurveData(1000./T,log(k));
% fitting parameters
reaction_rate_fittype=fittype('a+b*log(1000/x)-c*x','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'});
myfitoptions=fitoptions(reaction_rate_fittype);
myfitoptions.StartPoint=[log(1e14) 0 4];
% fit using matlab built-in function
fOUT=fit(x,y,reaction_rate_fittype,myfitoptions);
% Extraction of the a, b and c coefficients of the linear equation
a=fOUT.a;
b=fOUT.b;
c=fOUT.c;
% De-linearisation of the parameters
A0=exp(a);
n=b;
Ea=1000*kB*c;
end