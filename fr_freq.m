function [nu_f,free_rotor]=fr_freq(nu_i,path,T)
%==================================================================================================================================%
% fr_freq.m:    calculation of the pseudo-vibrational frequency corresponding to the rotational partition function of an adsorbate
%               like a CH3 on a surface (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.2 (25/08/2025) - Add a logical variable assessing whether the last frequency has changed
%       author: EYG
%==================================================================================================================================%
% args:
%   nu_i:   Array of frequencies
%   path:   Location of the directory where the free_rotor.dat file is stored (if no file is there, the function returns the input
%           array of frequencies
%   T:      Temperature in Kelvin
%==================================================================================================================================%
load('constant_fund.mat','kB','h')
% The output array starts as identical as the input array
nu_f=nu_i;
% If there exists a free_rotor.dat file in the directory specified by path, the first value is read as sigma and the second as the
% moment of inertia
if exist([path,'free_rotor.dat'],'file')
    data=load([path,'free_rotor.dat']);
    sigma=data(1);
    I=data(2);
    % Calculation of the 1-D rotational partition function  
    Q=sqrt(8*pi^3*I*kB*T/h^2)/sigma;
    % Calculation of the equivalent pseudo-frequency that produce the same partition function
    nu_f(end)=-kB*T*log(1-1/Q)/h;
    % If there is a modification of the last frequency, the hind_rotor variable is set to true. If not, it is set to false
    free_rotor=true;
else
    free_rotor=false;
end
end