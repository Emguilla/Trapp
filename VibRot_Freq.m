function [nu_f,nu_vib,nu_rot]=VibRot_Freq(nu_i,path,T)
%==================================================================================================================================%
% VibRot_Freq.m:    calculation of the pseudo-vibrational frequency corresponding to the rotational partition function of a 
%                   rotor-like adsorbate on a surface (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (19/02/2026) - Creation by merging hr_freq.m and fr_freq.m
%       author: EYG
%==================================================================================================================================%
% args:
%   nu_i:   Array of frequencies
%   path:   Location of the directory where the free_rotor.dat file is stored (if no file is there, the function returns the input
%           array of frequencies
%   T:      Temperature in Kelvin
%==================================================================================================================================%
load('constant_fund.mat','kB','h','e')
hind=false;
free=false;
% The output array starts as identical as the input array
nu_f=nu_i;
nu_vib=NaN;
nu_rot=NaN;
% If there exists a VibRotor.dat file in the directory specified by path, the first value is read as the index of the vibrational 
% mode to be replaced, the second value as the n-fold rotation number, the third as the moment of inertia and the fourth as the 
% rotational energy barrier
if exist([path,'VibRotor.dat'],'file')
    data=load([path,'VibRotor.dat']);
    mode=data(1);
    if mode==1
        idx=data(2);
        n=data(3);
        I=data(4);
        % Calculation of the 1-D rotational partition function  
        Q=sqrt(8*pi^3*I*kB*T/h^2)/n;
        % Calculation of the equivalent pseudo-frequency that produce the same partition function
        nu_vib=nu_f(length(nu_i)+1-idx);
        nu_f(length(nu_i)+1-idx)=-kB*T*log(1-1/Q)/h;
        nu_rot=nu_f(length(nu_i)+1-idx);
        % If there is a modification of a frequency, the free variable is set to true. If not, it stays to false
        free=true;
    elseif mode==2
        idx=data(2);
        n=data(3);
        I=data(4);
        W=data(5)*e;
        % Calculation of the roto-vibrational frequency (McClurg et al.)
        nu_f(idx)=n*sqrt(0.5*W/I)/(2*pi);
        % If there is a modification of a frequency, the hind variable is set to true. If not, it stays to false
        hind=true;
    end
end
end