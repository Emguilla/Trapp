function [nu_f,nu_vib,nu_rot]=VibRot_Freq(VIB,T,varargin)
%==================================================================================================================================%
% VibRot_Freq.m:    calculation of the pseudo-vibrational frequency corresponding to the rotational partition function of a 
%                   rotor-like adsorbate on a surface (v0.3.0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (19/02/2026) - Creation by merging hr_freq.m and fr_freq.m
%       author: EYG
%   version 0.2 (03/07/2026) - Use of the new VIB structure which are processed from the raw phonon calculations. In addition, the 
%       author: EYG             user can now override the type of rotor by setting mode to free or hindered for the whole
%                               calculation.
%   version 0.3 (06/09/2026) - The default behaviour now returns the frequency how the description that minimise the partition
%       author: EYG             function. However, there is still the option to override and inforce a free or hindered rotor.
%   version 0.3.0.1 (12/10/2026) - Bug fix: "mode_set" was spelled out without the underscore sign.
%       contrib: EYG
%==================================================================================================================================%
% args:
%   nu_i:   Array of frequencies
%   path:   Location of the directory where the free_rotor.dat file is stored (if no file is there, the function returns the input
%           array of frequencies
%   T:      Temperature in Kelvin
%==================================================================================================================================%
load('constant_fund.mat','kB','h','e')
mode=NaN;
mode_set=false;
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'mode'
                mode_set=true;
                switch varargin{p+1}
                    case {'free',1}
                        mode=1;
                    case {'hindered',2}
                        mode=2;
                    otherwise
                        mode_set=false;
                        warning('Unknown entry for "mode" optional argument! Default back to "mode_set=false".')
                end
        end
    end
end

% The output array starts as identical as the input array
nu_i=VIB.nu;
nu_f=VIB.nu;
nu_vib=NaN;
nu_rot=NaN;
% If there exists a VibRotor.dat file in the directory specified by path, the first value is read as the index of the vibrational 
% mode to be replaced, the second value as the n-fold rotation number, the third as the moment of inertia and the fourth as the 
% rotational energy barrier
if isstruct(VIB.Rotor)
    idx=VIB.Rotor.Vib_idx;
    n=VIB.Rotor.sigma;
    I=VIB.Rotor.I;
    W=VIB.Rotor.RotEnergyBarrier*e;
    nu_vib=nu_f(length(nu_i)+1-idx);
    Q=sqrt(8*pi^3*I*kB*T/h^2)/n;
    % Calculation of the 1-D rotational partition function  
    % Calculation of the equivalent pseudo-frequency that produce the same partition function
    nu_rot_free=-kB*T*log(1-1/Q)/h;
    % Calculation of the roto-vibrational frequency (McClurg et al.)
    nu_rot_hind=n*sqrt(0.5*W/I)/(2*pi);
    W_kBT=(nu_rot_hind*2*pi/n)^2*I/(2*kB);
    Qv_free=Qv(nu_rot_free,T);
    Qv_hind=Qv(nu_rot_hind,T);
    if mode_set
        switch mode
            case{1,'free'}
                nu_rot=nu_rot_free;
            case{2,'hindered'}
                nu_rot=nu_rot_hind;
            case{3,'auto'}
        end
    else
        if Qv_free<Qv_hind
            nu_rot=nu_rot_free;
        else
            nu_rot=nu_rot_hind;
        end
    end
    nu_f(length(nu_i)+1-idx)=nu_rot;
end
end