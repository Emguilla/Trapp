function [nu_f,nu_vib,nu_rot]=VibRot_Freq(VIB,T,varargin)
%==================================================================================================================================%
% VibRot_Freq.m:    calculation of the pseudo-vibrational frequency corresponding to the rotational partition function of a 
%                   rotor-like adsorbate on a surface (v0.4)
%==================================================================================================================================%
% Version history:
%   version 0.1 (19/02/2026) - Creation by merging hr_freq.m and fr_freq.m
%       author: EYG
%   version 0.2 (03/07/2026) - Use of the new VIB structure which are processed from the raw phonon calculations. In addition, the 
%       author: EYG             user can now override the type of rotor by setting mode to free or hindered for the whole
%                               calculation.
%   version 0.3 (06/09/2026) - The default behaviour now returns the frequency how the description that minimise the partition
%       author: EYG             function. However, there is still the option to override and enforce a free or hindered rotor.
%   version 0.3.0.1 (12/10/2026) - Bug fix: "mode_set" was spelled out without the underscore sign.
%       contrib: EYG
%   version 0.4 (31/10/2026) - Modification to accomodate the fact that frequencies are now listed in increasing order due to the 
%       author: EYG             diagonalisation being performed within the toolbox rather than relying on HIVE. In addition, the 
%                               automatic mode is now set to 0 and remain the default mode
%==================================================================================================================================%
% args:
%   VIB:    Vibration structure on which the substitution should be performed
%   T:      Temperature in Kelvin
% opt. args:            'mode', followed by the mode of substitution:
%                           'auto' (or 0) for automatic choice between hindered or free rotor 
%                           'free' (or 1) to choose a free rotor description
%                           'hindered' (or 2) to choose a hindered rotor description
%                           (default: 0)
%==================================================================================================================================%
load('constant_fund.mat','kB','h','e')
mode=0;
mode_set=true;
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'mode'
                mode_set=true;
                switch varargin{p+1}
                    case {'auto',0}
                        mode=0;
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
% If the Rotor field of the VIB structure is itself a structure, its own fields are extract to get the free and hindered rotational 
% partition function
if isstruct(VIB.Rotor)
    idx=VIB.Rotor.vibidx;
    n=VIB.Rotor.INPUT.sigma;
    I=VIB.Rotor.Iz;
    NEB=VIB.Rotor.NEB;
    W=max(NEB.energies-NEB.energies(1))*e;
    nu_vib=nu_f(idx);
    Q=sqrt(8*pi^3*I*kB*T/h^2)/n;
    % Calculation of the 1-D rotational partition function  
    % Calculation of the equivalent pseudo-frequency that produce the same partition function
    nu_rot_free=-kB*T*log(1-1/Q)/h;
    % Calculation of the roto-vibrational frequency (McClurg et al.)
    nu_rot_hind=n*sqrt(0.5*W/I)/(2*pi);
    Qv_free=Qv(nu_rot_free,T);
    Qv_hind=Qv(nu_rot_hind,T);
    switch mode
        case{1,'free'}
            nu_rot=nu_rot_free;
        case{2,'hindered'}
            nu_rot=nu_rot_hind;
        case{0,'auto'}
            if Qv_free<Qv_hind
                nu_rot=nu_rot_free;
            else
                nu_rot=nu_rot_hind;
            end
    end
    nu_f(idx)=nu_rot;
end
end