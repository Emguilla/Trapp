function [nu_f,hind_rotor]=hr_freq(nu_i,path)
%==================================================================================================================================%
% hr_freq.m:    calculation of the pseudo-vibrational frequency corresponding to the rotational partition function of an adsorbate
%               like a CH3 on a surface (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.1.1 (21/08/2025) - swap of the position of the moment of inertia and the rotational energy barrier. The reason comes
%                                   from the fact that n and I are preprocessed, whereas W is manually entered by the user
%                                   afterwards
%   version 0.2 (25/08/2025) - Add a logical variable assessing whether the last frequency has changed
%       author: EYG
%==================================================================================================================================%
% args:
%   nu_i:   Array of frequencies
%   path:   Location of the directory where the free_rotor.dat file is stored (if no file is there, the function returns the input
%           array of frequencies
%==================================================================================================================================%
% The output array starts as identical as the input array
nu_f=nu_i;
% If there exists a free_rotor.dat file in the directory specified by path, the first value is read as the n-fold rotation number, 
% the second as the energy barrier for rotation and the third as the moment of inertia
if exist([path,'hindered_rotor.dat'],'file')
    data=load([path,'hindered_rotor.dat']);
    n=data(1);
    I=data(2);
    W=data(3);
    % Calculation of the roto-vibrational frequency (McClurg et al.)
    nu_f(end)=n*sqrt(0.5*W/I)/(2*pi);
    % If there is a modification of the last frequency, the hind_rotor variable is set to true. If not, it is set to false
    hind_rotor=true;
else
    hind_rotor=false;
end
end