function [SIGERR,msg_err]=compatPOSCARs(P1,P2)
%==================================================================================================================================%
% compatPOSCAR.m:   Check compatibility between POSCAR structures for merge/comparison operations (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (27/11/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   P1, P2: POSCAR structures
%==================================================================================================================================%

% Initialisation of default parameters
SIGERR=0;
msg_err='';

% check POSCAR structures for number of atoms and types of atoms
if length(P1.Z)~=length(P2.Z)
    SIGERR=1;
    msg_err='POSCARs do not contain the same number of atoms.';
else
    if ~all(P1.Z==P2.Z)
        SIGERR=2;
        msg_err='POSCARs contain different types of atoms';
    end
end
