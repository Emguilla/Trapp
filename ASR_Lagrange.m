function H_corrected=ASR_Lagrange(H)
%==================================================================================================================================%
% ASR_Lagrange.m: enforce the acoustic sum rule (ASR) via Lagrange multipliers on a Hessian matrix (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (12/10/2026) - Creation
%==================================================================================================================================%
% args:
%   H:  Hessian matrix to which ASR should be enforced
%==================================================================================================================================%
Natoms=size(H,1)/3;
% Symmetrisation
H_corrected=0.5*(H+H');
% Construct the spatial projection operator M = J_N (x) I_3
J_N=ones(Natoms,Natoms)/Natoms;
M=kron(J_N,eye(3));
% Apply Lagrange projection: (I - M) * H_sym * (I - M)
I_3N=eye(3*Natoms);
P=I_3N-M;
H_corrected=P*H_corrected*P;
% Symmetrisation
H_corrected=0.5*(H_corrected + H_corrected');
end