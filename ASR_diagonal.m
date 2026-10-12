function H_corrected=ASR_diagonal(H)
%==================================================================================================================================%
% ASR_diagonal.m: enforce the acoustic sum rule (ASR) on the diagonal element of a Hessian matrix (v0.1)
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
for u=1:Natoms
    idx_u=(3*(u-1)+1):(3*u);
    % Compute residual sum over all columns for atom u
    residual_u=zeros(3,3);
    for v=1:Natoms
        idx_v=(3*(v-1)+1):(3*v);
        residual_u=residual_u+H_corrected(idx_u,idx_v);
    end
    % Subtract residual strictly from the self-block (u,u)
    H_corrected(idx_u,idx_u)=H_corrected(idx_u,idx_u)-residual_u;
end
% Symmetrisation
H_corrected=0.5*(H_corrected+H_corrected');
end