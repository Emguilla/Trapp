function proj=mode_projection(POSCARs,CONTCARs)
%==================================================================================================================================%
% mode_projection.m:    Computation of the projection of the trajectories contained in the POSCARs and CONTCARs structure arrays
%                       (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (28/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   POSCARs, CONTCARs:  Arrays of POSCAR structures
%==================================================================================================================================%

% set a threshold for what is considered a movement
tol_mov=1e-4;

for p=1:length(POSCARs)-1
    % Compute 3N displacement vectors
    [~,ds1]=ds_POSCAR(POSCARs(p),POSCARs(p+1));
    % Discart below threshold tol_mov
    ds1=ds1.*(abs(ds1)>tol_mov);
    % Normalise
    u1=ds1/norm(ds1);
    for q=1:length(CONTCARs)-1
        % Compute 3N displacement vectors
        [~,ds2]=ds_POSCAR(CONTCARs(q),CONTCARs(q+1));
        % Discart below threshold tol_mov
        ds2=ds2.*(abs(ds2)>tol_mov);
        % Normalise
        u2=ds2/norm(ds2);
        % Project
        proj(p,q)=abs(dot(u1(:),u2(:)))/norm(u1(:))/norm(u2(:));
    end
end
