function R_out=Reaction_interp(R_in,n,varargin)
%==================================================================================================================================%
% Reaction_interp.m:    Interpolation of a reaction structure between each images (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (01/09/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   R_in:   Reaction structure or POSCAR array
%   n:      number of intermediate POSCAR structures to generate between each POSCAR
%==================================================================================================================================%
subNEB=false;
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'idx'
                subNEB=true;
                idx=varargin{p+1};
                R_in.POSCAR=R_in.POSCAR(idx);
                R_in.CONTCAR=R_in.CONTCAR(idx);
                R_in.XDATCAR=R_in.XDATCAR(idx,:);
                R_in.energies=R_in.energies(idx,:);
                R_in.reaction_coordinates=R_in.reaction_coordinates(idx,:);
                k=0;
                for p=1:length(idx)
                    for q=1:length(R_in.energies(1,:))
                        Forces{p,q}=R_in.Forces{idx(p),q};
                    end
                end
                R_in.Forces=Forces;
                R_in.MaxForces=R_in.MaxForces(idx,:);
        end
    end
end

if ~isfield(R_in,'POSCAR')
    % set incremental counter for the number of images of the final array of POSCAR structures
    ki=1;
    R_out(1)=R_in(1);
    for k=1:length(R_in)-1
        % For successive POSCAR structures, get the displacement vector using ds_POSCAR.m
        [~,ds]=ds_POSCAR(R_in(k+1),R_in(k));
        % the positions between a set of position and the next are computed through an interpolation of the displacement vector
        for p=1:n+1
            ki=ki+1;
            R_out(ki)=R_in(k);
            R_out(ki).positions=R_in(k).positions-(p/(n+1))*ds;
        end
    end
else
    R_out.Title=R_in.Title;
    % Recursive call to get the POSCAR structure array case
    R_out.POSCAR=Reaction_interp(R_in.POSCAR,n);
    R_out.CONTCAR=Reaction_interp(R_in.CONTCAR,n);
    for q=1:length(R_in.energies(1,:))
        R_out.XDATCAR(:,q)=Reaction_interp(R_in.XDATCAR(:,q),n);
        ki=1;
        E(1,q)=R_in.energies(1,q);
        x(1,q)=R_in.reaction_coordinates(1,q);
        for p=1:length(R_in.POSCAR)-1
            dE=R_in.energies(p+1,q)-R_in.energies(p,q);
            dx=R_in.reaction_coordinates(p+1,q)-R_in.reaction_coordinates(p,q);
            for k=1:n+1
                ki=ki+1;
                E(ki,q)=R_in.energies(p,q)+(k/(n+1))*dE;
                x(ki,q)=R_in.reaction_coordinates(p,q)+(k/(n+1))*dx;
            end
        end
    end
    R_out.energies=E;
    R_out.reaction_coordinates=x;
    % Forces and MaxForces field are not to be interpolated. For one, it would be physically irrelevant, and second it allows to
    % recognize an interpolated reaction from an original
    R_out.Forces=NaN;
    R_out.MaxForces=NaN;
end

end