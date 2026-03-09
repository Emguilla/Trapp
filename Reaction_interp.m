function R_out=Reaction_interp(R_in,n,varargin)
%==================================================================================================================================%
% Reaction_interp.m:    Interpolation of a reaction structure between each images (v0.3)
%==================================================================================================================================%
% Version history:
%   version 0.1 (01/09/2025) - Creation
%       author: EYG
%   version 0.2 (27/01/2026) - The degree of freedom of the interpolated structures are now set to true if any POSCAR from the input
%       author: EYG             has that degree of freedom set to true.
%   version 0.3 (09/03/2026) - Add optional argument to repel atoms too close to each other using the "Proximity_check" function.
%       author: EYG
%==================================================================================================================================%
% args:
%   R_in:       Reaction structure or POSCAR array
%   n:          number of intermediate POSCAR structures to generate between each POSCAR
%   opt. args:          'idx', followed by the indices of the images to be used to perform the interpolation.
%                           (default: all indices are considered)
%                       'min', followed by the ratio wrt the maximum bond distance allowed between each pair of atoms. 
%                           (default: 0.82)
%                       'pushback', followed by true or false to specify whether pairs of too close atoms should be repelled.
%                           (default: false)
%==================================================================================================================================%

% Initialisation of the default parameters and handling of the optional arguments
subNEB=false;
pushback=false;
DBOND_MIN_ratio=0.82;

% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'min'
                DBOND_MIN_ratio=varargin{p+1};
            case 'pushback'
                pushback=varargin{p+1};
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

% Two cases are implemented: One where the input variable is a POSCAR array, the other a NEB structure. The latter is treated
% recursively to fall back on the former case
if ~isfield(R_in,'POSCAR')
    % set incremental counter for the number of images of the final array of POSCAR structures
    ki=1;
    R_out(1)=R_in(1);
    for k=1:length(R_in)
        if ~isnan(R_out(1).constraint(1))
            R_out(1).constraint=(R_out(1).constraint+R_in(k).constraint)>0;
        end
    end
    if any(isnan(R_out(1).constraint))
        R_out(1).constraint=ones(1,1)*NaN;
    end
    for k=1:length(R_in)-1
        % For successive POSCAR structures, get the displacement vector using ds_POSCAR.m
        [~,ds]=ds_POSCAR(R_in(k+1),R_in(k));
        % the positions between a set of position and the next are computed through an interpolation of the displacement vector
        for p=1:n+1
            ki=ki+1;
            R_out(ki)=R_in(k);
            R_out(ki).constraint=R_out(1).constraint;
            R_out(ki).positions=R_in(k).positions-(p/(n+1))*ds;
            % Call to the proximity_check function to repel, if relevant, the pairs of atoms that are too close
            if pushback
                [R_out(ki),~]=Proximity_check(R_out(ki),'min',DBOND_MIN_ratio,'pushback',pushback);
            end
            for p=1:sum(R_in(1).n_chemicals)
                R_out(ki).xred(p,:)=(1/R_out(ki).acell)*inv(R_out(ki).vec')*R_out(ki).positions(p,:)';
            end
        end
    end
else
    R_out.Title=R_in.Title;
    % Recursive call to get the POSCAR structure array case
    R_out.POSCAR=Reaction_interp(R_in.POSCAR,n,'min',DBOND_MIN_ratio,'pushback',pushback);
    R_out.CONTCAR=Reaction_interp(R_in.CONTCAR,n,'min',DBOND_MIN_ratio,'pushback',pushback);
    for q=1:length(R_in.energies(1,:))
        R_out.XDATCAR(:,q)=Reaction_interp(R_in.XDATCAR(:,q),n,'min',DBOND_MIN_ratio,'pushback',pushback);
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