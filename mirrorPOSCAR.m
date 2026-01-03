function POSCAR=mirrorPOSCAR(POSCAR,axis,axis_position)
%==================================================================================================================================%
% mirrorPOSCAR.m: Apply a mirror symmetry operation across a plane normal to the specified axis located at a certain position (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (03/01/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   POSCAR:         POSCAR structure, or path+filename of a POSCAR file
%   axis:           sets the orientation of the inversion plane. Possible arguments are 'x', 'y', 'z' or 'a', 'b', 'c'.
%   axis_position:  position (along the specified axis) where lies the inversion plane
%==================================================================================================================================%

% If the POSCAR input is a character string, the corresponding file is read into a POSCAR structure
if ischar(POSCAR)
    POSCAR=readPOSCAR(POSCAR);
end

% For each possible orientation of the inversion plane, either the cartesian or direct coordinates are shifted to put the location 
% of the inversion plane at the origin, then inverse the positions, then shifted the atoms back into place. The direct or cartesian
% coordinates are then updated accordingly.
switch lower(axis)
    case 'x'
        POSCAR.positions=POSCAR.positions-[axis_position 0 0];
        POSCAR.positions(:,1)=-POSCAR.positions(:,1);
        POSCAR.positions=POSCAR.positions+[axis_position 0 0];
        for p=1:sum(POSCAR.n_chemicals)
            POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
        end
    case 'y'
        POSCAR.positions=POSCAR.positions-[0 axis_position 0];
        POSCAR.positions(:,2)=-POSCAR.positions(:,2);
        POSCAR.positions=POSCAR.positions+[0 axis_position 0];
        for p=1:sum(POSCAR.n_chemicals)
            POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
        end
    case 'z'
        POSCAR.positions=POSCAR.positions-[0 0 axis_position];
        POSCAR.positions(:,3)=-POSCAR.positions(:,3);
        POSCAR.positions=POSCAR.positions+[0 0 axis_position];
        for p=1:sum(POSCAR.n_chemicals)
            POSCAR.xred(p,:)=(1/POSCAR.acell)*inv(POSCAR.vec')*POSCAR.positions(p,:)';
        end
    case 'a'
        POSCAR.xred=POSCAR.xred-[axis_position 0 0];
        POSCAR.xred(:,1)=-POSCAR.xred(:,1);
        POSCAR.xred=POSCAR.xred+[axis_position 0 0];
        for p=1:sum(POSCAR.n_chemicals)
            POSCAR.positions(p,:)=POSCAR.acell*POSCAR.vec'*POSCAR.positions(p,:)';
        end
    case 'b'
        POSCAR.xred=POSCAR.xred-[0 axis_position 0];
        POSCAR.xred(:,2)=-POSCAR.xred(:,2);
        POSCAR.xred=POSCAR.xred+[0 axis_position 0];
        for p=1:sum(POSCAR.n_chemicals)
            POSCAR.positions(p,:)=POSCAR.acell*POSCAR.vec'*POSCAR.positions(p,:)';
        end
    case 'c'
        POSCAR.xred=POSCAR.xred-[0 0 axis_position];
        POSCAR.xred(:,3)=-POSCAR.xred(:,3);
        POSCAR.xred=POSCAR.xred+[0 0 axis_position];
        for p=1:sum(POSCAR.n_chemicals)
            POSCAR.positions(p,:)=POSCAR.acell*POSCAR.vec'*POSCAR.positions(p,:)';
        end
end

