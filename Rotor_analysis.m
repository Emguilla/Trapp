function ROTOR=Rotor_analysis(POSCAR,idx_origin,idx_rotor,sigma,varargin)
%==================================================================================================================================%
% Rotor_analysis.m: Post-processing of a rotation calculation and creation of a ROTOR structure (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (28/10/2026) - Creation
%       author: EYG
%   version 0.1.1 (29/10/2026) - Bug fix: there was an error wrt isnan for "idx_n" check when user-defined with two indices.
%       contrib: EYG
%==================================================================================================================================%
% args:
%   POSCAR:         POSCAR structure or string pointing to the location of the POSCAR file to be read as a POSCAR structure
%   idx_origin:     index of the anchor/origin atom
%   idx_rotor:      index of the atom belonging to the rotor.
%   sigma:          number of equivalent positions around a 360 rotation of the rotor (e.g., sigma of CH3 is equal to 3)
%   opt. args:          'path', followed by the path to the NEB directory
%                           (default: current directory)
%                       'idx_n', followed by the index of the second atom that defines, with the anchor atoms, the rotation axis. 
%                           Alternatively, two indices can be provided, and the vector between the corresponding atoms will define
%                           the rotation axis
%                           (default: NaN)
%                       'n', followed by a 3x1 array, defines the rotation axis directly
%                           (default: [NaN NaN NaN])
%                       'save_data', followed by the filename where to save the "Rotor" structure.
%==================================================================================================================================%

% Set default and initial parameters
path='./';
idx_n=NaN;
n=[NaN NaN NaN];

% Cleanup of the inputs
idx_rotor=unique(idx_rotor);
idx_origin=unique(idx_origin);
if ~isscalar(idx_origin)
    error('More than one index provided for the origin (anchor) atom')
end

% If the POSCAR input is a character string, the corresponding file is read into a POSCAR structure
if ischar(POSCAR)
    POSCAR=readPOSCAR([path,POSCAR]);
end

% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'path'
                path=varargin{p+1};
                if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
                    path=[path,'/'];
                end
            % sets the index of the atom(s) defining the rotation axis
            case 'idx_n'
                idx_n=varargin{p+1};
            % sets the rotation axis directly
            case 'n'
                n=varargin{p+1};
            case 'save_data'
                save_data=true;
                filename=varargin{p+1};
            otherwise
                error('Unknown argument !')
        end
    end
end


% Check existence of the OUTCAR resulting from the static calculation
ldir=dir([path,'*/POSCAR']);
for p=1:length(ldir)
    if ~exist([ldir(p).folder,'/OUTCAR'],'file')
        warning(['Nothing detected in image ',num2str(p),'!'])
        error('Please generate first the rotational images and perform a static calculation')
    end
end

% Extract all energies during the rotation as if it were a NEB calculation (which it technically is, but a static one)
NEB=NEB_analysis('path',path);

% Extraction of the origin coordinates
origin=POSCAR.positions(idx_origin,:);

% extraction of the rotation axis vector. If neither "idx_n" nor "n" are set by the user, the default is to set the rotation axis 
% as the vector from the anchor/origin to the CM of the rotor (excluding the anchor/origin atom).
if any(isnan(idx_n))&&any(isnan(n))
    idx_rotor_short=idx_rotor(~(idx_rotor==idx_origin));
    CM_rotor_short=sum(POSCAR.positions(idx_rotor_short,:))/length(idx_rotor_short);
    n=CM_rotor_short-POSCAR.positions(idx_origin,:);
else
    if isscalar(idx_n)
        n=POSCAR.positions(idx_n,:)-POSCAR.positions(idx_origin,:);
    elseif length(idx_n)==2
        n=POSCAR.positions(idx_n(2),:)-POSCAR.positions(idx_n(1),:);
    else
        error('More than one index were provided for the origin (anchor) atom')
    end
end

% Extraction of the rotor from the POSCAR
Rotor_0=subPOSCAR(POSCAR,idx_rotor);
% Rotation of the rotor to match the rotational axis and the z-axis, and translation of the anchor to the origin
[Rotor_z,~,~,~]=RotorProjectionZ(origin,n,Rotor_0);
% Computation of the inertia moment
[~,~,Iz]=InertiaMoment(Rotor_z);

% Recording of the inputs and outputs in a Rotor structure 
ROTOR.INPUT.idx_origin=idx_origin;
ROTOR.INPUT.idx_rotor=idx_rotor;
ROTOR.INPUT.idx_n=idx_n;
ROTOR.INPUT.n=n;
ROTOR.INPUT.sigma=sigma;
ROTOR.POSCAR=POSCAR;
ROTOR.Iz=Iz;
ROTOR.NEB=NEB;

% Save Rotor structure if requested
if save_data
    save([path,filename],'ROTOR')
end
