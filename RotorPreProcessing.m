function RotorPreProcessing(path,POSCAR,vibidx,idx_origin,idx_n,idx_rotor,sigma,varargin)
%==================================================================================================================================%
% RotorPreProcessing.m: Creation of the rotor file to be read by "hr_freq.m" or "fr_freq.m" (v0.4)
%==================================================================================================================================%
% Version history:
%   version 0.1 (21/08/2025) - Creation
%       author: EYG
%   version 0.2 (25/08/2025) - Add the possibility to specify a vector instead of the index of the second atom on the rotation axis
%       author: EYG
%   version 0.3 (19/02/2026) - Removal of the mode argument. This argument is now an input from the user based on its analysis of
%       author: EYG             the relative height of the rotational energy barrier. The vibrational mode to be replaced by the 
%                               VibRot mode is provided via the new input vibidx. For sanity purposes, HIVE_analysis is called to
%                               generate a rendering of that particular mode.
%   version 0.4 (03/07/2026) - Add optional arguments, which can enforce the mode (free/hindered) to avoid user input. Other 
%       author: EYG             optional arguments relates to the rendering.
%==================================================================================================================================%
% args:
%   mode:       'hindered' -> create a "hindered_rotor.dat" file.
%               'free' -> create a "free_rotor.dat" file.
%   path:       path to the POSCAR file and where the file are to be written
%   vibidx:     Index of the vibrational mode of interest, starting from the lowest frequencies.
%   idx_origin: index (in the POSCAR structure) of the atom used as the anchor of the rotor adsorbate on the surface
%   idx_n:      integer scalar - indices (in the POSCAR structure) of the second atom that defines, along with the anchor atom at 
%                   the origin, the rotation axis of the adsorbate.
%   idx_rotor:      indices (in the POSCAR structure) of all the atoms that are part of the rotor, including anchor and origin.
%   sigma:          number of equivalent positions around a 360 rotation of the rotor (e.g., sigma of CH3 is equal to 3)
%   opt. args:          'quadview', followed by true or false to display 4 view angle at once
%                           (default: false)
%                       'rotate_view_angle', followed by a vector [theta,phi] to rotate the camera from its default direction
%                           (default: the model is shown from the spherical coordinates [19.5,9])
%                       'mode', followed by "free" or "hindered".
%                           (default to user input)
%                       'xyzlim', followed by a 3x2 matrix to give the limits of the rendering (in direct coordinates).
%                           (default: [0 1;0 1;0 1]).
%==================================================================================================================================%

% set default and initial parameters
rotate_view_angle=[0 0];
default_view_angle=[19.5 9];
quadview=false;
mode_set=false;
XYZLim=[0 1;0 1;0 1];

% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'quadview'
                quadview=varargin{p+1};
            case 'rotate_view_angle'
                default_view_angle=[0 0];
                rotate_view_angle=varargin{p+1};
            case 'mode'
                mode_set=true;
                switch varargin{p+1}
                    case {'free',1}
                        mode=1;
                    case {'hindered',2}
                        mode=2;
                    otherwise
                        mode_set=false;
                        warning('Unknown entry for "mode" optional argument! Default back to "mode_set=false".')
                end
            case 'xyzlim'
                set_lim=true;
                XYZLim=varargin{p+1};
        end
    end
end


if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end
% If the POSCAR input is a character string, the corresponding file is read into a POSCAR structure
if ischar(POSCAR)
    POSCAR=readPOSCAR([path,POSCAR]);
end
idx=unique(idx_rotor);
Rotor_0=subPOSCAR(POSCAR,idx_rotor);
aLim=[min(Rotor_0.xred(:,1)) max(Rotor_0.xred(:,1))];
bLim=[min(Rotor_0.xred(:,2)) max(Rotor_0.xred(:,2))];
cLim=[min(Rotor_0.xred(:,3)) max(Rotor_0.xred(:,3))];

if exist([path,'/DYNMAT'],'file')
    POSCARs=HIVE_analysis(vibidx,'path',path);
    for p=1:length(POSCARs)
        CONTCARs(p)=subPOSCAR(POSCARs(p),idx_rotor);
        aLim=[min([aLim(1) min(CONTCARs(p).xred(:,1))]) max([aLim(2) max(CONTCARs(p).xred(:,1))])];
        bLim=[min([bLim(1) min(CONTCARs(p).xred(:,2))]) max([bLim(2) max(CONTCARs(p).xred(:,2))])];
        cLim=[min([cLim(1) min(CONTCARs(p).xred(:,3))]) max([cLim(2) max(CONTCARs(p).xred(:,3))])];
    end
    ReactionRendering(CONTCARs,aLim,bLim,cLim,'quadview',true);
end
% extract the origin and rotation axis vector from POSCAR structure. if the "idx_n" variable is a vector, this vector is taken as
% the rotation axis
origin=POSCAR.positions(idx_origin,:);
if isscalar(idx_n)
    n=POSCAR.positions(idx_n,:)-POSCAR.positions(idx_origin,:);
elseif length(idx_n)==3
    n=idx_n;
end
% Creation of a sub-POSCAR structure containing only the rotor
idx=unique(idx_rotor);
Rotor_0=subPOSCAR(POSCAR,idx);
% Rotation of the rotor to match the rotational axis and the z-axis, and translation of the anchor to the origin
[Rotor_z,~,~,~]=RotorProjectionZ(origin,n,Rotor_0);
% Computation of the inertia moment
[~,~,Iz]=InertiaMoment(Rotor_z);

ldir=dir([path,'Rotation/*/POSCAR']);
if length(ldir)>1
    NEB=NEB_analysis('path',[path,'Rotation/']);
    ReactionRendering(NEB,XYZLim(1,:),XYZLim(2,:),XYZLim(3,:),'save_data','Infinite','rotate_view_angle',rotate_view_angle+default_view_angle,'quadview',quadview)
else
    error('Please generate first the rotational images and perform a static calculation')
end
for p=1:length(ldir)
    if ~exist([ldir(p).folder,'/OUTCAR'],'file')
        warning('Nothing detected in image ',num2str(p),'!')
        error('Please generate first the rotational images and perform a static calculation')
    end
end
W=max(NEB.energies)-min(NEB.energies);
if ~mode_set
    fprintf(['The rotational energy barrier seems to be about ',num2str(W*1000),' meV. Do you consider that to be a free or hindered rotation?\n 1) Free\n2) Hindered\n'])
    OK=false;
    while ~OK
        mode=input('==> ');
        switch mode
            case {1,2}
                OK=true;
            otherwise
                fprintf('Please enter 1 or 2.')
        end
    end
end
% For either mode, the value of the index of the vibrational mode, the symmetry number and inertia moment are written.
fid=fopen([path,'VibRotor.dat'],'w');
fprintf(fid,'%i',mode);
fprintf(fid,'   %04i',vibidx);
fprintf(fid,'    %03i',sigma);
fprintf(fid,'    %10.6e',Iz);
fprintf(fid,'    %10.6e',W);
fclose(fid);