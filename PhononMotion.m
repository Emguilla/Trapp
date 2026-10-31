function [ds,POSCARs]=PhononMotion(V,vibidx,varargin)
%==================================================================================================================================%
% PhononMotion.m:   Calculation of the displacement associated with the vibidx-th of a vibrational calculation (v0.1) 
%==================================================================================================================================%
% Version history:
%   version 0.1 (31/10/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   V:                  Matrix of the eigenvectors of the dynamical matrix
%   vibidx:             index of the mode for which the user want to get the delta_s
%   opt. args:          'POSCAR', followed by the POSCAR structure (or the path leading to it)
%                           (default: NaN, and POSCARs output is set to NaN)
%                       'nframes', number of elements in the POSCARs output.
%                           (default: 36)
%                       'amplitude', followed by the absolute value of the amplitude for the displacement of the atoms 
%                           corresponding to the mode vibidx.
%                           (default: 1)
%==================================================================================================================================%

% Set default and initial parameters
POSCAR=NaN;
nframes=36;
amplitude=1;

% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'poscar'
                POSCAR=varargin{p+1};
            case 'nframes'
                nframes=varargin{p+1};
            case 'amplitude'
                amplitude=varargin{p+1};
            otherwise
                error('Unknown argument !')
        end
    end
end

% Extraction of the displacement vector
ds=V(:,vibidx)/norm(V(:,vibidx));
ds=reshape(ds,3,length(V(1,:))/3)';
ds=ds/norm(ds);

% If the POSCAR input is a character string, the corresponding file is read into a POSCAR structure
if ischar(POSCAR)
    POSCAR=readPOSCAR(POSCAR);
end

% If the initial POSCAR is provided, the motion related to the phonon is computed and stored in an array of POSCAR structures.
if isstruct(POSCAR)
    for p=1:nframes
        phase=2*pi*(p-1)/nframes;
        POSCARs(p)=POSCAR;
        POSCARs(p).positions=POSCAR.positions+amplitude*sin(phase)*ds;
    end
else
    POSCARs=NaN;
end
end