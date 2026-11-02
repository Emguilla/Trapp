function writeDYNMAT(DYNMAT,filename,varargin)
%==================================================================================================================================%
% writeDYNMAT.m: write a DYNMAT MatLab structure in a file following VASP format (v0.1.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (01/11/2026) - Creation
%       author: EYG
%   version 0.1.1 (01/11/2026) - The default value for the mask is now set to the value of the mask contained in the DYNMAT
%       author: EYG                 structure.
%==================================================================================================================================%
% args:
%   DYNMAT:     DYNMAT structure to be written
%   filename:   path + name of the file in which the DYNMAT structure must be stored
%   opt. args:  'idx', followed by the atoms on which finite differences are applied.
%                   (default: all atoms are considered)
%               'mask', followed by a Nx3 logical array to specify which DOF are considered.
%                   (default: mask of the DYNMAT)
%==================================================================================================================================%

% Initialisation of the default parameters
idx_enabled=false;
mask_enabled=false;
mask=DYNMAT.mask;

% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'idx'
                idx_enabled=true;
                idx=varargin{p+1};
            case 'mask'
                mask_enabled=true;
                mask=varargin{p+1};
        end
    end
end

% If unspecified, generation of the mask of the DOF
if mask_enabled&&idx_enabled
    error('"idx" and "mask" arguments are mutually exclusive.')
elseif idx_enabled
    mask=mask*false;
    for p=idx
        mask(p,1:3)=true;
    end
end

% File opening
fid=fopen(filename,'w');

% Writing of the header
fprintf(fid,'%5d',DYNMAT.Ntype);
fprintf(fid,'%5d',DYNMAT.Natoms);
fprintf(fid,'%5d\n',DYNMAT.NDOF);
for p=1:length(DYNMAT.masses)
    fprintf(fid,'%8.3f',DYNMAT.masses(p));
end
fprintf(fid,'\n');

% Writing of the finite differences
for p=1:DYNMAT.Natoms
    for q=1:3
        if mask(p,q)
            fprintf(fid,'%5d%5d',[p,q]);
            fprintf(fid,'%8.4f%8.4f%8.4f\n',DYNMAT.dh{p,q});
            for r=1:DYNMAT.Natoms
                fprintf(fid,'%10.6f %10.6f %10.6f\n',DYNMAT.dF{p,q}(r,:));
            end
        end
    end
end
fclose(fid);