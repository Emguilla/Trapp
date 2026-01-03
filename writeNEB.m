function writeNEB(POSCARs,dir)
%==================================================================================================================================%
% readPOSCAR.m: write a POSCAR MatLab structure in a file following VASP format (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (04/01/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   POSCAR:     Array of POSCAR structures to be written
%   dir:        name of the directory where the NEB must be written
%==================================================================================================================================%
n_images=length(POSCARs);
mkdir(dir)
format_folder=['%0',num2str(ceil(log10(n_images+1))+1),'i'];
for p=1:n_images
    mkdir([dir,'/',num2str(p-1,format_folder)])
    writePOSCAR(POSCARs(p),[dir,'/',num2str(p-1,format_folder),'/POSCAR'])
end