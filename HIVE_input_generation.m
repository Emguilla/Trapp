function HIVE_input_generation(path,idx,varargin)
%==================================================================================================================================%
% HIVE_input_generation.m:  Generate input file for HIVE post-processing
%==================================================================================================================================%
% Version history:
%   version 0.1 (17/07/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   path:       path to the directory where the "inputphonon.dat" file must be written
%   idx:        index of the mode to analyse (in decreasing order of frequency magnitude)
%   opt. args:  'n_images', followed by the number of intermediate structures along the direction of the mode analysed 
%                   (default: n_images=36)
%               'clean3D', followed by true or false to specify whether a cleanup of the rotational mode must be performed
%                   (default: false)
%==================================================================================================================================%

% set default parameters
clean3D=false;
n_images=36;

% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'n_images'
                n_images=varargin{p+1};
            case 'clean3d'
                clean3D=varargin{p+1};
        end
    end
end

% Writing operation
fid=fopen([path,'/inputphonon.dat'],'w');
fprintf(fid,'general:\n');
fprintf(fid,'  datasource = VASP\n');
fprintf(fid,'  ncores = 10\n');
fprintf(fid,'  energy unit = THz\n');
fprintf(fid,'  supercell = .false.\n');
if clean3D
    fprintf(fid,'  cleanup = rot3D\n');
    fprintf(fid,'  rotationprojection = smallrot\n');
    fprintf(fid,'  ModifiedGramSchmidtAll = .TRUE.\n');
end
fprintf(fid,'  GammaPhonon = .TRUE.\n');
fprintf(fid,'SpecialQpoint\n');
fprintf(fid,'  singleQ = .TRUE.\n');
fprintf(fid,'  singleQcoord = 0.0  0.0  0.0\n');
if ischar(idx)
    fprintf(fid,['  writeQmodes = ',char(idx),'\n']);
else
    fprintf(fid,['  writeQmodes = ',num2str(idx),'\n']);
end
fprintf(fid,['  images = ',num2str(n_images),'\n']);
fprintf(fid,'  writeDisplacement = .FALSE.\n');
fclose(fid);
