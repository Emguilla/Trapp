function [POSCARs,Freq]=HIVE_analysis(idx,varargin)
%==================================================================================================================================%
% HIVE_analysis.m:  Use of the HIVE program to get accurate frequencies and the POSCAR structures corresponding to the displacement
%                   in the direction of a specific mode. (v0.2.3)
%==================================================================================================================================%
% Version history:
%   version 0.1 (27/08/2025) - Creation
%       author: EYG
%   version 0.2 (28/08/2025) - Add 3D rendering (and its recording) options
%       author: EYG
%   version 0.2.1 (11/09/2025) - The frequencies are now extracted from the "Freq_HIVE.dat" file to accomodate for the changes in
%       contrib: EYG                the HIVE_analysis.m script
%   version 0.2.2 (08/12/2025) - For some reason, the index was considered an optional argument. This has been fixed.
%       contrib: EYG
%   version 0.2.3 (06/01/2026) - The mounting of the drive has been commented. The possibility to mount will be added as a proper
%       contrib: EYG                optional argument, or explained in the docs.
%   version 0.2.4 (26/02/2026) - Correction to the path where the output files are located. There was a conflict between the 
%       contrib: EYG                directory change to the path and the call to the file in the same path, which doesn't exist at 
%                                   the location of the target folder
%==================================================================================================================================%
% args:
%   idx:        index of the mode to analyse (in decreasing order of frequency magnitude)
%   opt. args:  'n_images', followed by the number of intermediate structures along the direction of the mode analysed 
%                   (default: n_images=36)
%               'path', followed by the path to the directory that contain the DYNMAT, OUTCAR and POSCAR files
%                   (default: '.')
%               'save_data', followed by the playing direction of the gif file. Possible options are 'Fwd', 'Bwd', 'Round', or 
%                   infinite
%                   (default: false)
%               'rendering', followed by the boundaries of the rendering in a 2x3 array
%                   (default: false)
%==================================================================================================================================%
n_images=36;
path='./';
rendering=false;
save_data=false;
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'n_images'
                n_images=varargin{p+1};
            case 'path'
                path=varargin{p+1};
                if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
                    path=[path,'/'];
                end
            case 'save_data'
                save_data=true;
                gif_direction=varargin{p+1};
            case 'rendering'
                rendering=true;
                XLim=varargin{p+1}(1,:);
                YLim=varargin{p+1}(2,:);
                ZLim=varargin{p+1}(3,:);
            otherwise
                warning(['Unknown argument: "',varargin{p},'"'])
        end
    end
end
curr_dir=pwd;

% In the event the hive3.exe file is not located on your Windows path (this is only a temporary change to your Windows path, it will
% be reset upon restarting/shutting off).
drive_letter='E'; % Pick the letter corresponding to that assigned to the USB key by your Windows installation (e.g. "E")
WSL_distro='Ubuntu-22.04'; % name of your linux distro
old_path=getenv('PATH');
if ~strcmpi(old_path(end),';')
    old_path(end+1)=';';
end
setenv('PATH',old_path);
add_path='C:\Users\emgui\OneDrive-UNamur\Cluster;'; % Path to your local copy of the hive3.exe file, do not forget the ending ";"
if ~exist([add_path(1:end-1),'/hive3.exe'],'file')
    error('hive3.exe cannot be found in the specified path')
end
% If the Windows path already includes the path to your local copy of the hive3.exe file, nothing is added to the path.
if contains(old_path,add_path(2:end-1))
    new_path=old_path;
else
    new_path=[old_path,add_path];
end
setenv('PATH',new_path);
%==================================================================================================================================%
% Type the following command (in your WSL installation) to allow any program to mount a drive without password privilege
% THIS IS A RISKY PROCEDURE, PROCEED WITH CAUTION !!! AND DO NOT PRECEDE THE FOLLOWING WITH SUDO !
% echo "`whoami` ALL=(ALL) NOPASSWD:ALL" | sudo tee /etc/sudoers.d/`whoami` && sudo chmod 0440 /etc/sudoers.d/`whoami`
%==================================================================================================================================%
cd(path)

% In case the analysis must be performed on another drive than the one which hosts the WSL installation (e.g. USB key), the drive
% must be mounted inside WSL. 
% [~,~]=system(['wsl -d ',WSL_distro,' sudo mount -t drvfs ',upper(drive_letter),': /mnt/',lower(drive_letter)]);
% it should look like : "wsl -d Ubuntu-22.04 sudo mount -t drvfs E: /mnt/e"

% Check that all necessary files are present in the current folder (i.e., check that the calculation ran properly)
if ~exist('POSCAR','file')
    error('POSCAR file is missing!')
else
    POSCAR=readPOSCAR('POSCAR');
end
if ~exist('OUTCAR','file')
    error('OUTCAR file is missing!')
end
if ~exist('DYNMAT','file')
    error('DYNMAT file is missing!')
end
if ~exist('idx','var')
    idx=sum(POSCAR.n_chemicals);
end

% Remove older HIVE_analysis directory, if any
if exist('HIVE_analysis','dir')
    system(['wsl -d ',WSL_distro,' rm -rf HIVE_analysis']);
end

% Create a new HIVE_analysis directory and copy the files inside
[~,~]=system(['wsl -d ',WSL_distro,' mkdir HIVE_analysis']);
[~,~]=system(['wsl -d ',WSL_distro,' cp POSCAR DYNMAT OUTCAR HIVE_analysis/.']);

cd HIVE_analysis
% Create the inputphonon.dat file to set all the HIVE parameters
fid=fopen('inputphonon.dat','w');
fprintf(fid,'general:\n');
fprintf(fid,'  datasource = VASP\n');
fprintf(fid,'  ncores = 10\n');
fprintf(fid,'  energy unit = THz\n');
fprintf(fid,'  supercell = .false.\n');
fprintf(fid,'  GammaPhonon = .TRUE.\n');
fprintf(fid,'SpecialQpoint\n');
fprintf(fid,'  singleQ = .TRUE.\n');
fprintf(fid,'  singleQcoord = 0.0  0.0  0.0\n');
fprintf(fid,['  writeQmodes = ',num2str(idx),'\n']);
fprintf(fid,['  images = ',num2str(n_images),'\n']);
fprintf(fid,'  writeDisplacement = .FALSE.\n');
fclose(fid);

% Run the linux HIVE executable from MatLab via WSL and extract frequencies from PHONONOUT.hive
[~,~]=system(['wsl -d ',WSL_distro,' hive3.exe phonons inputphonon.dat']);
grep('PHONONOUT.hive',')    ','prt','Freq_HIVE.dat');
% The value of the frequencies are extracted from the Freq_HIVE.dat file
k=1;
Data=readlines('Freq_HIVE.dat');
while ~isempty(Data{k})
    Freq(k)=str2double(Data{k}(11:21))*1e12;
    k=k+1;
end
copyfile PHONONOUT.hive ..

% HIVE standard output for the phonon calculation includes the POSCAR of the system along the direction of the vibrational mode
% specified by idx. Each POSCAR_XXXXX_N.vasp is read and stored in an array of POSCAR structures. Note that since there is not
% a standard format for the XXXXX part of the file names, a search from the end for the last underscore is first performed.
ldir=dir('*.vasp');
l=length(ldir(1).name);
while ~strcmpi(ldir(1).name(l),'_')
    l=l-1;
end
for p=1:n_images
    POSCARs(p)=readPOSCAR([ldir(p).name(1:l),num2str(p),'.vasp']);
end

cd ..
% clean the temporary folder
[~,~]=system(['wsl -d ',WSL_distro,' rm -rf HIVE_analysis']);
% choose to render the vibrational mode (and save_data its rendering if requested)
if rendering
    if save_data
        ReactionRendering(POSCARs,XLim,YLim,ZLim,'save_data',gif_direction);
    else
        ReactionRendering(POSCARs,XLim,YLim,ZLim);
    end
end
cd(curr_dir)