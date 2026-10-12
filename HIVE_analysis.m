function [POSCARs,Freq]=HIVE_analysis(idx,varargin)
%==================================================================================================================================%
% HIVE_analysis.m:  Use of the HIVE program to get accurate frequencies and the POSCAR structures corresponding to the displacement
%                   in the direction of a specific mode. (v0.3.1)
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
%   version 0.2.4 (16/02/2026) - Correction to the path where the output files are located. There was a conflict between the 
%       contrib: EYG                directory change to the path and the call to the file in the same path, which doesn't exist at 
%                                   the location of the target folder
%   version 0.2.5 (19/02/2026) - Correction to the limits of the rendering. Previously the atoms could get out of the field of view
%       contrib: EYG                during the vibration displacement
%   version 0.2.6 (17/07/2026) - Add the possibility to require all modes to be written out. The input file for HIVE is now written 
%       contrib: EYG                by a separate function, that uses the new "clean3D" optional argument to clean rotational 
%                                   modes. If left empty, idx defaults to "1" instead of the previous value of the 3N-th mode.
%   version 0.3 (06/09/2026) - Add the possibility to save the POSCARs in a dedicated "POSCARs.mat" save file. Additional changes 
%       author: EYG             were made to accomodate for the use of HIVE4 instead of HIVE3. 
%   version 0.3.1 (12/10/2026) - The program now avoids to compute the displacement when the user enters "NaN" as the "idx" input.
%       contrib: EYG
%==================================================================================================================================%
% args:
%   idx:        index of the mode to analyse (in decreasing order of frequency magnitude), or idx="all" to write all modes.
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
clean3D=false;
n_images=36;
path='./';
rendering=false;
save_POSCAR=false;
save_data=false;
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'clean3d'
                clean3D=varargin{p+1};
            case 'n_images'
                n_images=varargin{p+1};
            case 'path'
                path=varargin{p+1};
                if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
                    path=[path,'/'];
                end
            case 'save_poscar'
                save_POSCAR=true;
                if ischar(idx)
                    if ~strcmpi(idx,'all')
                        error('The "save_POSCAR" argument can only be used when analysing all modes for the first time.')
                    end
                else
                    error('The "save_POSCAR" argument can only be used when analysing all modes for the first time.')
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
WSL_distro='Ubuntu'; % name of your linux distro
old_path=getenv('PATH');
if ~strcmpi(old_path(end),';')
    old_path(end+1)=';';
end
setenv('PATH',old_path);

% Path to your local copy of the hive3.exe file, do not forget the ending ";"
add_path='C:\Users\emgui\OneDrive-UNamur\Bureau\Research\HIVE4.x-tools-dev\source\CommandLineInterface\bin\Release;'; 

% % source variable within ubuntu
% system(['wsl -d ',WSL_distro,' /home/emguilla/.bash_profile']);

if ~exist([add_path(1:end-1),'/hive4.exe'],'file')
    error('hive4.exe cannot be found in the specified path')
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

% Remove older HIVE_analysis directory, if any
if exist('HIVE_analysis','dir')
    system(['wsl -d ',WSL_distro,' rm -rf HIVE_analysis']);
end

if ischar(idx)
    if exist('POSCARs.mat')&&save_POSCAR&&strcmpi(idx,'all')
        error('A prior analysis has been performed in this directory. Please manually remove the "POSCARs.mat" file if this move was really intended')
    end
end

% Create a new HIVE_analysis directory and copy the files inside
[~,~]=system(['wsl -d ',WSL_distro,' mkdir HIVE_analysis']);
[~,~]=system(['wsl -d ',WSL_distro,' cp POSCAR DYNMAT OUTCAR HIVE_analysis/.']);
cd HIVE_analysis
% Create the inputphonon.dat file to set all the HIVE parameters
if ischar(idx)||any(isnan(idx))
    HIVE_input_generation('.',1,'clean3D',clean3D,'n_images',n_images)
else
    HIVE_input_generation('.',idx(1),'clean3D',clean3D,'n_images',n_images)
end
% Run the linux HIVE executable from MatLab via WSL and extract frequencies from PHONONOUT.hive
[~,~]=system(['wsl -d ',WSL_distro,' bash -c "source /home/emguilla/.bashrc && hive4.exe phonons inputphonon.dat"']);
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
if ischar(idx)
    if strcmpi(idx,'all')
        clear idx
        idx=1:length(Freq);
        HIVE_input_generation('.','all','clean3D',clean3D,'n_images',n_images)
    else
        error('Unknown char value for "idx" argument. Possible options: idx="all"')
    end
    [~,~]=system(['wsl -d ',WSL_distro,' bash -c "source /home/emguilla/.bashrc && hive4.exe phonons inputphonon.dat"']);
else
    for q=1:length(idx)
        tic
        HIVE_input_generation('.',idx(q),'clean3D',clean3D,'n_images',n_images)
        [~,~]=system(['wsl -d ',WSL_distro,' bash -c "source /home/emguilla/.bashrc && hive4.exe phonons inputphonon.dat"']);
        toc
    end
end
q=1;
if ~isnan(idx)
    nbytes = fprintf(['Reading SQ_PHONON_',num2str(length(Freq)-idx(q)+1)],' out of %d',length(idx));
    for q=1:length(idx)
        fprintf(repmat('\b',1,nbytes))
        POSCARs(q,:)=readSQPHONON(['SQ_PHONON_',num2str(length(Freq)-idx(q)+1),'.hive']);
    end
end
%system(['wsl -d ',WSL_distro,' rm -rf POSCARs']);
cd ..
if save_POSCAR&&length(idx)==length(Freq)
    disp('Writing POSCARs ...')
    save('POSCARs.mat','POSCARs')
end
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