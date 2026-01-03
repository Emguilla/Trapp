function EnergyPathway=NEB_analysis(varargin)
%==================================================================================================================================%
% NEB_analysis.m:   Post-processing of a (c)NEB calculation (v0.1.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (02/09/2025) - Creation using bits and pieces from my thesis works
%       author: EYG
%   version 0.1.1 (16/12/2025) - Correction of the call to readPOSCAR to account for the removal of optional argument 'XDATCAR'
%       contrib: EYG
%   version 0.1.2 (03/01/2026) - The recursive call only occurs on folders that starts with "subNEB" instead of "sub".
%       contrib: EYG
%==================================================================================================================================%
% args:
%   opt. args:          'path', followed by the path to the NEB directory
%                           (default: current directory)
%                       'save_data', followed by the filename where the NEB structure must be saved
%                           (default: nothing is saved)
%                       'coordinates_mapping', followed by the type of calculation of the reaction coordinates (possible choices 
%                           are 'parametric' or 'projected')
%                           (default: parametric)
%                       'projected_resolution', followed by the resolution of the projection onto the interpolation between states
%                           (default: 1e-3)
%                       'subNEB', followed by true or false to allow for recursive search of subNEB calculations
%                           (default: false)
%==================================================================================================================================%
% set default and initial parameters
save_data=false;
parametric_coordinates=true;
projected_coordinates=false;
projected_resolution=0.001;
subNEB=false;
idx_up=[NaN NaN];
filename='EP.mat';
curr_dir=pwd;
halted_calculation=false;
path='./';
% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'path'
                path=varargin{p+1};
                if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
                    path=[path,'/'];
                end
            case 'save_data'
                save_data=true;
                filename=varargin{p+1};
            case 'coordinates_mapping'
                if strcmpi(varargin{p+1},'projected')
                    parametric_coordinates=false;
                    projected_coordinates=true;
                elseif strcmpi(varargin{p+1},'parametric')
                    parametric_coordinates=true;
                    projected_coordinates=false;
                else
                    error('Unknown value for coordinates_mapping !')
                end
            case 'projected_resolution'
                projected_resolution=varargin{p+1};
            case 'subneb'
                subNEB=varargin{p+1};
            otherwise
                error('Unknown argument !')
        end
    end
end

cd(path)
% Read the title of the VASP calculation and find format of directory names (i.e. number of zero padding)
for p=1:10 % Format limitation: no more than 1e9 images (should be plenty enough)
    for q=0:9
        if exist(num2str(q,['%0',num2str(p),'i']),'dir')
            Title=grep('INCAR','SYSTEM','fwd',0);
            Title_idx=find(Title{1}=='=');
            EnergyPathway.Title=strip(Title{1}((Title_idx(1)+1):end));
            format_dir=['%0',num2str(p),'i'];
        end
    end
end

% Find all directories containing endpoints and images, and reading of the energies, positions and forces for each iteration and
% each image
ldir=dir('0*');
n_images=length(ldir)-2;
n_iter=[];
for p=1:n_images+2
    tmp_energies{p}=readEnergy([num2str(p-1,format_dir),'/']);
    n_iter=max([n_iter length(tmp_energies{p})]);
    if ~(p==1||p==n_images+2)&&length(tmp_energies{p})~=length(tmp_energies{2})
        warning('This calculation has not been successful, some images were still running when VASP stopped')
        halted_calculation=true;
    end
    tmp_Forces{p}=readForces([num2str(p-1,format_dir),'/']);
    tmp_MaxForces{p}=readMaxForces([num2str(p-1,format_dir),'/']);
    tmp_XDATCAR{p}=readPOSCAR([num2str(p-1,format_dir),'/XDATCAR']);
end
n_atoms=sum(tmp_XDATCAR{1}(1).n_chemicals);

% Sometimes NEB calculations stops (e.g. due to time limit) and usually all images did not reach the same number of iteration
if halted_calculation
    n_iter=n_iter-1;
end

for p=1:n_images+2
    POSCAR(p)=tmp_XDATCAR{p}(1);
    CONTCAR(p)=tmp_XDATCAR{p}(end);
    for q=1:n_iter
        if p==1||p==n_images+2
            XDATCAR(p,q)=tmp_XDATCAR{p}(1);
            energies(p,q)=tmp_energies{p}(1);
            Forces{p,q}=tmp_Forces{p}{1};
            MaxForces(p,q)=tmp_MaxForces{p}(1);
        else
            XDATCAR(p,q)=tmp_XDATCAR{p}(q);
            energies(p,q)=tmp_energies{p}(q);
            Forces{p,q}=tmp_Forces{p}{q};
            MaxForces(p,q)=tmp_MaxForces{p}(q);
        end
    end
end

% Recursive call to NEB_analysis to get subNEBs:
% All subNEB calculations should be put in a folder named "subNEB_" followed by the letter corresponding to the interval within
% which the subNEB was performed (e.g. "subNEB_C" if the subNEB ran between images in directories "02" and "03")
subdir=dir('subNEB_*');
subdir=subdir([subdir(:).isdir]);
if subNEB
    if ~isempty(subdir)
        for p=1:length(subdir)
            cd(subdir(p).name)
            if length(subdir(p).name)==8
                idx_up(p,:)=[double(subdir(p).name(end))-64 double(subdir(p).name(end))-63];
            elseif length(subdir(p).name)==9
                idx_up(p,:)=[double(subdir(p).name(end-1))-64 double(subdir(p).name(end))-63];
            end
            if save_data
                if parametric_coordinates
                    subrun(p)=NEB_analysis('save_data',filename,'subNEB',subNEB);
                elseif projected_coordinates
                    subrun(p)=NEB_analysis('save_data',filename,'coordinates_mapping','projected','projected_resolution',projected_resolution,'subNEB',subNEB);
                end
            else
                if parametric_coordinates
                    subrun(p)=NEB_analysis('subNEB',subNEB);
                elseif projected_coordinates
                    subrun(p)=NEB_analysis('coordinates_mapping','projected','projected_resolution',projected_resolution,'subNEB',subNEB);
                end
            end
            cd(curr_dir)
        end
    else
        subrun=NaN;
    end
end

x=coordinate_mapping(XDATCAR,parametric_coordinates,projected_coordinates,projected_resolution);

% Restructuring of the relevant data into a structure
EnergyPathway.POSCAR=POSCAR;
EnergyPathway.CONTCAR=CONTCAR;
EnergyPathway.XDATCAR=XDATCAR;
EnergyPathway.energies=energies;
EnergyPathway.reaction_coordinates=x;
EnergyPathway.Forces=Forces;
EnergyPathway.MaxForces=MaxForces;

if subNEB
    EnergyPathway.subrun=subrun;
    for p=1:length(subdir)
        EnergyPathway.subrun(p).idx_up=idx_up(p,:);
    end
end

clearvars -except EnergyPathway subNEB subdir projected_coordinates parametric_coordinates projected_resolution save_data curr_dir

% Reconstruction of the main reaction structure
k=0;
kp=0;
if subNEB
    while kp<length(EnergyPathway.POSCAR)
        kp=kp+1;
        k=k+1;
        POSCAR(k)=EnergyPathway.POSCAR(kp);
        CONTCAR(k)=EnergyPathway.CONTCAR(kp);
        XDATCAR(k,1)=EnergyPathway.XDATCAR(kp,end);
        energies(k,1)=EnergyPathway.energies(kp,end);
        Forces{k,1}=EnergyPathway.Forces{kp,end};
        MaxForces(k,1)=EnergyPathway.MaxForces(kp,end);
        for q=1:length(subdir)
            if EnergyPathway.subrun(q).idx_up(1)==kp
                for r=2:length(EnergyPathway.subrun(q).POSCAR)-1
                    k=k+1;
                    POSCAR(k)=EnergyPathway.subrun(q).POSCAR(r);
                    CONTCAR(k)=EnergyPathway.subrun(q).CONTCAR(r);
                    XDATCAR(k,1)=EnergyPathway.subrun(q).XDATCAR(r,end);
                    energies(k,1)=EnergyPathway.subrun(q).energies(r,end);
                    Forces{k,1}=EnergyPathway.subrun(q).Forces{r,end};
                    MaxForces(k,1)=EnergyPathway.subrun(q).MaxForces(r,end);
                end
                kp=kp+EnergyPathway.subrun(q).idx_up(2)-EnergyPathway.subrun(q).idx_up(1)-1;
            end
        end
    end
    EnergyPathway.POSCAR=POSCAR;
    EnergyPathway.CONTCAR=CONTCAR;
    EnergyPathway.XDATCAR=XDATCAR;
    EnergyPathway.energies=energies;
    EnergyPathway.reaction_coordinates=coordinate_mapping(XDATCAR,parametric_coordinates,projected_coordinates,projected_resolution);
    EnergyPathway.Forces=Forces;
    EnergyPathway.MaxForces=MaxForces;
    EnergyPathway=rmfield(EnergyPathway,'subrun');
    if isfield('idx_up',EnergyPathway)
        EnergyPathway=rmfield(EnergyPathway,'idx_up');
    end
end

% save if requested
if save_data
    save(filename,'EnergyPathway')
end
cd(curr_dir)

