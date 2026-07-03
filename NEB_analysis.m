function EnergyPathway=NEB_analysis(varargin)
%==================================================================================================================================%
% NEB_analysis.m:   Post-processing of a (c)NEB calculation (v0.3)
%==================================================================================================================================%
% Version history:
%   version 0.1 (02/09/2025) - Creation using bits and pieces from my thesis works
%       author: EYG
%   version 0.1.1 (16/12/2025) - Correction of the call to readPOSCAR to account for the removal of optional argument 'XDATCAR'
%       contrib: EYG
%   version 0.1.2 (03/01/2026) - The recursive call only occurs on folders that starts with "subNEB" instead of "sub".
%       contrib: EYG
%   version 0.1.3 (26/01/2026) - Sanity check has been included to ensure continuity of the DFT settings across the NEB calculation.
%       contrib: EYG                Current parameters include NELECT, NUPDOWN, NKPTS and IVDW.
%   version 0.2 (16/02/2026) - Addition of sanity check over the frozen DOF, plus on the actual magnetisation of the system. There
%       author: EYG             is now a separate visualisation of the results for single-iteration calculations. In addition,
%                               titles of POSCAR structures are chosen to reflect the number of the images and the subNEB
%                               calculation it originated from when relevant.
%   version 0.2.1 (18/02/2026) - Deletion of the field "str_xticklab" from the EnergyPathway structure and correction of a typo in  
%       contrib: EYG                a warning.
%   version 0.2.2 (02/03/2026) - Reading of files now occurs in the folders in the "ldir" list.
%       contrib: EYG
%   version 0.2.3 (12/06/2026) - Fix of the sanity check to ensure that all atoms moving during the reaction are free to move in 
%       contrib: EYG                each image.
%   version 0.3 (01/07/2026) - This function is no longer recursive, the list of directories is generated recursively with the 
%       author: EYG             "subdir_ordering" function. In addition, the images have now a depth (i.e. the number of subNEB
%                               that were used to find them). This allows for an explicit discrimination between parent and child
%                               images when plotting the energy profile.
%   version 0.3.1 (03/07/2026) - Minor fix, the subdir_ordering function was called with a target to 'path', but after the cd 
%       contrib: EYG                command. Now it points to "./".
%==================================================================================================================================%
% args:
%   opt. args:          'path', followed by the path to the NEB directory
%                           (default: current directory)
%                       'save_data', followed by the filename where the NEB structure must be saved
%                           (default: nothing is saved)
%                       'verbose', followed by true of false to specify whether text should be shown to indicate NEB parameters
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
visual_only=false;
save_data=false;
parametric_coordinates=true;
projected_coordinates=false;
projected_resolution=0.001;
subNEB=false;
idx_up=[NaN NaN];
filename='EP.mat';
curr_dir=[pwd,'/'];
verbose=false;
postprocessing=false;
halted_calculation=false;
path='./';
% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        if isstruct(varargin{1})
            visual_only=true;
            postprocessing=true;
            EnergyPathway=varargin{1};
        else
            switch lower(varargin{p})
                case 'path'
                    path=varargin{p+1};
                    if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
                        path=[path,'/'];
                    end
                case 'save_data'
                    save_data=true;
                    filename=varargin{p+1};
                case 'visual'
                    visual_only=true;
                    postprocessing=true;
                    EnergyPathway=varargin{p+1};
                case 'verbose'
                    verbose=varargin{p+1};
                case 'postprocessing'
                    postprocessing=varargin{p+1};
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
end
if ~visual_only
    cd(path)
    % Read the title of the VASP calculation and find format of directory names (i.e. number of zero padding)
    for p=1:10 % Format limitation: no more than 1e9 images (should be plenty enough)
        for q=0:9
            if exist(num2str(q,['%0',num2str(p),'i']),'dir')
                if exist('INCAR','file')
                    Title=grep('INCAR','SYSTEM','fwd',0);
                    Title_idx=find(Title{1}=='=');
                    EnergyPathway.Title=strip(Title{1}((Title_idx(1)+1):end));
                else
                    EnergyPathway.Title='NEB';
                end
                format_dir=['%0',num2str(p),'i'];
            end
        end
    end
    % Find all directories containing endpoints and images, and reading of the energies, positions and forces for each iteration and
    % each image. This happen using the "subdir_ordering" function when subNEB is enabled.
    if subNEB
        [ldir,~,lsname]=subdir_ordering('./');
    else
        ldir=dir('0*');
        for p=1:length(ldir)
            image_depth{p}=0;
            lsname{p}=ldir(p).name;
        end
        [ldir.image_depth]=image_depth{:};
        clear image_depth
    end
    Title=grep([ldir(1).name,'/OUTCAR'],'SYSTEM','fwd',0);
    Title_idx=find(Title{1}=='=');
    EnergyPathway.Title=strip(Title{1}((Title_idx(1)+1):end));
    n_images=length(ldir)-2;
    n_iter=[];
    for p=1:n_images+2
        tmp_energies{p}=readEnergy(ldir(p).name);
        n_iter=max([n_iter length(tmp_energies{p})]);
        if ~(p==1||p==n_images+2)&&length(tmp_energies{p})~=length(tmp_energies{2})
            warning('This calculation has not been successful, some images were still running when VASP stopped')
            halted_calculation=true;
        end
        tmp_Forces{p}=readForces(ldir(p).name);
        tmp_MaxForces{p}=readMaxForces(ldir(p).name);
        tmp_XDATCAR{p}=readPOSCAR([ldir(p).name,'/XDATCAR']);
        tmp_POSCAR(p)=readPOSCAR([ldir(p).name,'/POSCAR']);
        tmp_CONTCAR(p)=readPOSCAR([ldir(p).name,'/CONTCAR']);
    end
    n_atoms=sum(tmp_XDATCAR{1}(1).n_chemicals);
    
    % Sanity check: ensure that the most important parameters are the same across the NEB calculation. If verbose is activated, shows
    % the value of these parameters (currently, the NKPTS, IVDW and NUPDOWN settings are shown)
    for p=1:n_images+2
        grep_tmp=grep([ldir(p).name,'/OUTCAR'],'NKPTS');
        nkpts(p)=str2num(grep_tmp{1}(31:36));
        grep_tmp=grep([ldir(p).name,'/OUTCAR'],'IVDW');
        ivdw(p)=str2num(grep_tmp{1}(end-1:end));
        grep_tmp=grep([ldir(p).name,'/OUTCAR'],'NUPDOWN');
        nupdown(p)=str2num(grep_tmp{end}(13:24));
        grep_tmp=grep([ldir(p).name,'/OUTCAR'],'NELECT');
        nelect(p)=str2num(grep_tmp{1}(13:24));
        if p~=1&&p~=n_images+2
            grep_tmp=grep([ldir(p).name,'/OUTCAR'],'stopping-criterion for IOM');
            EDIFFG=str2num(grep_tmp{1}(13:19));
        end
    end
    if any(nkpts~=nkpts(1))
        warning('Severe problem found: not all images were computed with the same BZ sampling!')
    elseif nkpts(1)==1&&verbose
        disp('Gamma only calculation')
    elseif verbose
        disp([num2str(nkpts(1)),' k-points are considered in the IBZ.'])
    end
    if any(ivdw~=ivdw(1))
        warning('Severe problem found: not all images were computed with the same vdW correction scheme!')
    elseif verbose
        disp(['van der Waals correction scheme: ',num2str(ivdw(1))])
    end
    if any(nupdown~=nupdown(1))
        warning('NUPDOWN/MAGMON constraint are not identical for all images!')
        for p=1:n_images+2
            grep_tmp=grep([ldir(p).name,'/OUTCAR'],'number of electron ');
            effective_nupdown(p)=str2num(grep_tmp{end}(55:65));
        end
        if any(effective_nupdown~=effective_nupdown(1))
            warning('Not all images were computed with the same number of up and down electrons!')
        end
    elseif nupdown(1)==-1&&verbose
        disp('No magnetisation enforced')
    elseif verbose
        disp(['Total magnetisation forced to ',num2str(nupdown(1))])
    end
    if any(nelect~=nelect(1))
        error('Severe problem found: not all images contains the same number of electrons!')
    elseif nupdown(1)~=-1&&(mod(nelect(1),2)==0&&mod(nupdown(1),2)==1)
        error('Wrong magnetisation enforced: even number of electrons, but odd difference of up and down electrons!')
    elseif nupdown(1)~=-1&&(mod(nelect(1),2)==1&&mod(nupdown(1),2)==0)
        error('Wrong magnetisation enforced: odd number of electrons, but even difference of up and down electrons!')
    end
    
    % Sometimes NEB calculations stops (e.g. due to time limit) and usually all images did not reach the same number of iteration
    if halted_calculation
        n_iter=n_iter-1;
    end
    constraint=true(n_atoms,3);
    for p=1:n_images+2
        POSCAR(p)=tmp_POSCAR(p);
        CONTCAR(p)=tmp_CONTCAR(p);
        for q=1:n_iter
            if p==1||p==n_images+2
                XDATCAR(p,q)=tmp_XDATCAR{p}(1);
                energies(p,q)=tmp_energies{p}(1);
                Forces{p,q}=tmp_Forces{p}{1};
                MaxForces(p,q)=tmp_MaxForces{p}(1);
            else
                constraint=constraint.*POSCAR(p).constraint;
                XDATCAR(p,q)=tmp_XDATCAR{p}(q);
                energies(p,q)=tmp_energies{p}(q);
                Forces{p,q}=tmp_Forces{p}{q};
                MaxForces(p,q)=tmp_MaxForces{p}(q);
            end
        end
        tmp_str_xticklab{p}=num2str(p-1,format_dir);
    end
        
    % Additional sanity check on the constraint over the degree of freedom of moving atoms:
    if parametric_coordinates
        [x,ds_vec]=coordinate_mapping(XDATCAR,parametric_coordinates,projected_coordinates,projected_resolution);
        ds_tot_atoms=zeros(n_atoms,1);
        for p=1:length(ds_vec)
            ds_tot_atoms=ds_tot_atoms+vecnorm(ds_vec{p}')';
        end
        for p=1:n_atoms
            if any(~constraint(p,:))&&ds_tot_atoms(p)~=0
                warning(sprintf(['Severe problem found: position of (moving) atom ',num2str(p),' has been frozen along at least one direction! Yet the atom is displaced by ',num2str(ds_tot_atoms(p),'%6.2e'),' Angstrom during the reaction!']))
            end
        end
    else
        [x,~]=coordinate_mapping(XDATCAR,parametric_coordinates,projected_coordinates,projected_resolution);
    end

    % Restructuring of the relevant data into a structure
    EnergyPathway.POSCAR=POSCAR;
    EnergyPathway.CONTCAR=CONTCAR;
    EnergyPathway.XDATCAR=XDATCAR;
    EnergyPathway.energies=energies;
    EnergyPathway.reaction_coordinates=x;
    EnergyPathway.Forces=Forces;
    EnergyPathway.MaxForces=MaxForces;
    for p=1:length(ldir)
        image_depth(p)=ldir(p).image_depth;
    end
    EnergyPathway.ImDepth=image_depth;
    
    % Reconstruction of the main reaction structure
    k=0;
    kp=0;

    if exist('EDIFFG','var')
        EnergyPathway.EDIFFG=EDIFFG;
    else
        EnergyPathway.EDIFFG=NaN;
    end
    % save if requested
    if save_data
        if isfield('str_xticklab',EnergyPathway)
            EnergyPathway=rmfield(EnergyPathway,'str_xticklab');
        end
        save(filename,'EnergyPathway')
    end
    cd(curr_dir)
end
subNEB=false;
cmap=colororder;

% Generate labels for each image
for p=1:length(lsname)
    xlsname{p}=['x(',lsname{p},')'];
end

% To show the scale, the endpoint of the parent NEB are shown as equal to 0 and 1 respectively
xlsname{1}=[xlsname{1},'=0'];
xlsname{end}=[xlsname{end},'=1'];
if postprocessing&&~subNEB&&length(EnergyPathway.XDATCAR(1,:))~=1
    % Figure 1: Evolution of the energy and reaction coordinate throughout the calculation
    figure;
    plot(EnergyPathway.reaction_coordinates(2:end-1,:)',EnergyPathway.energies(2:end-1,:)'-EnergyPathway.energies(1,end))
    hold on
    grid on
    for p=1:length(EnergyPathway.POSCAR)-2
        plot(EnergyPathway.reaction_coordinates(p+1,end)',EnergyPathway.energies(p+1,end)'-EnergyPathway.energies(1,end),'o','color',cmap(p,:))
    end
    xlim([0 1])
    xlabel('Reaction coordinate')
    ylabel('Energy')
    for p=1:length(EnergyPathway.POSCAR)-2
        str_legend{p}=['Image ',num2str(p)];
    end
    legend(str_legend,'Location','best')
    set(gca,'fontsize',12,'fontname','cambria math')

    % Figure 2: Maximum force acting on the atoms of the systems for each ionic iteration
    figure;
    semilogy(EnergyPathway.MaxForces(2:end-1,:)')
    hold on
    grid on
    h=semilogy([1 length(EnergyPathway.MaxForces(2,:)')],abs([EnergyPathway.EDIFFG EnergyPathway.EDIFFG]),'--','color',[1 1 1]*0.5);
    xlim([1 length(EnergyPathway.MaxForces(2,:)')])
    xlabel('Number of ionic iteration')
    ylabel('Maximum force acting on atoms')
    for p=1:length(EnergyPathway.POSCAR)-2
        str_legend{p}=['Image ',num2str(p)];
    end
    legend(str_legend,'Location','best')
    set(gca,'fontsize',12,'fontname','cambria math')
    uistack(h, 'bottom')

    % Figure 3: Energy difference between consecutive ionic iterations
    figure;
    plot(EnergyPathway.energies(2:end-1,2:end)'-EnergyPathway.energies(2:end-1,1:end-1)')
    hold on
    grid on
    xlim([1 length(EnergyPathway.MaxForces(2,:)')])
    xlabel('Number of ionic iteration')
    ylabel('\Delta E (eV)')
    for p=1:length(EnergyPathway.POSCAR)-2
        str_legend{p}=['Image ',num2str(p)];
    end
    legend(str_legend,'Location','best')
    set(gca,'fontsize',12,'fontname','cambria math')
elseif postprocessing&&~subNEB&&isscalar(EnergyPathway.XDATCAR(1,:)) % Case of a single iteration calculation
    % Figure 1: Energy along the reaction coordinate. Parent images are shown in red, their child in yellow, their own child
    % in cmap(EnergyPathway.ImDepth(1)+2,:)), etc ...
    figure;
    plot(EnergyPathway.reaction_coordinates',EnergyPathway.energies'-EnergyPathway.energies(1),'-','LineWidth',4)
    hold on
    for p=1:length(ldir)
        plot(EnergyPathway.reaction_coordinates(p,end)',EnergyPathway.energies(p,end)'-EnergyPathway.energies(1),'.','MarkerFaceColor',cmap(EnergyPathway.ImDepth(p)+2,:),'MarkerEdgeColor',cmap(EnergyPathway.ImDepth(p)+2,:),'MarkerSize',32)
    end
    xlim([0 1])
    grid on
    xlabel('Reaction coordinate')
    clear str_xticklab
    for p=1:length(EnergyPathway.POSCAR)
        str_xticklab{p}=EnergyPathway.POSCAR(p).Title;
    end
    hold on
    grid on
    xlabel('Reaction coordinate x')
    xticks(EnergyPathway.reaction_coordinates')
    xticklabels(xlsname)
    ylabel('Energy')
    set(gca,'fontsize',12,'fontname','cambria math')

    % Figure 2: Max force acting on the atoms of the system. Parent images are shown in red, their child in yellow, their own child
    % in cmap(EnergyPathway.ImDepth(1)+2,:)), etc ...
    figure;
    h=bar(0,EnergyPathway.MaxForces(1)','FaceColor',cmap(EnergyPathway.ImDepth(1)+2,:),'EdgeColor',cmap(EnergyPathway.ImDepth(1)+2,:));
    hold on
    for p=2:length(ldir)
        bar(p-1,EnergyPathway.MaxForces(p)','FaceColor',cmap(EnergyPathway.ImDepth(p)+2,:),'EdgeColor',cmap(EnergyPathway.ImDepth(p)+2,:));
    end
    set(get(h,'Parent'), 'YScale', 'log')
    ylim([min(EnergyPathway.MaxForces)/2 max(EnergyPathway.MaxForces)*2])
    clear str_xticklab
    for p=1:length(EnergyPathway.POSCAR)
        str_xticklab{p}=EnergyPathway.POSCAR(p).Title;
    end
    hold on
    grid on
    xlabel('Image number')
    xticks(0:length(EnergyPathway.MaxForces)-1)
    xticklabels(lsname)
    ylabel('Maximum force acting on atoms')
    set(gca,'fontsize',12,'fontname','cambria math')
end
if isfield('str_xticklab',EnergyPathway)
    EnergyPathway=rmfield(EnergyPathway,'str_xticklab');
end
cd(curr_dir)