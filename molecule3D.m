function axout = molecule3D(POSCAR,aLim,bLim,cLim,varargin)
%==================================================================================================================================%
% molecule3D.m: Draw 3D molecules (v2.1) (based on André Ludwig's - aludwig@phys.ethz.ch - 2020 molecule3D script)
%               (https://www.mathworks.com/matlabcentral/fileexchange/55231-molecule3d). Retrieved December 3, 2020.
%==================================================================================================================================%
% Version history:
%   version 1.2 (03/12/2020) - Creation
%       author: André Ludwig (aludwig@phys.ethz.ch)
%   version 1.3 (20/08/2025) - Added many bond length between elements, along with additional colors
%       author: EYG
%   version 2.0 (26/08/2025) - Complete overhaul of the code. Properties of atoms (colours, radii, bond length, ...) are now stored
%       author: EYG             in arrays and matrices which indices points to the atomic number of the elements. The input is now 
%                               a POSCAR structure. Many optional arguments have been implemented, such as the sphere resolution, 
%                               the choice between van der Waals or covalent radii, ... (see below). Visuals can be drawn directly 
%                               from this function instead of using the previous VESTA wrapper.
%   version 2.0.1 (01/09/2025) - The light position modification in the different viewing modes are now conditional to the use of
%       contrib: EYG                lights.
%   version 2.0.2 (04/09/2025) - Modification of the boundary duplication, previously too many cell were added in the positive
%       contrib: EYG                direction of all axes. In addition, the pair checking loops indices "p" and "q" only runs on the
%                                   Z that are present in the "Z" field of the POSCAR structure, instead of all atomic number 
%                                   (ranging from 1 to 118).
%   version 2.0.3 (24/09/2025) - Modification of the x/y/zlim in the case of non-periodical rendering IF there is no lattice vector
%       contrib: EYG                displayed: the limits are set to the extreme positions of the atoms in the POSCAR structure.
%   version 2.0.4 (29/11/2025) - Add bond length of Si, Ca, Pb and S and their colour according to the Jmol CPK convention. In
%       contrib: EYG                addition, the x/y/zlim are now set by the vertices of the lattices, or the extremes atoms
%                                   positions when requesting to display the lattice vectors.
%   version 2.1 (02-12-2025) - Add vdW bonding to the rendering, but only between O and H atoms, when the distance is above the OH
%       author: EYG            distance for covalent bonding and below the maximum van der Waals bond distance set to 2 Angstrom
%==================================================================================================================================%
% args:
%   POSCAR:             POSCAR structure, or path+filename of a POSCAR file
%   aLim, bLim, cLim:   Limits of the structure (in direct coordinates) along the a, b and c vectors
%   opt. args:          'save', followed by the format to which the 3D model must be saved (e.g. "png", "svg", ...)
%                           (default: nothing is saved)
%                       'periodic', followed by true or false to enable or disable periodic folding of the model into the lattice. 
%                           Note: turning the periodicity off allows to show atoms outside the unit cell, as they appear on the
%                           "positions" field of the POSCAR structure.
%                           (default: true)
%                       'rotate_view_angle', followed by a vector [theta,phi] to rotate the camera from its default direction
%                           (default: the model is shown from the spherical coordinates [19.5,9])
%                       'multi_view', followed by true or false to display additional camera direction in the normal plane to x, y
%                           and z.
%                           (default: false)
%                       'verbose', followed by "on" or "off" to display information about the rendering (timing and number of atoms)
%                           (default: "off")
%                       'lattice_vec', followed by 2D or 3D to display the limit of the cell as defined in the POSCAR structure
%                           Note: The 2D case implicitely assumes the system to be parallel to the xOy plane !
%                           (default: no lattice is shown)
%                       'bond_radius', followed by the diameter of the cylinder that represents the bond.
%                           Note: if the 'style' is set to 'liquorice', this tag also sets the radii of the atoms
%                           (default: RC=0.1)
%                       'smoothness', followed by the number of faces used to model the spheres and cylinders
%                           (default: NS=50, NB=50)
%                       'lights', followed by "on" or "off" to turn on/off the lights. 
%                           (default: "on")
%                       'style', followed by the style of representation of the atoms and bonds. Possible choices are:
%                           --> "ballstick" or "covalent": usual representation using balls and sticks, using the covalent radii 
%                               of the chemical elements.
%                           --> "van der Waals": use the van der Waals radii in a ball-and-stick representation
%                           --> "liquorice": modeling using only sticks, no balls shown
%                           --> "large": radii are twice as big as with the "ballstick" option 
%==================================================================================================================================%
% Parameters for the modeling : maximum interatomic distance, atoms colours and radii
%==================================================================================================================================%
% By default, if the interatomic distance is not set between a specific pair of atom, the bond will not be shown
DBOND=zeros(118,118)*NaN;
% Hydrogen
DBOND(1,1)=0.8;DBOND(1,2:118)=1.25;
% Lithium
DBOND(3,3)=2.5;DBOND(3,8)=2.3;DBOND(3,22)=2.5;
% Boron
DBOND(5,6)=1.7;DBOND(5,7)=1.6;
% Carbon
DBOND(6,6)=1.8;DBOND(6,7)=1.7;DBOND(6,8)=1.8;DBOND(6,15)=1.7;DBOND(6,17)=1.8;
% Oxygen
DBOND(8,8)=1.5;DBOND(8,14)=2.6;DBOND(8,16)=1.8;DBOND(8,20)=2.6;DBOND(8,22)=2.6;
% Sodium
DBOND(11,11)=4;DBOND(11,17)=3;
% Silicon
DBOND(14,14)=2.6;
% Sulfur
DBOND(16,42)=2.7;DBOND(16,82)=3.4;
% Titanium
DBOND(22,22)=2.5;
% Molybdenum
DBOND(42,42)=2.7;

% Set maximum bond distance 
DBOND_MAX=max(DBOND(:)); 
% Set maximum van der Waals bond distance between O and H
DBOND_vdW_OH_MAX=2;

% For practicality, make the distance matrix symetrical
for p=1:118
    for q=p:118
        DBOND(q,p)=DBOND(p,q);
    end
end

% Set colour for atoms and bonds (default for unknown atoms is purple-ish). The index of the line indicate the atomic number of the
% element.
atcol=ones(118,3).*[0.9 0.5 1.0];
atcol(1,:)=[0.95 0.95 0.95];
atcol(2,:)=[0.2 1.0 1.0];
atcol(3,:)=[0.5 0.1 1.0];
atcol(4,:)=[0.1 0.5 0.1];
atcol(5,:)=[58 163 45]/255;
atcol(6,:)=[0.3 0.3 0.3];
atcol(7,:)=[176 185 230]/255;
atcol(8,:)=[1.0 0.1 0.1];
atcol(9,:)=[0.2 0.9 0.2];
atcol(10,:)=atcol(2,:);
atcol(11,:)=atcol(3,:);
atcol(12,:)=atcol(4,:);
atcol(14,:)=[240 200 160]/255;
atcol(15,:)=[1.0 0.6 0.2];
atcol(16,:)=[0.9 0.9 0.2];
atcol(17,:)=atcol(9,:);
atcol(18,:)=atcol(2,:);
atcol(19,:)=atcol(3,:);
atcol(20,:)=[0.7 0.7 0.7];
atcol(22,:)=[0.6 0.6 0.6];
atcol(26,:)=[0.9 0.5 0.1];
atcol(35,:)=[0.6 0.1 0.1];
atcol(36,:)=atcol(2,:);
atcol(37,:)=atcol(3,:);
atcol(38,:)=atcol(4,:);
atcol(42,:)=[115 173 193]/255;
atcol(53,:)=[0.4 0.1 0.7];
atcol(54,:)=atcol(2,:);
atcol(55,:)=atcol(3,:);
atcol(56,:)=atcol(4,:);
atcol(82,:)=[87,89,97]/255;
atcol(88,:)=atcol(4,:);

% Set radius of atoms, considering either the van der Waals or covalent radius
covrad=ones(118,1)*0.5;
covrad(1)=0.2;
covrad(5)=0.3;
covrad(6)=0.3;
covrad(7)=0.4;
covrad(8)=0.35;
covrad(11)=0.25;
covrad(16)=0.3;
covrad(17)=0.4;
covrad(20)=0.5;
covrad(42)=0.4;
covrad(53)=0.6;
vdwrad=ones(118,1)*1.8;
vdwrad(1)=1.2;
vdwrad(6)=1.7;
vdwrad(8)=1.52;
vdwrad(16)=1.52;
%==================================================================================================================================%
% Handling of the mandatory parameters
%==================================================================================================================================%
% If the POSCAR input is a character string, the corresponding file is read into a POSCAR structure
if ischar(POSCAR)
    POSCAR=readPOSCAR(POSCAR);
end
% In the event an atom sits on the boundary defined by XLim, YLim and ZLim, it and its periodic copy are shown. This is possible if
% the windows of XLim, YLim and ZLim is slightly increased by a tolerance set by the tolpos variable.
tolpos=1e-3;
aLim=[aLim(1)-tolpos aLim(2)+tolpos];
bLim=[bLim(1)-tolpos bLim(2)+tolpos];
cLim=[cLim(1)-tolpos cLim(2)+tolpos];
%==================================================================================================================================%
% Initialisation of the default parameters and handling of the optional arguments
%==================================================================================================================================%
save=false;
periodic=true;
rotate_view_angle=[0 0];
default_view_angle=[19.5 9];
multi_view=false;
verbose=false;
lattice_vec3D=false;
lattice_vec2D=false;
RC = 0.1; % Bond cylinder radius
NS = 50; % Number of faces on the spheres
NB = 50; % Number of faces on the cylinders
lights=true;
style='ballstick';
% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'save'
                save=true;
                save_format=varargin{p+1};
            case 'periodic'
                if ~varargin{p+1}
                    periodic=false;
                end
            case 'rotate_view_angle'
                default_view_angle=[0 0];
                rotate_view_angle=varargin{p+1};
            case 'multi_view'
                if varargin{p+1}
                    multi_view=true;
                end
            case 'verbose'
                if strcmpi(varargin{p+1},'on')
                    verbose=true;
                    tic
                elseif strcmpi(varargin{p+1},'off')
                    verbose=false;
                else
                    warning('unknown argument for "verbose" parameter (choices are "on" or "off") --> reset to default ("off")')
                end
            case 'lattice_vec'
                if ~isnan(varargin{p+1})
                    if strcmpi(varargin{p+1},'3D')
                        lattice_vec3D=true;
                    elseif strcmpi(varargin{p+1},'2D')
                        lattice_vec2D=true;
                    else
                        error('unknown argument for "lattice_vec" parameter (choices are "2D" or "3D")')
                    end
                end
            case 'bond_radius'
                RC=varargin{p+1};
            case 'smoothness'
                NS=varargin{p+1};
                NB=varargin{p+1};
            case 'lights'
                if strcmpi(varargin{p+1},'on')
                    lights=true;
                elseif strcmpi(varargin{p+1},'off')
                    lights=false;
                else
                    error('unknown argument for "lights" parameter (choices are "on" or "off")')
                end
            case 'style'
                style=varargin{p+1};
            otherwise 
                error('unknown optional argument')
        end
    end
end
%==================================================================================================================================%
% Computation of the position of the atoms within the boundary set by XLim, YLim, ZLim
%==================================================================================================================================%
% Duplication of the positions of the atoms in the neighboring periodical copies of the unit cell (with one unit cell outside of
% the boundary set by XLim, YLim and ZLim). Doing so "brings back" into the unit cell all the atoms.
POSCAR_init=POSCAR;
if periodic
    for p=floor(aLim(1)):ceil(aLim(2))-1
        for q=floor(bLim(1)):ceil(bLim(2))-1
            for r=floor(cLim(1)):ceil(cLim(2))-1
                if ~(p==0&&q==0&&r==0)
                    POSCAR_tmp=POSCAR_init;
                    POSCAR_tmp.positions=POSCAR_tmp.positions+p*POSCAR.vec(1,:)+q*POSCAR.vec(2,:)+r*POSCAR.vec(3,:);
                    POSCAR_tmp.xred=POSCAR.xred+p+q+r;
                    POSCAR=appendPOSCAR(POSCAR,POSCAR_tmp);
                end
            end
        end
    end
end
% Deletion of all the atoms that are beyond the boundaries set by XLim, YLim and ZLim
count=0;
max_count=sum(POSCAR.n_chemicals);
k=1;
while count<max_count
    if POSCAR.xred(k,1)<aLim(1)||POSCAR.xred(k,1)>aLim(2)
        POSCAR=delPOSCAR(POSCAR,k);
        k=k-1;
    elseif POSCAR.xred(k,2)<bLim(1)||POSCAR.xred(k,2)>bLim(2)
        POSCAR=delPOSCAR(POSCAR,k);
        k=k-1;
    elseif POSCAR.xred(k,3)<cLim(1)||POSCAR.xred(k,3)>cLim(2)
        POSCAR=delPOSCAR(POSCAR,k);
        k=k-1;
    end
    k=k+1;
    count=count+1;
end
% Display of the number of atoms that are in the POSCAR structure and the number of atoms actually shown
if verbose
    disp(['Rendering ',num2str(sum(POSCAR.n_chemicals)),' atoms from a ',num2str(sum(POSCAR_init.n_chemicals)),' atoms system'])
end
%==================================================================================================================================%
% Detect atom pairs to be displayed
%==================================================================================================================================%
% Extract useful info from POSCAR structure for efficiency
xyz=POSCAR.positions;
Z=POSCAR.Z;
Natoms=sum(POSCAR.n_chemicals);

% There is no point in showing a single (or no) atom
if Natoms<=1
    error('There is only one or zero atom in the rendering !')
end

% get all combinations of atoms
pairs = nchoosek(1:Natoms,2);
if isempty(pairs)
    warning('No atom pair detected !')
end

% get all distance between pairs, then restrict the list to the distance below a threshold (here the maximum distance in the DBOND
% matrix). This is only a first clean-up of all the possible pairs.
ds=sqrt(sum((xyz(pairs(:,1),:)-xyz(pairs(:,2),:)).^2,2)); 
ks=find(ds<DBOND_MAX)';
% Browsing through the list, check the symbols and keep the pair only if the distance is lower than the one defined in the DBOND
% matrix.
k=1;
kv=1;
for p_ks=1:length(ks)
    for pZ=unique(POSCAR.Z)
        for qZ=unique(POSCAR.Z)
            if Z(pairs(ks(p_ks),1))==pZ&&Z(pairs(ks(p_ks),2))==qZ
                if ds(ks(p_ks))<DBOND(pZ,qZ)
                    ks_n(k)=ks(p_ks);
                    k=k+1;
                elseif ((pZ==1&&qZ==8)||(pZ==8&&qZ==1))&&ds(ks(p_ks))<DBOND_vdW_OH_MAX
                    ksv_n(kv)=ks(p_ks);
                    kv=kv+1;
                end
            end
        end
    end
end
if isempty(pairs)
    warning('No atom pair detected below the DBOND thresholds ! Check lattice parameter !')
else
    ks=ks_n;
    if exist('ksv_n','var')
        kv=ksv_n;
    end
end
%==================================================================================================================================%
% Drawing of the spheres and cylinders
%==================================================================================================================================%
ax=[];
ax=newplot(ax);
axes(ax);
set(ax,'Visible','off','DataAspectRatio',[1 1 1])
% draw sphere with adapted radius for each atom of the POSCAR structure
for k = 1:size(xyz,1)
    switch lower(style)
        case {'ballstick','covalent'}
            radii(k)=covrad(POSCAR.Z(k));
            offset=0.5;
        case 'licorice'
            radii(k)=RC;
            offset=0.2;
        case 'vdw'
            radii(k)=vdwrad(POSCAR.Z(k));
            offset=2;
        case 'large'
            radii(k)=2*covrad(POSCAR.Z(k));
            offset=1;
        otherwise 
            warning('unknown style --> ballstick will be used')
    end
    [sx,sy,sz]=sphere(NS);
    sx=xyz(k,1)+radii(k)*sx;
    sy=xyz(k,2)+radii(k)*sy;
    sz=xyz(k,3)+radii(k)*sz;
    surface(sx,sy,sz,'FaceColor',atcol(POSCAR.Z(k),:),'EdgeColor','none','FaceLighting','gouraud')
end

% draw cylinders for the pair which index k is included in the restricted list ks
for p_ks=ks
    pos(1,:)=xyz(pairs(p_ks,1),:); % coordinates atom 1
    pos(2,:)=xyz(pairs(p_ks,2),:); % coordinates atom 2
    % bond vector in spherical coordinates
    n=(pos(2,:)-pos(1,:))/norm(pos(2,:)-pos(1,:));
    phi=atan2(n(2),n(1));
    theta=-asin(n(3));
    % bond distance minus sphere radii to get the distance between the surface of the spheres
    bd=ds(p_ks)-radii(pairs(p_ks,1))-radii(pairs(p_ks,2));
    delta_s(1)=radii(pairs(p_ks,1))+bd/2;
    delta_s(2)=radii(pairs(p_ks,2))+bd/2;
    for p=1:2
        % prototype (unitary length) cylinders for bond aligned along the x-axis
        [cylz,cyly,cylx]=cylinder(RC,NB);
        % first half-bond rescaling
        cylx(2,:)=cylx(2,:)*delta_s(p);
        % rotation of the cylinders to match bond vector n and their translation at the location of their respective anchor atom
        for k=1:length(cylx(:))
            cylr=rotz(phi)*roty(theta)*[cylx(k);cyly(k);cylz(k)];
            cylx(k)=cylr(1);
            cyly(k)=cylr(2);
            cylz(k)=cylr(3);
        end
        cylx=pos(p,1)+(-1)^(p+1)*cylx;
        cyly=pos(p,2)+(-1)^(p+1)*cyly;
        cylz=pos(p,3)+(-1)^(p+1)*cylz;
        surface(cylx,cyly,cylz,'FaceColor',atcol(POSCAR.Z(pairs(p_ks,p)),:),'EdgeColor','none','FaceLighting','gouraud')
    end
end
if exist('ksv_n','var')
    for p_ks=kv
        pos(1,:)=xyz(pairs(p_ks,1),:); % coordinates atom 1
        pos(2,:)=xyz(pairs(p_ks,2),:); % coordinates atom 2
        line([pos(1,1) pos(2,1)],[pos(1,2) pos(2,2)],[pos(1,3) pos(2,3)],'LineStyle','--','color',[0.5 0.5 0.5])
    end
end
%==================================================================================================================================%
% Polishing and finishing
%==================================================================================================================================%
% convert aLim, bLim and cLim to XLim, YLim and ZLim
XLim(1)=[aLim(1) bLim(1) cLim(1)]*POSCAR(1).vec(:,1);
YLim(1)=[aLim(1) bLim(1) cLim(1)]*POSCAR(1).vec(:,2);
ZLim(1)=[aLim(1) bLim(1) cLim(1)]*POSCAR(1).vec(:,3);
XLim(2)=[aLim(2) bLim(2) cLim(2)]*POSCAR(1).vec(:,1);
YLim(2)=[aLim(2) bLim(2) cLim(2)]*POSCAR(1).vec(:,2);
ZLim(2)=[aLim(2) bLim(2) cLim(2)]*POSCAR(1).vec(:,3);
% There is an offset of 0.2 to 1 Å (depending on the rendering style) in all direction to display the full sphere, even when the 
% atom is on the boundary
xlim([min([min(xyz(:,1)) XLim(1)]) max([max(xyz(:,1)) XLim(2)])]+[-1 1]*offset)
ylim([min([min(xyz(:,2)) YLim(1)]) max([max(xyz(:,2)) YLim(2)])]+[-1 1]*offset)
zlim([min([min(xyz(:,3)) ZLim(1)]) max([max(xyz(:,3)) ZLim(2)])]+[-1 1]*offset)
if ~periodic&&~lattice_vec3D&&~lattice_vec2D
    xlim([min(xyz(:,1)) max(xyz(:,1))]+[-1 1]*offset)
    ylim([min(xyz(:,2)) max(xyz(:,2))]+[-1 1]*offset)
    zlim([min(xyz(:,3)) max(xyz(:,3))]+[-1 1]*offset)
end
% Draw the unit cell if the lattice_vec parameter is set.
if lattice_vec3D
    a=POSCAR.vec(1,:);
    b=POSCAR.vec(2,:);
    c=POSCAR.vec(3,:);
    lattice_vertices_x=[0 a(1);0 b(1);0 c(1);a(1) a(1)+b(1);a(1) a(1)+c(1);b(1) a(1)+b(1);b(1) b(1)+c(1);c(1) a(1)+c(1);c(1) b(1)+c(1);a(1)+b(1) a(1)+b(1)+c(1);a(1)+c(1) a(1)+b(1)+c(1);b(1)+c(1) a(1)+b(1)+c(1)];
    lattice_vertices_y=[0 a(2);0 b(2);0 c(2);a(2) a(2)+b(2);a(2) a(2)+c(2);b(2) a(2)+b(2);b(2) b(2)+c(2);c(2) a(2)+c(2);c(2) b(2)+c(2);a(2)+b(2) a(2)+b(2)+c(2);a(2)+c(2) a(2)+b(2)+c(2);b(2)+c(2) a(2)+b(2)+c(2)];
    lattice_vertices_z=[0 a(3);0 b(3);0 c(3);a(3) a(3)+b(3);a(3) a(3)+c(3);b(3) a(3)+b(3);b(3) b(3)+c(3);c(3) a(3)+c(3);c(3) b(3)+c(3);a(3)+b(3) a(3)+b(3)+c(3);a(3)+c(3) a(3)+b(3)+c(3);b(3)+c(3) a(3)+b(3)+c(3)];	
    line(lattice_vertices_x(1,:),lattice_vertices_y(1,:),lattice_vertices_z(1,:),'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5)
    line(lattice_vertices_x(2,:),lattice_vertices_y(2,:),lattice_vertices_z(2,:),'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5)
    line(lattice_vertices_x(3,:),lattice_vertices_y(3,:),lattice_vertices_z(3,:),'color',[0.3010, 0.7450, 0.9330],'LineWidth',1.5)
    for p=4:12
        line(lattice_vertices_x(p,:),lattice_vertices_y(p,:),lattice_vertices_z(p,:),'color',[1 1 1]*0.25,'LineWidth',1,'LineStyle','--')
    end
    xl=xlim;
    xlim([min([lattice_vertices_x(:);xl(1)]) max([lattice_vertices_x(:);xl(2)])])
    yl=ylim;
    ylim([min([lattice_vertices_y(:);yl(1)]) max([lattice_vertices_y(:);yl(2)])])
    zl=zlim;
    zlim([min([lattice_vertices_z(:);zl(1)]) max([lattice_vertices_z(:);zl(2)])])
elseif lattice_vec2D
    a=POSCAR.vec(1,:);
    b=POSCAR.vec(2,:);
    avg_h=mean(POSCAR.positions(:,3));
    line([   0      a(1)],[   0      a(2)],[avg_h avg_h],'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5)
    line([   0      b(1)],[   0      b(2)],[avg_h avg_h],'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5)
    line([a(1) a(1)+b(1)],[a(2) a(2)+b(2)],[avg_h avg_h],'color',[1 1 1]*0.25,'LineWidth',1,'LineStyle','--')
    line([b(1) a(1)+b(1)],[b(2) a(2)+b(2)],[avg_h avg_h],'color',[1 1 1]*0.25,'LineWidth',1,'LineStyle','--')
end
% Add lights
if lights
    cl=camlight(0,0);
end
% Set camera direction (and optionally save)
if ~multi_view
    view(default_view_angle+rotate_view_angle)
    if save
        saveas(gcf,['standard.',save_format])
    end
    if lights
        camlight(cl,"left");
    end
elseif save
    view([0,0]+rotate_view_angle)
    if lights
        camlight(cl,"left");
    end
    saveas(gcf,['a_view.',save_format])
    view([90,0]+rotate_view_angle)
    if lights
        camlight(cl,"left");
    end
    saveas(gcf,['b_view.',save_format])
    view([0,90]+rotate_view_angle)
    if lights
        camlight(cl,"left");
    end
    saveas(gcf,['c_view.',save_format])
    view(default_view_angle+rotate_view_angle)
    if lights
        camlight(cl,"left");
    end
    saveas(gcf,['standard.',save_format])
else
    view([0,0]+rotate_view_angle)
    if lights
        camlight(cl,"left");
    end
    pause(1)
    view([90,0]+rotate_view_angle)
    if lights
        camlight(cl,"left");
    end
    pause(1)
    view([0,90]+rotate_view_angle)
    if lights
        camlight(cl,"left");
    end
    pause(1)
    view(default_view_angle+rotate_view_angle)
    if lights
        camlight(cl,"left");
    end
end
% get timing of rendering
if verbose
    toc
end
% if output arguments specified in the call to the currently executing function, outputs the axes
if nargout > 0
    axout = ax;
end