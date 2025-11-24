function frames=ReactionRendering(R_in,aLim,bLim,cLim,varargin)
%==================================================================================================================================%
% ReactionRendering.m:  Render an animated or static 3D model of an atomic structure, with or without the corresponding energy 
%                       profile (v0.3.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (28/08/2025) - Creation
%       author: EYG
%   version 0.1.1 (01/09/2025) - The call to molecule3D.m now happens wo the lights on, to be turned on in this script.
%       contrib: EYG
%   version 0.2 (02/09/2025) - The POSCAR structure input has been replaced by a reaction structure input. However, if there is no
%       author: EYG             POSCAR field in the input structure, the script assumes an array of POSCAR structures. In addition,
%                               the absolute limits of the rendering have been redefined to account for non-orthorombic cells.
%   version 0.3 (03/09/2025) - Add the possibility to render the lattice vector
%       author: EYG
%   version 0.3.1 (10/09/2025) - Add an auto-updating information line information as a sort of progression bar, and also the
%                                   background color of the figure is set to white
%==================================================================================================================================%
% args:
%   R_in:               Array of reaction structures (or POSCAR structures)
%   aLim, bLim, cLim:   Limits of the structure (in direct coordinates) along the a, b and c vectors
%   opt. args:          'save_data', followed by the playing direction of the GIF. Possible choices are 'Fwd', 'Bwd', 'Round', 'Infinite' 
%                           (default: nothing is saved)
%                       'static', followed by true or false to request a static rendering.
%                           (default: false)
%                       'rotate_view_angle', followed by a vector [theta,phi] to rotate the camera from its default direction
%                           (default: the model is shown from the spherical coordinates [19.5,9])
%                       Note: in that case, the reaction structure at the last position of the array must concern a single 
%                               configuration.
%                       'lattice_vec', followed by either "2D" or "3D" to specify the type of lattice vector to be displayed
%==================================================================================================================================%
% set default parameters
EnergyProfile=true;
save_data=false;
gif_direction='infinite';
lcolor=colororder;
static=false;
rotate_view_angle=[0 0];
default_view_angle=[19.5 9];
lattice_vec=NaN;
% Reading of the optional argument
if exist('varargin','var')
    for p=1:2:length(varargin)
        switch lower(varargin{p})
            case 'save_data'
                save_data=true;
                gif_direction=varargin{p+1};
            case 'static'
                static=true;
                idxTS=varargin{p+1};
            case 'rotate_view_angle'
                default_view_angle=[0 0];
                rotate_view_angle=varargin{p+1};
            case 'lattice_vec'
                lattice_vec=varargin{p+1};
        end
    end
end
% check for consistency between x and E
if ~isfield(R_in,'POSCAR')
    POSCARs=R_in;
    EnergyProfile=false;
    n_images=length(POSCARs);
    figure('color','white')
else
    POSCARs=R_in(end).CONTCAR;
    E0=R_in(1).energies(1,end);
    n_images=length(R_in(end).POSCAR);
    n_dataset=length(R_in)-1;
    % plot the static part of the energy profile only once at the beginning
    figure('Position', [200 120 1200 425],'color','white')
    subplot(1,2,2);
    plot(R_in(1).reaction_coordinates(:,end),R_in(1).energies(:,end)-E0,'.-','LineWidth',3,'color',lcolor(1,:),'MarkerSize',20,'MarkerFaceColor',lcolor(2,:),'MarkerEdgeColor',lcolor(2,:))
    hold on
    box on
    grid on
    for p=2:n_dataset
        plot(R_in(p).reaction_coordinates(:,end),R_in(p).energies(:,end)-E0,'.','MarkerSize',20,'color',lcolor(p+1,:))
    end
    xlabel('Reaction coordinate')
    ylabel('Energy (eV)')
    set(gca,'fontsize',12,'fontname','cambria math')
    xticks(linspace(R_in(1).reaction_coordinates(1,end),R_in(1).reaction_coordinates(end,end),11))
    xticklabels({'Reactants','','','','','','','','','','Products'})
    if ~static
        moving_pt=plot(R_in(end).reaction_coordinates(1,end),R_in(end).energies(1,end)-E0,'.','MarkerSize',32,'color',lcolor(6,:));
    else
        plot(R_in.reaction_coordinates(idxTS,end),R_in.energies(idxTS,end)-E0,'.','MarkerSize',20,'MarkerFaceColor',lcolor(3,:),'MarkerEdgeColor',lcolor(3,:))
    end
    xlim([R_in(1).reaction_coordinates(1,end),R_in(1).reaction_coordinates(end,end)])
end

% convert aLim, bLim and cLim to XLim, YLim and ZLim
XLim(1)=[aLim(1) bLim(1) cLim(1)]*POSCARs(1).vec(:,1);
YLim(1)=[aLim(1) bLim(1) cLim(1)]*POSCARs(1).vec(:,2);
ZLim(1)=[aLim(1) bLim(1) cLim(1)]*POSCARs(1).vec(:,3);
XLim(2)=[aLim(2) bLim(2) cLim(2)]*POSCARs(1).vec(:,1);
YLim(2)=[aLim(2) bLim(2) cLim(2)]*POSCARs(1).vec(:,2);
ZLim(2)=[aLim(2) bLim(2) cLim(2)]*POSCARs(1).vec(:,3);
% find the absolute limits of the rendering across all POSCARs structure to avoid resizing between frames
xl=XLim;
yl=YLim;
zl=ZLim;
for p=1:n_images
    xyz=POSCARs(p).positions;
    xyz(find(~prod((POSCARs(p).xred>[aLim(1) bLim(1) cLim(1)]).*(POSCARs(p).xred<[aLim(2) bLim(2) cLim(2)]),2)),:)=NaN;
    xl=[min([xl(1) min(xyz(:,1))]) max([xl(2) max(xyz(:,1))])];
    yl=[min([yl(1) min(xyz(:,2))]) max([yl(2) max(xyz(:,2))])];
    zl=[min([zl(1) min(xyz(:,3))]) max([zl(2) max(xyz(:,3))])];
end
CartLim=[xl;yl;zl]+[-0.5 0.5];
% Rendering of the frame(s)
if static
    if EnergyProfile
        subplot(1,2,1);
        molecule3D(R_in(end).XDATCAR(idxTS,end),aLim,bLim,cLim,'rotate_view_angle',default_view_angle+rotate_view_angle)
    else
        error(sprintf(['The static option is only meant to display an energy profile along with a 3D rendering of the model.\n...' ...
            'If you want a simple rendering of the model, please use molecule3D instead.']))
    end
else
    nbytes = fprintf('Rendering model 0 out of %d', n_images);
    for p=1:n_images
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('Rendering model %d out of %d\n', p, n_images);
        if EnergyProfile
            subplot(1,2,1);
        end
        molecule3D(POSCARs(p),aLim,bLim,cLim,'lights','off','rotate_view_angle',default_view_angle+rotate_view_angle,'lattice_vec',lattice_vec)
        % Add lights
        cl=camlight(0,0);
        % use the absolute limits of the rendering
        xlim(CartLim(1,:))
        ylim(CartLim(2,:))
        zlim(CartLim(3,:))
        % only the "moving_pt" graphical object is updated, to track the energy of the 3D model currently rendered
        if EnergyProfile
            subplot(1,2,2);
            set(moving_pt, 'XData', R_in(end).reaction_coordinates(p,end), 'YData', R_in(end).energies(p,end)-E0);
        end
        % save current figure as a frame in a frame array
        M(p)=getframe(gcf);
    end
    if save_data
        delay_time=0.05;
        % conversion of the frames into a GIF file
        frames2gif(M,'ReactionRendering',gif_direction,delay_time)
    end
end
% optionnally, the frames are set as the output of the function
if nargout > 0
    frames=M;
end
end
