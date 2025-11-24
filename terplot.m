function plt=terplot(mode,c1,c2,c3,varargin)
%==================================================================================================================================%
% terplot.m: ternary diagram plot function (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (30/07/2025) - Creation
%       author: EYG
%   version 0.1.1 (31/07/2025) - Specify the colormap type explicitly to obtain the correct colorbar scale
%       contrib: EYG
%   version 0.1.2 (04/08/2025) - Add a cap for the colorscale to keep the scale between 0 and 255
%       contrib: EYG
%   version 0.2 (18/09/2025) - The colormap is now a tiling of triangle that uses the same colorscale as the hexagons that used to
%       author: EYG             tile the ternary diagram. The values to be entered as c1, c2 and c3 are matrices, which can be seen
%                               as a list of intervals in the three different concentrations.
%==================================================================================================================================%
% args:
%   mode:   'init' -> initialise the ternary plot function
%               c1,c2 and c3 must be character strings that label each axis. Call to this mode is MANDATORY to start the plot.
%           'append' -> Data visualisation
%               c1 and c2 are the concentration of the first two species specified during initialisation and the third is 
%               calculated as 1-c1-c2. To give a certain height for each point, use the c3 argument. Otherwise, the argument c3 
%               must be a NaN.
%               opt. args:  'color', followed by a RGB triplet.
%                           '-','--' or '-.' to choose the type of line, followed by the width of the line (in pxl)
%                           '.' or '^' followed by the size (in pxl) of the dot or triangle, multiplied by 4 (matlab thing)
%           'append_surf' -> 2-D heightmap visualisation
%               c1, c2 and c3 are the intervals of concentration of the different species specified during initialisation. Although 
%               optional, the 'height' argument is here mandatory as it provides the intensity associated with the area defined by 
%               the three intervals.
%               opt. args:  'height', followed by the intensity value associated with each area. !!! Mandatory !!!
%                           'colormap', followed by the type of the desired color palette (default is 'jet')
%                           'clim', followed by the lower and upper boundary of the color scale.
%==================================================================================================================================%
z_axis=true;
if strcmpi(mode,'init')
    h=sqrt(3)/2;
    A=[0 0];
    B=[0.5 h];
    C=[1 0];
    
    AB=[linspace(A(1),B(1),11)' linspace(A(2),B(2),11)'];
    BC=[linspace(B(1),C(1),11)' linspace(B(2),C(2),11)'];
    CA=[linspace(C(1),A(1),11)' linspace(C(2),A(2),11)'];
    triangle=nsidedpoly(3,'center',[0.5 h/3],'SideLength',1);
    plt=figure;
    plot(triangle,'FaceColor',[1 1 1],'EdgeColor','none') %background is a white triangle
    hold on
    d1=[cos(0) sin(0)];
    d2=[cos(pi/6) sin(pi/6)];
    d3=[cos(pi/3) sin(pi/3)];
    d4=[cos(pi/2) sin(pi/2)];
    d5=[cos(2*pi/3) sin(2*pi/3)];
    d6=[cos(5*pi/6) sin(5*pi/6)];
    
    % Generation of the tick labels along the three axis, with a rotation of the text to align with the triangular grid
    text(AB(1,1) -d1(1)*0.1,AB(1,2) -d1(2)*0.1,'    0','rotation',0,'EdgeColor','none','fontname','cambria math')
    text(AB(2,1) -d1(1)*0.1,AB(2,2) -d1(2)*0.1,'  10','rotation' ,0,'EdgeColor','none','fontname','cambria math')
    text(AB(3,1) -d1(1)*0.1,AB(3,2) -d1(2)*0.1,'  20','rotation' ,0,'EdgeColor','none','fontname','cambria math')
    text(AB(4,1) -d1(1)*0.1,AB(4,2) -d1(2)*0.1,'  30','rotation' ,0,'EdgeColor','none','fontname','cambria math')
    text(AB(5,1) -d1(1)*0.1,AB(5,2) -d1(2)*0.1,'  40','rotation' ,0,'EdgeColor','none','fontname','cambria math')
    text(AB(6,1) -d1(1)*0.1,AB(6,2) -d1(2)*0.1,'  50','rotation' ,0,'EdgeColor','none','fontname','cambria math')
    text(AB(7,1) -d1(1)*0.1,AB(7,2) -d1(2)*0.1,'  60','rotation' ,0,'EdgeColor','none','fontname','cambria math')
    text(AB(8,1) -d1(1)*0.1,AB(8,2) -d1(2)*0.1,'  70','rotation' ,0,'EdgeColor','none','fontname','cambria math')
    text(AB(9,1) -d1(1)*0.1,AB(9,2) -d1(2)*0.1,'  80','rotation' ,0,'EdgeColor','none','fontname','cambria math')
    text(AB(10,1)-d1(1)*0.1,AB(10,2)-d1(2)*0.1,'  90','rotation' ,0,'EdgeColor','none','fontname','cambria math')
    text(AB(11,1)-d1(1)*0.1,AB(11,2)-d1(2)*0.1,'100','rotation'  ,0,'EdgeColor','none','fontname','cambria math');
    
    text(BC(1,1) +d3(1)*0.025,BC(1,2) +d3(2)*0.025,'    0','rotation',60,'EdgeColor','none','fontname','cambria math')
    text(BC(2,1) +d3(1)*0.025,BC(2,2) +d3(2)*0.025,'  10','rotation' ,60,'EdgeColor','none','fontname','cambria math')
    text(BC(3,1) +d3(1)*0.025,BC(3,2) +d3(2)*0.025,'  20','rotation' ,60,'EdgeColor','none','fontname','cambria math')
    text(BC(4,1) +d3(1)*0.025,BC(4,2) +d3(2)*0.025,'  30','rotation' ,60,'EdgeColor','none','fontname','cambria math')
    text(BC(5,1) +d3(1)*0.025,BC(5,2) +d3(2)*0.025,'  40','rotation' ,60,'EdgeColor','none','fontname','cambria math')
    text(BC(6,1) +d3(1)*0.025,BC(6,2) +d3(2)*0.025,'  50','rotation' ,60,'EdgeColor','none','fontname','cambria math')
    text(BC(7,1) +d3(1)*0.025,BC(7,2) +d3(2)*0.025,'  60','rotation' ,60,'EdgeColor','none','fontname','cambria math')
    text(BC(8,1) +d3(1)*0.025,BC(8,2) +d3(2)*0.025,'  70','rotation' ,60,'EdgeColor','none','fontname','cambria math')
    text(BC(9,1) +d3(1)*0.025,BC(9,2) +d3(2)*0.025,'  80','rotation' ,60,'EdgeColor','none','fontname','cambria math')
    text(BC(10,1)+d3(1)*0.025,BC(10,2)+d3(2)*0.025,'  90','rotation' ,60,'EdgeColor','none','fontname','cambria math')
    text(BC(11,1)+d3(1)*0.025,BC(11,2)+d3(2)*0.025,'100','rotation'  ,60,'EdgeColor','none','fontname','cambria math');
    
    text(CA(1,1) -d4(1)*0.025,CA(1,2) -d4(2)*0.025,'    0','rotation',-60,'EdgeColor','none','fontname','cambria math')
    text(CA(2,1) -d4(1)*0.025,CA(2,2) -d4(2)*0.025,'  10','rotation' ,-60,'EdgeColor','none','fontname','cambria math')
    text(CA(3,1) -d4(1)*0.025,CA(3,2) -d4(2)*0.025,'  20','rotation' ,-60,'EdgeColor','none','fontname','cambria math')
    text(CA(4,1) -d4(1)*0.025,CA(4,2) -d4(2)*0.025,'  30','rotation' ,-60,'EdgeColor','none','fontname','cambria math')
    text(CA(5,1) -d4(1)*0.025,CA(5,2) -d4(2)*0.025,'  40','rotation' ,-60,'EdgeColor','none','fontname','cambria math')
    text(CA(6,1) -d4(1)*0.025,CA(6,2) -d4(2)*0.025,'  50','rotation' ,-60,'EdgeColor','none','fontname','cambria math')
    text(CA(7,1) -d4(1)*0.025,CA(7,2) -d4(2)*0.025,'  60','rotation' ,-60,'EdgeColor','none','fontname','cambria math')
    text(CA(8,1) -d4(1)*0.025,CA(8,2) -d4(2)*0.025,'  70','rotation' ,-60,'EdgeColor','none','fontname','cambria math')
    text(CA(9,1) -d4(1)*0.025,CA(9,2) -d4(2)*0.025,'  80','rotation' ,-60,'EdgeColor','none','fontname','cambria math')
    text(CA(10,1)-d4(1)*0.025,CA(10,2)-d4(2)*0.025,'  90','rotation' ,-60,'EdgeColor','none','fontname','cambria math')
    text(CA(11,1)-d4(1)*0.025,CA(11,2)-d4(2)*0.025,'100','rotation'  ,-60,'EdgeColor','none','fontname','cambria math');
    
    % Creation of the axis label and rotation
    t=text(0,0,c1,'rotation' ,0,'EdgeColor','none','fontname','cambria math');
    pos_t(:)=CA(6,:)-d4*0.125-d1*t.Extent(3)/2;
    set(t,'Position',[pos_t(1) pos_t(2)])
    t=text(0,0,c2,'rotation' ,0,'EdgeColor','none','fontname','cambria math');
    pos_t(1:2)=AB(6,:)+d6*0.125-d3*t.Extent(3)/2;
    set(t,'Position',[pos_t(1) pos_t(2)])
    set(t,'Rotation',60)
    t=text(0,0,c3,'rotation' ,0,'EdgeColor','None','fontname','cambria math');
    pos_t(1:2)=BC(6,:)+d2*0.125+d5*t.Extent(3)/2;
    set(t,'Position',[pos_t(1) pos_t(2)])
    set(t,'Rotation',-60)

    % Drawing of the grid in gray
    nsub=10;
    for p=1:(nsub-1)
        plot([AB(p+1,1) BC(nsub+1-p,1)],[AB(p+1,2) BC(nsub+1-p,2)],'--','color',[0.85 0.85 0.85])
        plot([BC(p+1,1) CA(nsub+1-p,1)],[BC(p+1,2) CA(nsub+1-p,2)],'--','color',[0.85 0.85 0.85])
        plot([CA(p+1,1) AB(nsub+1-p,1)],[CA(p+1,2) AB(nsub+1-p,2)],'--','color',[0.85 0.85 0.85])
    end
    
    % Drawing of the axis delimitation of the plot
    plot(AB(:,1),AB(:,2),'k','LineWidth',2);
    pbaspect([1 1 1])
    axis off
    plot(BC(:,1),BC(:,2),'k','LineWidth',2);
    plot(CA(:,1),CA(:,2),'k','LineWidth',2);
    
    xlim([0 1])
    ylim([h/2-0.5 h/2+0.5])
    
    % Set the font and fontsize
    set(gca,'fontsize',12,'fontname','cambria math')
elseif strcmpi(mode,'append')
    if isnan(c3)
        z_axis=false;
    end
    line=false;
    dot=false;
    clrt_defined=false;
    % Check if there are optional arguments (line/dot and color are possible)
    if exist('varargin','var')
        for p=1:2:length(varargin)
            switch varargin{p}
                case '-'
                    line=true;
                    dot=false;
                    LW_type=varargin{p};
                    LW=varargin{p+1};
                case '--'
                    line=true;
                    dot=false;
                    LW_type=varargin{p};
                    LW=varargin{p+1};
                case '-.'
                    line=true;
                    dot=false;
                    LW_type=varargin{p};
                    LW=varargin{p+1};
                case '.'
                    line=false;
                    dot=true;
                    MS_type=varargin{p};
                    MS=varargin{p+1};
                case '^'
                    line=false;
                    dot=true;
                    MS_type=varargin{p};
                    MS=varargin{p+1};
                case 'color'
                    color=varargin{p+1};
                    clrt_defined=true;
            end
        end
    end

    % remapping of the positions from cartesian to ternary coordinates
    x=1-c1-c2/2;
    y=sqrt(3)/2*c2;
    % drawing of the lines/dots in the specified color, if defined
    if z_axis
        if line
            if clrt_defined
                plot3(x,y,c3,LW_type,'LineWidth',LW,'color',color)
            else
                plot3(x,y,c3,LW_type,'LineWidth',LW)
            end
        elseif dot
            if clrt_defined
                plot3(x,y,c3,MS_type,'MarkerSize',MS,'color',color)
            else
                plot3(x,y,c3,MS_type,'MarkerSize',MS)
            end
        else
            if clrt_defined
                plot3(x,y,c3,'.','MarkerSize',8,'color',color)
            else
                plot3(x,y,c3,'.','MarkerSize',8)
            end
        end
    else
        if line
            if clrt_defined
                plot(x,y,LW_type,'LineWidth',LW,'color',color)
            else
                plot(x,y,LW_type,'LineWidth',LW)
            end
        elseif dot
            if clrt_defined
                plot(x,y,MS_type,'MarkerSize',MS,'color',color)
            else
                plot(x,y,MS_type,'MarkerSize',MS)
            end
        else
            if clrt_defined
                plot(x,y,'.','MarkerSize',8,'color',color)
            else
                plot(x,y,'.','MarkerSize',8)
            end
        end
    end
elseif strcmpi(mode,'append_surf')
    clim_defined=false;
    clrmap_defined=false;
    dx=0.01;
    load('colormaps.mat');
    clrmap=c_jet;
    colormap jet
    % Check if there are optional arguments (type of color palette and limits of the color scale)
    if exist('varargin','var')
        for p=1:2:length(varargin)
            switch varargin{p}
                case 'colormap'
                    switch varargin{p+1}
                        case 'parula'
                            clrmap=c_parula;
                            colormap parula
                        case 'turbo'
                            clrmap=c_turbo;
                            colormap turbo
                        case 'hsv'
                            clrmap=c_hsv;
                            colormap hsv
                        case 'hot'
                            clrmap=c_hot;
                            colormap hot
                        case 'cool'
                            clrmap=c_cool;
                            colormap cool
                        case 'spring'
                            clrmap=c_spring;
                            colormap spring
                        case 'summer'
                            clrmap=c_summer;
                            colormap summer
                        case 'autumn'
                            clrmap=c_autumn;
                            colormap autumn
                        case 'winter'
                            clrmap=c_winter;
                            colormap winter
                        case 'gray'
                            clrmap=c_gray;
                            colormap gray
                        case 'bone'
                            clrmap=c_bone;
                            colormap bone
                        case 'copper'
                            clrmap=c_copper;
                            colormap copper
                        case 'pink'
                            clrmap=c_pink;
                            colormap pink
                        case 'sky'
                            clrmap=c_sky;
                            colormap sky
                        case 'abyss'
                            clrmap=c_abyss;
                            colormap abyss
                        case 'jet'
                            clrmap=c_jet;
                            colormap jet
                    end
                case 'clim'
                    clim_defined=true;
                    clim=varargin{p+1};
                case 'discretisation'
                    dx=varargin{p+1};
                case 'height'
                    c4=varargin{p+1};
            end
        end
    end
    if ~clrmap_defined
        clrmap=c_jet;
    end
    if clim_defined
        c_min=clim(1);
        c_max=clim(2);
    else
        c_min=min(c4);
        c_max=max(c4);
    end
    % computation of the color scale based on the most extreme values of the height argument
    c_axis=round((c4-c_min)/(c_max-c_min)*255);
    for p=1:length(c1(:,1))
        for q=1:length(c_axis)
            if c_axis(q)<0
                c_axis(q)=0;
            elseif c_axis(q)>255
                c_axis(q)=255;
            end
        end
        if ~isnan(c4(p))
            reverted=false;
            normal=false;
            % Extract the interval for the current line of the list of intervals
            c1min=min(c1(p,:));
            c1max=max(c1(p,:));
            c2min=min(c2(p,:));
            c2max=max(c2(p,:));
            c3min=min(c3(p,:));
            c3max=max(c3(p,:));
            % The area to be drawn is a triangle formed by three points, themselves being the only three-way intersection of the
            % lines of the interval boundary in the ternary diagram.
            % Two possibility exists: either the triangle rests on its basis, or on its tip. To discriminate the two cases, we
            % resolve analytically the x-position of the intersection of the lines. For a specific area, only one of the two
            % following matrices has three times the same column, which tells us whether the triangle is on its tip or its base.
            reverted_mat=[  1-c1min-c2max/2 c2max/2+c3max 0.5*(1-c1min+c3max);
                            1-c1max-c2min/2 c2min/2+c3max 0.5*(1-c1max+c3max);
                            1-c1max-c2max/2 c2max/2+c3min 0.5*(1-c1max+c3min)];
            normal_mat=[    1-c1min-c2min/2 c2min/2+c3max 0.5*(1-c1min+c3max);
                            1-c1min-c2max/2 c2max/2+c3min 0.5*(1-c1min+c3min);
                            1-c1max-c2min/2 c2min/2+c3min 0.5*(1-c1max+c3min)];
            if sum((reverted_mat-reverted_mat(:,1))<-eps)+sum((reverted_mat-reverted_mat(:,1))>eps)==0
                reverted=true;
            end
            if sum((normal_mat-normal_mat(:,1))<-eps)+sum((normal_mat-normal_mat(:,1))>eps)==0
                normal=true;
            end
            % dummy triangle
            triangle=nsidedpoly(3,'center',[0 0],'SideLength',dx);
            % re-definition of the vertices of the triangle depending on its orientation
            if reverted
                triangle.Vertices=[1-c1min-c2max/2 sqrt(3)*c2max/2; 1-c1max-c2min/2 sqrt(3)*c2min/2; 1-c1max-c2max/2 sqrt(3)*c2max/2];
            elseif normal
                vertices=[1-c1min-c2min/2 sqrt(3)*c2min/2; 1-c1min-c2max/2 sqrt(3)*c2max/2; 1-c1max-c2min/2 sqrt(3)*c2min/2];
                [~,idx]=sort(vertices(:,2));
                triangle.Vertices=vertices(idx,:);
            end
            % each datapoint is assigned a color according to the specified palette, and is represented by a triangle.
            plot(triangle,'FaceColor',clrmap(c_axis(p)+1,:),'EdgeColor','none','FaceAlpha',1)
        end
    end
end


