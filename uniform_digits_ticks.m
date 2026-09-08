function uniform_digits_ticks(ax,xyz)
%==================================================================================================================================%
% uniform_digits_ticks.m:   Uniformisation of the number of significant digits in plots (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (08/09/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   ax:     GCA of the plot for which the modifications should be applied
%   xyz:    axes on which the change should be applied
%==================================================================================================================================%

for p=1:length(xyz)
    % For each directions specified by the user, the ticklabels are retrieved
    axis=get(ax,[upper(xyz(p)),'Axis']);
    tl=axis.TickLabels();
    % Each of the ticklabels are analysed to find what is their format ('fN.d' formatting)
    for p=1:length(tl)
        fN(p)=0;
        fd(p)=NaN;
        for q=1:length(tl{p})
            if strcmpi(tl{p}(q),'-')
                fN(p)=0;
            elseif strcmpi(tl{p}(q),'.')
                fd(p)=0;
            elseif isnan(fd(p))
                fN(p)=fN(p)+1;
            else
                fd(p)=fd(p)+1;
            end
        end
    end
    % The format of the longest ticklabels (excluding the minus sign) is used as a template for all the ticklabels.
    set(axis,'TickLabelFormat',['%.',num2str(max(fd)),'f'])
end