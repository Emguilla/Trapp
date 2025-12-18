function toggleVisibility(checkbox, plotHandle)
%==================================================================================================================================%
% toggleVisibility.m:   Handle whether a the visible field of plotHandle is on or off dynamically according to the state of a
%                       checkbox (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (17/12/2025) - Creation
%       author: EYG
%==================================================================================================================================%
if checkbox.Value
    plotHandle.Visible = 'on';
else
    plotHandle.Visible = 'off';
end
