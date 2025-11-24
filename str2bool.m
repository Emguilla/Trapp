function [out]=str2bool(string)
%==================================================================================================================================%
% str2bool.m:   conversion of a character array of "T" and "F" into boolean matrix (v0.5)
%==================================================================================================================================%
% Version history:
%   version 0.1 (10/09/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   string: character array containing only "T"s and "F"s
%==================================================================================================================================%
szstr=size(string);
for p=1:szstr(1)
    for q=1:szstr(2)
        if strcmpi(string(p,q),'t')
            out(p,q)=true;
        elseif strcmpi(string(p,q),'f')
            out(p,q)=false;
        else
            error('Input must be a character vector including only "T" or "F".')
        end
    end
end