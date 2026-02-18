function E=ZPE(nu)
%==================================================================================================================================%
% ZPE.m:    computation of the zero-point energy of a system using the frequencies of its vibrational mode (v0.1.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.1.1 (18/02/2026) - Added an exception resulting in E=0 when the frequency array is empty.
%       contrib: EYG
%==================================================================================================================================%
% args:
%   nu: Array of frequencies
%==================================================================================================================================%
load('constant_fund.mat','h')
if ~isempty(nu)
    E=0.5*h*sum(nu);
else
    E=0;
end