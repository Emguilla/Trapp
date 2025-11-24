function Q=Qv_div(Q1,Q2)
%==================================================================================================================================%
% Qv_div.m: Division of vibration partition function arrays (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   Q1, Q2: vibration partition function arrays as obtained from Qv.m
%==================================================================================================================================%
% First, a term-by-term division is applied
for p=1:min([length(Q1) length(Q2)])
    Qtmp(p)=Q1(p)/Q2(p);
end
% Second, the remaining unpaired terms are included
for p=min([length(Q1) length(Q2)])+1:max([length(Q1) length(Q2)])
    if length(Q1)>length(Q2)
        Qtmp(p)=Q1(p);
    elseif length(Q1)<length(Q2)
        Qtmp(p)=1/Q2(p);
    end
end
% The partition function is the product of each term-by-term division
Q=prod(Qtmp);
end