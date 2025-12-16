function Freq=readFrequencies(path,varargin)
%==================================================================================================================================%
% readEnergy.m: Extraction of the frequencies of a system from the output of HIVE (v0.2.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.2 (11/09/2025) - The recording of the grep search is now optional, and default is no save
%       author: EYG
%   version 0.2.1 (12/09/2025) - Add checking of the last character of the path against 
%==================================================================================================================================%
% args:
%   path:   Location of the directory where the PHONONOUT.hive file is stored
%==================================================================================================================================%
if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
    path=[path,'/'];
end
save=false;
% read optional arguments
if exist('varargin')
    for p=1:2:length(varargin)
        switch varargin{p}
            case 'save'
                save=varargin{p+1};
        end
    end
end
% check if the OUTCAR file has already been processed. If not, use the home-made grep function for MatLab
if ~exist([path,'PHONONOUT.hive'])
    HIVE_analysis('path',path);
end
grep([path,'PHONONOUT.hive'],')    ','prt',[path,'Freq_HIVE.dat']);

k=1;
Data=readlines([path,'Freq_HIVE.dat']);
while ~isempty(Data{k})
    Freq(k)=str2double(Data{k}(11:21))*1e12;
    % Due to numerical instabilities, "true" zero frequencies are often shown as very small frequencies (either positive or
    % negative). As such a threshold of 1e8 Hz is applied.
    if abs(Freq(k))<1e8
        Freq(k)=0;
    end
    k=k+1;
end
if ~save
    delete([path,'Freq_HIVE.dat'])
end
end