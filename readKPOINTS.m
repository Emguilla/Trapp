function KPOINTS=readKPOINTS(filename,varargin)
%==================================================================================================================================%
% readPOSCAR.m: Read of a KPOINTS file from VASP and creation of a KPOINTS structure file in MatLab (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (16/12/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   filename:   path + name of the file to be read as a KPOINTS (works for IBZKPT as well)
%==================================================================================================================================%
% Initialisation of the default parameters
KPOINTS.Title='';
KPOINTS.mode='';
KPOINTS.coord=NaN;
path='./';
if nargin==0
    filename='KPOINTS';
end

% Reading of the optional argument
if exist('varargin')
    for p=1:2:length(varargin)
        switch varargin{p}
            case 'path'
                path=varargin{p+1};
                if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
                    path=[path,'/'];
                end
        end
    end
end

% file opening
fid=fopen([path,filename],'r');

% Reading of the title
Data=fgetl(fid);
KPOINTS.Title=Data;
% Check if it is in an IBZKPT file format
if strcmpi('Auto',KPOINTS.Title(1:4))||strcmpi(filename(end-5:end),'IBZKPT')
    KPOINTS.mode='Explicit';
    % Reading of the number of k-points
    Data=textscan(fid,'%d',1,'commentStyle','%');
    KPOINTS.N=Data{1}(1);
    % Reading of the types of coordinates
    Data=fgetl(fid);
    Data=fgetl(fid);
    KPOINTS.coord=lower(Data(1));
    Data=textscan(fid,'%f %f %f %d\n',KPOINTS.N,'commentStyle','%');
    KPOINTS.k=[Data{1} Data{2} Data{3}];
    KPOINTS.weights=Data{4};
else
    % Reading of the number of k-points
    Data=textscan(fid,'%d',1,'commentStyle','%');
    KPOINTS.N=Data{1}(1);
    % Reading of the types of k-point generation
    Data=fgetl(fid);
    Data=fgetl(fid);
    if strcmpi(Data(1),'l') % k-points along lines for band structure calculations
        KPOINTS.mode='Line';
        % Reading of the types of coordinates
        Data=fgetl(fid);
        KPOINTS.coord=lower(Data(1));
        % Reading of the High-Symmetry point defining each line
        Data=textscan(fid,'%f',3,'commentStyle','%');
        fgetl(fid);
        k=0;
        while length(Data{1})==3
            k=k+1;
            KPOINTS.k(k,:)=Data{1};
            Data=textscan(fid,'%f',3,'commentStyle','%');
            fgetl(fid);
        end
    else % k-points mesh according to Monkhorst and Pack
        if strcmpi(Data(1),'m')
            KPOINTS.mode='Monkhorst-Pack';
        elseif strcmpi(Data(1),'g')
            KPOINTS.mode='Gamma-centered';
        end
        % Reading of the sampling
        Data=textscan(fid,'%d',3,'commentStyle','%');
        KPOINTS.N=Data{1};
        % Reading of the offset, if any
        Data=textscan(fid,'%f',3,'commentStyle','%');
        if length(Data{1})==3
            KPOINTS.offset=Data{1};
        else
            KPOINTS.offset=NaN;
        end
    end
end
