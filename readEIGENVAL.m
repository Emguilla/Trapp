function EIGENVAL=readEIGENVAL(filename,varargin)
%==================================================================================================================================%
% readEIGENVAL.m: Read of an EIGENVAL file from VASP (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (15/12/2025) - Creation
%       author: EYG
%   version 0.2 (16/12/2025) - Handling of spin-polarised calculations
%       author: EYG
%==================================================================================================================================%
% args:
%   filename:   path + name of the file to be read as an EIGENVAL. Technically optionnal, and in such case it defaults to the file 
%               EIGENVAL in the current directory
%==================================================================================================================================%
% Initialisation of the default parameters
path='./';
if nargin==0
    filename='EIGENVAL';
end
EFermi=NaN;

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

% Reading of the header information: spin-polarisation
Data=textscan(fid,'%d',4,'commentStyle','%');
Ispin=Data{1}(4)==2;

% Reading of the header information: Volume of cell, length of lattice vectors, POTIM
Data=textscan(fid,'%f',5,'commentStyle','%');
Vcell=Data{1}(1);
norm_a=Data{1}(2);
norm_b=Data{1}(3);
norm_c=Data{1}(4);
POTIM=Data{1}(5);

% Reading of the header information: Initial temperature of the MD calculation
Data=textscan(fid,'%f',1,'commentStyle','%');
T_init=Data{1};

% Reading of the header information: Sanity check, there must be the word "CAR" on line 5
Data=textscan(fid,'%s',1,'commentStyle','%');
if ~strcmpi('car',Data{1})
    error('Something went wrong with your DOSCAR file. Line 5 should only contain the word "CAR"')
end

% Reading of the header information: Title of the system
Data=textscan(fid,'%s',1,'commentStyle','%');
Title=Data{1};

% 
Data=textscan(fid,'%d',3,'commentStyle','%');
nbands=Data{1}(3);
nkpts=Data{1}(2);
if Ispin
    for p=1:nkpts
        fgetl(fid);
        Data=textscan(fid,'%f',4,'commentStyle','%');
        k(p,:)=Data{1}(1:3);
        for q=1:nbands
            Data=textscan(fid,'%d %f %f %f %f',1,'commentStyle','%');
            E_u(p,q)=Data{2};
            E_d(p,q)=Data{3};
            Occ_u(p,q)=Data{4};
            Occ_d(p,q)=Data{5};
        end
        EFermi=max([max(E_u(p,:).*(Occ_u(p,:)>0.5)) max(E_d(p,:).*(Occ_d(p,:)>0.5)) EFermi]); % def foireuse
    end
else
    for p=1:nkpts
        fgetl(fid);
        Data=textscan(fid,'%f',4,'commentStyle','%');
        k(p,:)=Data{1}(1:3);
        for q=1:nbands
            Data=textscan(fid,'%d %f %f',1,'commentStyle','%');
            E(p,q)=Data{2};
            Occ(p,q)=Data{3};
        end
        EFermi=max([max(E(p,:).*(Occ(p,:)>0.5)) EFermi]); % def foireuse
    end
end
EIGENVAL.k=k;
EIGENVAL.Ispin=Ispin;
if Ispin
    EIGENVAL.E_u=E_u;
    EIGENVAL.E_d=E_d;
    EIGENVAL.Occ_u=Occ_u;
    EIGENVAL.Occ_d=Occ_d;
else
    EIGENVAL.E=E;
    EIGENVAL.Occ=Occ;
end
EIGENVAL.EFermi=EFermi;