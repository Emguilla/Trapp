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
EIGENVAL.Title='';
EIGENVAL.Ispin=false;
EIGENVAL.EFermi=NaN;
path='./';
if nargin==0
    filename='EIGENVAL';
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

% Reading of the header information: spin-polarisation
Data=textscan(fid,'%d',4,'commentStyle','%');
EIGENVAL.Ispin=Data{1}(4)==2;

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
EIGENVAL.Title=Data{1};

EIGENVAL.EFermi=NaN;
% Reading of the header information: Number of bands and k-points
Data=textscan(fid,'%d',3,'commentStyle','%');
nbands=Data{1}(3);
nkpts=Data{1}(2);

if EIGENVAL.Ispin% If ISPIN is on, the data contains the spin up and down energies of the bands and their occupation
    for p=1:nkpts
        fgetl(fid);
        Data=textscan(fid,'%f',4,'commentStyle','%');
        EIGENVAL.k(p,:)=Data{1}(1:3);
        for q=1:nbands
            Data=textscan(fid,'%d %f %f %f %f',1,'commentStyle','%');
            EIGENVAL.E_u(p,q)=Data{2};
            EIGENVAL.E_d(p,q)=Data{3};
            EIGENVAL.Occ_u(p,q)=Data{4};
            EIGENVAL.Occ_d(p,q)=Data{5};
        end
        % As a proxy, the Fermi energy can be thought as the limit for which the band are less than half occupied
        EIGENVAL.EFermi=max([max(EIGENVAL.E_u(p,:).*(EIGENVAL.Occ_u(p,:)>0.5)) max(EIGENVAL.E_d(p,:).*(EIGENVAL.Occ_d(p,:)>0.5)) EIGENVAL.EFermi]); % def foireuse
    end
else% If ISPIN is off, the data contains the energies of the bands and their occupation irrespective of the spin
    for p=1:nkpts
        fgetl(fid);
        Data=textscan(fid,'%f',4,'commentStyle','%');
        EIGENVAL.k(p,:)=Data{1}(1:3);
        for q=1:nbands
            Data=textscan(fid,'%d %f %f',1,'commentStyle','%');
            EIGENVAL.E(p,q)=Data{2};
            EIGENVAL.Occ(p,q)=Data{3};
        end
        % As a proxy, the Fermi energy can be thought as the limit for which the band are less than half occupied
        EIGENVAL.EFermi=max([max(EIGENVAL.E(p,:).*(EIGENVAL.Occ(p,:)>0.5)) EIGENVAL.EFermi]); % def foireuse
    end
end
