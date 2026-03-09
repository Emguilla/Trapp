%==================================================================================================================================%
% ptable_generator.m:  Creation of a matlab module containing periodic table data (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (09/03/2026) - Creation
%       author: EYG
%==================================================================================================================================%
% The point of this script being to save variable into a module, everything is cleaned first
clear all
% Periodic table data
ptable.symbol={...
%     1    2    3                                                                          4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
    'H' ,                                                                                                                                                      'He',...
    'Li','Be',                                                                                                                        'B' ,'C' ,'N' ,'O' ,'F' ,'Ne',...
    'Na','Mg',                                                                                                                        'Al','Si','P' ,'S' ,'Cl','Ar',...
    'K' ,'Ca','Sc',                                                                      'Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',...
    'Rb','Sr','Y' ,                                                                      'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe',...
    'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'};
ptable.mass=[
  1.007825 ;                                                                                                                                                                                                                                                                                                                4.0026;
  6.94     ;   9.0122;                                                                                                                                                                                                                                                 10.81  ;  12.01074 ;  14.007 ;  15.999 ;  18.998 ;  20.180 ;
 22.990    ;  24.305 ;                                                                                                                                                                                                                                                 26.982 ;  28.085   ;  30.974 ;  32.06  ;  35.45  ;  39.948 ;
 39.098    ;  40.078 ;  44.956 ;                                                                                                                                             47.867 ;  50.942 ;  51.996 ;  54.938 ;  55.845 ;  58.933 ;  58.693 ;  63.546 ;  65.38  ;  69.723 ;  72.630   ;  74.922 ;  78.971 ;  79.904 ;  83.798 ;
 85.468    ;  87.62  ;  88.906 ;                                                                                                                                             91.224 ;  92.906 ;  95.95  ;  98     ; 101.07  ; 102.91  ; 106.42  ; 107.87  ; 112.41  ; 114.82  ; 118.71    ; 121.76  ; 127.60  ; 126.90  ; 131.29  ;
132.91     ; 137.33  ; 138.91  ; 140.12  ; 140.91  ; 144.24  ; 145    ; 150.36  ; 151.96  ; 157.25  ; 158.93  ; 162.50  ; 164.93  ; 167.26  ; 168.93  ; 173.05  ; 174.97  ; 178.49  ; 180.95  ; 183.84  ; 186.21  ; 190.23  ; 192.22  ; 195.08  ; 196.97  ; 200.59  ; 204.38  ; 207.2     ; 208.98  ; 209     ; 210     ; 222     ;
223        ; 226     ; 227     ; 232.04  ; 231.04  ; 238.03  ; 237    ; 244     ; 243     ; 247     ; 247     ; 251     ; 252     ; 257     ; 258     ; 259     ; 266     ; 267     ; 268     ; 269     ; 270     ; 277     ; 278     ; 281     ; 282     ; 285     ; 286     ; 289       ; 290     ; 293     ; 294     ; 294];

% By default, if the interatomic distance is not set between a specific pair of atom, the bond will not be shown
DBOND=zeros(118,118)*NaN;
DBOND(1,1)=0.8;DBOND(1,2:118)=1.25;
DBOND(3,3)=2.5;DBOND(3,8)=2.3;DBOND(3,22)=2.5;
DBOND(5,6)=1.7;DBOND(5,7)=1.6;
DBOND(6,6)=1.8;DBOND(6,7)=1.7;DBOND(6,8)=1.8;DBOND(6,15)=1.7;DBOND(6,17)=1.8;
DBOND(8,8)=1.5;DBOND(8,14)=2.6;DBOND(8,16)=1.8;DBOND(8,20)=2.6;DBOND(8,22)=2.6;
DBOND(11,11)=4;DBOND(11,17)=3;
DBOND(14,14)=2.6;
DBOND(16,42)=2.7;DBOND(16,82)=3.4;
DBOND(22,22)=2.5;
DBOND(42,42)=2.7;
for p=1:118
    for q=p:118
        DBOND(q,p)=DBOND(p,q);
    end
end
ptable.bond_length=DBOND;

% Set colour for atoms and bonds (default for unknown atoms is purple-ish). The index of the line indicate the atomic number of the
% element.
atcol=ones(118,3).*[0.9 0.5 1.0];
atcol(1,:)=[0.95 0.95 0.95];
atcol(2,:)=[0.2 1.0 1.0];
atcol(3,:)=[0.5 0.1 1.0];
atcol(4,:)=[0.1 0.5 0.1];
atcol(5,:)=[58 163 45]/255;
atcol(6,:)=[0.3 0.3 0.3];
atcol(7,:)=[176 185 230]/255;
atcol(8,:)=[1.0 0.1 0.1];
atcol(9,:)=[0.2 0.9 0.2];
atcol(10,:)=atcol(2,:);
atcol(11,:)=atcol(3,:);
atcol(12,:)=atcol(4,:);
atcol(14,:)=[240 200 160]/255;
atcol(15,:)=[1.0 0.6 0.2];
atcol(16,:)=[0.9 0.9 0.2];
atcol(17,:)=atcol(9,:);
atcol(18,:)=atcol(2,:);
atcol(19,:)=atcol(3,:);
atcol(20,:)=[0.7 0.7 0.7];
atcol(22,:)=[0.6 0.6 0.6];
atcol(26,:)=[0.9 0.5 0.1];
atcol(35,:)=[0.6 0.1 0.1];
atcol(36,:)=atcol(2,:);
atcol(37,:)=atcol(3,:);
atcol(38,:)=atcol(4,:);
atcol(42,:)=[115 173 193]/255;
atcol(53,:)=[0.4 0.1 0.7];
atcol(54,:)=atcol(2,:);
atcol(55,:)=atcol(3,:);
atcol(56,:)=atcol(4,:);
atcol(82,:)=[87,89,97]/255;
atcol(88,:)=atcol(4,:);
ptable.color=atcol;

% Set radius of atoms, considering either the van der Waals or covalent radius
covrad=ones(118,1)*0.5;
covrad(1)=0.2;
covrad(5)=0.3;
covrad(6)=0.3;
covrad(7)=0.4;
covrad(8)=0.35;
covrad(11)=0.25;
covrad(16)=0.3;
covrad(17)=0.4;
covrad(20)=0.5;
covrad(42)=0.4;
covrad(53)=0.6;
ptable.covrad=covrad;
vdwrad=ones(118,1)*1.8;
vdwrad(1)=1.2;
vdwrad(6)=1.7;
vdwrad(8)=1.52;
vdwrad(16)=1.52;
ptable.vdwrad=vdwrad;

% recording the variable into a module "constant_fund.mat", to be placed in the matlab path for easy access
clearvars -except ptable
save('ptable.mat')