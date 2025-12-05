%function readDOSCAR
%==================================================================================================================================%
% readDOSCAR.m: Read of a DOSCAR file from VASP (v0.1). Limited to non-spin polarised calculations
%==================================================================================================================================%
% Version history:
%   version 0.1 (24/11/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   filename: path + name of the file to be read as a DOSCAR
%==================================================================================================================================%
clear
close all
clc

filename='DOSCAR';
path='./';
fid=fopen([path,filename],'r');
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
Data=textscan(fid,'%f',5,'commentStyle','%');
Emin=Data{1}(2);
Emax=Data{1}(1);
nE=Data{1}(3);
EFermi=Data{1}(4)
for p=1:nE
    Data=textscan(fid,'%f',3,'commentStyle','%');
    E(p)=Data{1}(1);
    DOS(p)=Data{1}(2);
end
Data=textscan(fid,'%f',5,'commentStyle','%');
for p=1:nE
    Data=textscan(fid,'%f',10,'commentStyle','%');
    E1(p)=Data{1}(1);
    DOS_A(p)=Data{1}(2);
    DOS_B(p)=Data{1}(3);
    DOS_C(p)=Data{1}(4);
    DOS_D(p)=Data{1}(5);
    DOS_E(p)=Data{1}(6);
    DOS_F(p)=Data{1}(7);
    DOS_G(p)=Data{1}(8);
    DOS_H(p)=Data{1}(9);
    DOS_I(p)=Data{1}(10);
end
fclose(fid)
plot(E-EFermi,DOS)
hold on
% plot(E1-EFermi,DOS_A+DOS_B+DOS_C+DOS_D+DOS_E+DOS_F+DOS_G+DOS_H+DOS_I)
figure;
plot(E1-EFermi,DOS_A)
hold on
plot(E1-EFermi,DOS_B)
plot(E1-EFermi,DOS_C)
plot(E1-EFermi,DOS_D)
plot(E1-EFermi,DOS_E)
plot(E1-EFermi,DOS_F)
plot(E1-EFermi,DOS_G)
plot(E1-EFermi,DOS_H)
plot(E1-EFermi,DOS_I)
