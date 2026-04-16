function SC=POSCARsupercell(UC,nx,ny,nz)
%==================================================================================================================================%
% POSCARsupercell.m:    duplicate a unit cell to create a supercell (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (25/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   UC:             POSCAR structure
%   nx, ny, nz:     number of replicate in each directions
%==================================================================================================================================%
% If the POSCAR input is a character string, the corresponding file is read into a POSCAR structure
if ischar(UC)
    UC=readPOSCAR(UC);
end
% Repetition of the POSCAR along "x" to generate a "line" of cells
SC_x=UC;
for px=1:nx-1
    POSCAR=UC;
    POSCAR.positions=POSCAR.positions+px*POSCAR.vec(1,:);
    SC_x=appendPOSCAR(SC_x,POSCAR);
end
% Repetition of the POSCAR along "y" to generate a "plane" of cells
SC_xy=SC_x;
for py=1:ny-1
    POSCAR=SC_x;
    POSCAR.positions=POSCAR.positions+py*POSCAR.vec(2,:);
    SC_xy=appendPOSCAR(SC_xy,POSCAR);
end
% Repetition of the POSCAR along "z" to generate a "volume" of cells, i.e. the supercell
SC=SC_xy;
for pz=1:nz-1
    POSCAR=SC_xy;
    POSCAR.positions=POSCAR.positions+pz*POSCAR.vec(3,:);
    SC=appendPOSCAR(SC,POSCAR);
end
% Correction of the lattice vector to reflect the duplications of the cells
SC.vec(1,:)=SC.vec(1,:)*nx;
SC.vec(2,:)=SC.vec(2,:)*ny;
SC.vec(3,:)=SC.vec(3,:)*nz;
% write to a POSCAR file
writePOSCAR(SC,['POSCAR_',num2str(nx),'x',num2str(ny),'x',num2str(nz)])