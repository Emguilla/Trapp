function Q=Qt(molecule,V,T)
%==================================================================================================================================%
% Qt.m: Calculation of the translational partition function (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.2 (20/08/2025) - the input variable "molecule" can now be a POSCAR structure, since these structure now include the
%       author: EYG             atomic masses
%==================================================================================================================================%
% args:
%   molecule:   - character string: type of radical/molecule investigation (only H, H2, CH3 and CH4 are implemented)
%               - POSCAR structure
%   V:          Volume of interest (1e-6 cubic meter, i.e. 1 cubic centimeter), this value determines the unit in which the final
%               reaction rate coefficient is expressed
%   T:          Temperature in K
%==================================================================================================================================%
load('constant_fund.mat','kB','h','uma','ptable')
if ischar(molecule) % Shortcut for common molecules
    switch molecule
        case 'H'
            m_g=1*ptable.mass(1)*uma;
        case 'H2'
            m_g=2*ptable.mass(1)*uma;
        case 'CH3'
            m_g=3*ptable.mass(1)*uma+1*ptable.mass(6)*uma;
        case 'CH4'
            m_g=4*ptable.mass(1)*uma+1*ptable.mass(6)*uma;
    end
elseif isstruct(molecule) % General case for a cluster of atoms defined in a POSCAR structure
    m_g=sum(molecule.mass);
end
Q=V*(2*pi*m_g*kB*T/(h^2))^(3/2);
end