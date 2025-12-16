function Q=Qr(molecule,T)
%==================================================================================================================================%
% Qt.m: Calculation of the rotational partition function (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (14/08/2025) - Creation
%       author: EYG
%   version 0.2 (20/08/2025) - Modification of the case where the argument is a POSCAR structure to use the mass now stored in the
%       author: EYG            POSCAR structure
%==================================================================================================================================%
% args:
%   molecule:   - character string:   type of radical/molecule investigation (only H2, CH3 and CH4 are implemented)
%               - POSCAR structure:   Can be used for larger CH molecules /!\ experimental /!\
%   T:          Temperature in K
%==================================================================================================================================%
load('constant_fund.mat','kB','uma','hbar','ptable');
% average mass of the hydrogen and carbon atoms
mH=ptable.mass(1)*uma;
mC=ptable.mass(6)*uma;
if ischar(molecule) % Case of a well-defined molecule
    switch molecule
        case 'H2' % Calculation of the rotational partition function as rigid rotor
            a_HH=0.75049e-10;
            Q=0;
            sigma=2;
            I=2*mH*(0.5*a_HH)^2;
            B=(hbar)^2/(2*I);
            J=-1;
            dQ=1;
            while abs(dQ)>1e-12
                J=J+1;
                dQ=((2*J+1)/sigma)*exp(-B*J*(J+1)/(kB*T));
                Q=Q+dQ;
            end
        case 'CH3' % Calculation of the rotational partition function as a planar top
            a_CH=1.08621e-10;
            sigma=6;
            I_A=3*mH*(a_CH)^2;
            I_B=2*mH*(sqrt(3)/2*a_CH)^2;
            I_C=mH*(a_CH)^2+2*mH*(a_CH/2)^2;
            A=hbar^2/(2*I_A);
            B=hbar^2/(2*I_B);
            C=hbar^2/(2*I_C);
            Q=(sqrt(pi)/sigma)*sqrt((kB*T)^3/(A*B*C));
        case 'CH4' % Calculation of the rotational partition function as a spherically symetrical top
            a_CH=1.0981e-10;
            sigma=12;
            I=8/3*mH*(a_CH)^2;
            B=hbar^2/(2*I);
            Q=(sqrt(pi)/sigma)*sqrt((kB*T/(B))^3);
    end
elseif isstruct(molecule) % Case of an asymetrical molecule
    % Determination of the centre of mass
    x_CM=sum(molecule.positions(:,1).*molecule.mass);
    y_CM=sum(molecule.positions(:,2).*molecule.mass);
    z_CM=sum(molecule.positions(:,3).*molecule.mass);
    r_CM=[x_CM y_CM z_CM];
    % Calculation of the inertia moment wrt the three cartesian axis
    for p=1:sum(molecule.n_chemicals)
        I_x(p)=molecule.mass(p)*norm(molecule.positions(p,[2 3])-r_CM([2 3]))^2;
        I_y(p)=molecule.mass(p)*norm(molecule.positions(p,[1 3])-r_CM([1 3]))^2;
        I_z(p)=molecule.mass(p)*norm(molecule.positions(p,[1 2])-r_CM([1 2]))^2;
    end
    I=[sum(I_x) sum(I_y) sum(I_z)];
    % Calculation of the rotational constant (spectroscopy stuff)
    A=hbar^2/(2*I(1));
    B=hbar^2/(2*I(2));
    C=hbar^2/(2*I(3));
    % Calculation of the rotational partition function
    Q=(sqrt(pi)/3)*sqrt((kB*T)^3/(A*B*C));
end
end
