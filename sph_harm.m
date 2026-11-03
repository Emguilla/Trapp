function Ylm=sph_harm(l,m,theta,phi)
    Plm=legendre(l,cos(theta));
    Plm=Plm(abs(m)+1);
    k=sqrt((2*l+1)*factorial(l-abs(m))/(4*pi*factorial(l+abs(m))));
    Ylm=(-1)^(sign(m)+1)*k*Plm*exp(1i*m*phi);
end
