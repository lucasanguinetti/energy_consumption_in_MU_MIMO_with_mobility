%function [var_vec, phin_vec] = var_analytic1(t_vec, beta, radius)
function var_vec = var_analytic_2D(t_vec, beta, radius, D)

% This function calculates the variance of the energy, as a function of
% z=T/R^2.
% beta=4,6 (I suspect it generalizes for even integers)
% const is the coefficient in the front given by
% const = 2*T*\rho^2*R^{2*beta+2}/((1-c)^2*K) which is given
%
if ((beta~=4)&&(beta~=6)&&(beta~=2))
    w = 'beta takes wrong value' %#ok<NASGU,NOPRT>
    return
end
nmax = 100;
phin_vec = zeros(nmax,1);
z_vec = D*t_vec/radius^2;
var_vec = zeros(size(z_vec));
zlen = length(z_vec);
for zz=1:zlen
    var_temp = 0;
    zval = z_vec(zz);
    kzero_vec = besselzero(1,nmax,1);
    for kk =1:nmax
        kzero_temp = kzero_vec(kk);
        
        phin = phin_fun(beta, kzero_temp);
        
        phin_vec(kk) = phin;
        zeff = zval*kzero_temp^2;
        int_temp = 1- 1.5/zeff + 2*exp(-zeff)/zeff - exp(-2*zeff)/2/zeff;
        var_temp = var_temp + 4*phin^2/kzero_temp^2/besselj(0,kzero_temp)^2*int_temp/D;
    end
    var_vec(zz) = var_temp;
end
var_vec = 2*radius^(2*beta+2)*t_vec.*var_vec;


end


function phi_out = phin_fun(beta, kn)
%xtemp = kn^2/4;
switch beta
    case 2
%    hypergeom1 = 16*(0.5*kn-besselj(1,kn))/kn^3;
    hypergeom1 = 8/kn^2; %besselj(1,kn)=0    
%    hypergeom2 = 64*(-1+kn^2/4 +besselj(0,kn))/kn^4;
    case 4
%    hypergeom1 = 6*(-2*sqrt(xtemp)+xtemp^(3/2)+2*besselj(1,2*sqrt(xtemp)))/xtemp^(5/2);
    hypergeom1 = 24*(kn^2-8)/kn^4; %besselj(1,kn)=besselj(1,2sqrt(x))=0 
%    hypergeom2 = 9*((xtemp-2)^2-4*besselj(0,2*sqrt(xtemp)))/xtemp^3;   
    otherwise %beta=6
%    hypergeom1 = 12*(12*sqrt(xtemp) - 6*sqrt(xtemp^3)+xtemp^(5/2)-12*besselj(1,2*sqrt(xtemp)))/xtemp^(7/2);
    hypergeom1 = 48*(kn^4 - 24*kn^2 + 192)/kn^6; 
%    hypergeom2 = 16*(-36+36*xtemp-9*xtemp^2+xtemp^3+36*besselj(0,2*sqrt(xtemp)))/xtemp^4;
    
end
%phi_out = besselj(0,kn)/(beta+2)*hypergeom1 + kn*besselj(1,kn)/(beta+2)^2*hypergeom2;
phi_out = besselj(0,kn)*hypergeom1/(beta+2); %besselj(1,kn)=0
end