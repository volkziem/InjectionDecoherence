% Jmrs_sigma.m, equation for J(m,x_r*x_s) for general beams
% offset X
function [out,Xhat]=Jmrs_sigma(m,mu,kappa,sigma,X)
dd=eye(2)+2i*m*kappa*sigma;
denom=sqrt(det(dd));
ddinv=inv(dd);
psi=-1i*m*kappa*X'*X-2*m^2*kappa^2*X'*ddinv*sigma*X;
Y=ddinv*X;
factor=exp(-1i*m*mu+psi)/denom;        % sign of psi fixed, 220424
out=factor*(ddinv*sigma+Y*Y') ;        % eq.34
Xhat=factor*(Y(1)+1i*Y(2));    % eq.18, if m=n
end
