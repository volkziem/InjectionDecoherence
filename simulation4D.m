% simulation4D.m, V. Ziemann 220424 
clear; close all
degree=pi/180;
nn=1:300;        % number of turns
mux=2*pi*0.028;  % horizontal phase advance per turn  
muy=2*pi*0.041;  % vertical phase advance per turn  
kappaxx=0.001;   % amplitude dependent tune shift 
kappaxy=0.001/2;
kappayy=0.002; 
kappax=diag([kappaxx,kappaxx,kappaxy,kappaxy]);
kappay=diag([kappaxy,kappaxy,kappayy,kappayy]);
X=[1;0;0;0];     % injection steering error
%.....Twiss parameters at injection
epsx=10;                                        % horizontal 
betax=3; alphax=0; gammax=(1+alphax^2)/betax;   
sigmax=epsx*[betax,-alphax;-alphax,gammax];
epsy=1;                                         % vertical 
betay=3; alphay=0; gammay=(1+alphay^2)/betay;   
sigmay=epsy*[betay,-alphay;-alphay,gammay];
eta=30*degree;                                  % coupling, roll angle
C=-sin(eta)*eye(2); g=sqrt(1-det(C)); Cplus=inv(C)*det(C);
Tinv=[g*eye(2),C;-Cplus,g*eye(2)];                    % eq.52
sigma=Tinv*[sigmax,zeros(2);zeros(2),sigmay]*Tinv';   % eq.51
epsilon0=diag([epsx,epsx,epsy,epsy]);

%-------------------------no need to change below
Bmag=Bmags(sigmax,eye(2))
Cmag=Bmags(sigmay,eye(2))     % eq.59
emittance_growth=[
  epsx*Bmag*cos(eta)^2+epsy*sin(eta)^2*Cmag+0.5*(X(1)^2+X(2)^2)
  epsy*Cmag*cos(eta)^2+epsx*sin(eta)^2*Bmag+0.5*(X(3)^2+X(4)^2)]

data=zeros(length(nn),10);     % allocate space
datay=data;
[JJ0,~]=Jmrs_sigma4D(0,mux,kappax,sigma,X);   % constant, taken from the loop
[JJ0y,~]=Jmrs_sigma4D(0,muy,kappay,sigma,X);  % vertical
for n=1:length(nn) 
  [~,Xhat,~]=Jmrs_sigma4D(n,mux,kappax,sigma,X);      % eq.18 
  XX1=real(Xhat); % Xhat1
  XX2=imag(Xhat); % Xhat2
  data(n,1)=XX1; 
  data(n,2)=XX2; 
  JJn=Jmrs_sigma4D(2*n,mux,kappax,sigma,X); % the second moments, eq.34, sign flippped
%   xx11=0.5*(JJ0(1,1)+real(JJn(1,1)))+imag(JJn(1,2))+0.5*(JJ0(2,2)-real(JJn(2,2)));
%   xx12=-0.5*imag(JJn(1,1))+real(JJn(1,2))+0.5*imag(JJn(2,2));
%   xx22=0.5*(JJ0(1,1)-real(JJn(1,1)))-imag(JJn(1,2))+0.5*(JJ0(2,2)+real(JJn(2,2)));
  xx11=0.5*(JJ0(1,1)+real(JJn(1,1)))-imag(JJn(1,2))+0.5*(JJ0(2,2)-real(JJn(2,2)));
  xx12=0.5*imag(JJn(1,1))+real(JJn(1,2))-0.5*imag(JJn(2,2));
  xx22=0.5*(JJ0(1,1)-real(JJn(1,1)))+imag(JJn(1,2))+0.5*(JJ0(2,2)+real(JJn(2,2)));
  data(n,3)=xx11;      % second moments
  data(n,4)=xx12;
  data(n,5)=xx22;
  sig11=xx11-XX1^2;    % centroid subtracted
  sig12=xx12-XX1*XX2;
  sig22=xx22-XX2^2;
  emit=sqrt(sig11*sig22-sig12^2);  % emittance
  data(n,6)=sig11;
  data(n,7)=sig12;
  data(n,8)=sig22;
  data(n,9)=emit;  
  %........vertical plane
  [~,~,Yhat]=Jmrs_sigma4D(n,muy,kappay,sigma,X);
  YY1=real(Yhat); % Xhat1
  YY2=imag(Yhat); % Xhat2
  datay(n,1)=YY1;
  datay(n,2)=YY2;
  JJn=Jmrs_sigma4D(2*n,muy,kappay,sigma,X); % the second moments, eq.34, sign flippped
  xx33=0.5*(JJ0y(3,3)+real(JJn(3,3)))-imag(JJn(3,4))+0.5*(JJ0y(4,4)-real(JJn(4,4)));
  xx34=0.5*imag(JJn(3,3))+real(JJn(3,4))-0.5*imag(JJn(4,4));
  xx44=0.5*(JJ0y(3,3)-real(JJn(3,3)))+imag(JJn(3,4))+0.5*(JJ0y(4,4)+real(JJn(4,4)));
  datay(n,3)=xx33;      % second moments
  datay(n,4)=xx34;
  datay(n,5)=xx44;
  sig33=xx33-YY1^2;    % centroid subtracted
  sig34=xx34-YY1*YY2;
  sig44=xx44-YY2^2;
  emity=sqrt(sig33*sig44-sig34^2);  % vertical emittance
  datay(n,6)=sig33;
  datay(n,7)=sig34;
  datay(n,8)=sig44;
  datay(n,9)=emity; 
  %......................correlation sigma13
  % I use +n and flip the signs of all imaginary parts, compared to the
  % appendix, because then we end up on the same branch of the square root
  JJp=Jmrs_sigma4D(n,mux+muy,kappax+kappay,sigma,X); % flipped sign of n
  JJm=Jmrs_sigma4D(n,mux-muy,kappax-kappay,sigma,X);
  xx13=real(JJp(1,3))+real(JJm(1,3))-imag(JJp(1,4))+imag(JJm(1,4)) ...
    -imag(JJp(2,3))-imag(JJm(2,3))-real(JJp(2,4))+real(JJm(2,4));
  xx13=0.5*xx13;
  sig13=xx13-XX1*YY1;
  datay(n,10)=sig13;
end

%..................only display below
set(gcf,'Position',[3200,100,1200,800])
subplot(3,1,1)
plot(nn,data(:,1),'k',nn,datay(:,1),'r-.','LineWidth',2); 
xlabel('Number of turns n'); 
ylabel('$\hat X_1, \hat X_3$','interpreter','latex'); 
legend('Horizontal position $\hat X_1$','Vertical position $\hat X_3$','interpreter','latex'); 
set(gca,'Fontsize',16)
subplot(3,1,2)
plot(nn,data(:,6),'k',nn,datay(:,6),'r',nn,datay(:,10),'b--','LineWidth',2);  
xlabel('Number of turns n'); 
ylabel('$\hat\sigma_{11}, \hat\sigma_{33}, \hat\sigma_{13}$','interpreter','latex');
legend('Horizontal beam size $\hat\sigma_{11}$','Vertical beam size $\hat\sigma_{33}$',...
  'Correlation $\hat\sigma_{13}$','interpreter','latex');
set(gca,'Fontsize',16)
subplot(3,1,3)
plot(nn,data(:,9),'k',nn,datay(:,9),'r','LineWidth',2);  
ylim([0.0*min(data(:,9)),1.05*max(data(:,9))])
xlabel('Number of turns n'); ylabel('$\hat\epsilon_x, \hat\epsilon_y$','interpreter','latex');
legend('Horizontal emittance $\hat\epsilon_x$','Vertical emittance $\hat\epsilon_y$','Location','SouthEast','interpreter','latex')
set(gca,'Fontsize',16)
saved_values_in_data_at_the_end=data(end,:)

%..................................Jmrs_sigma4D
% Jmrs_sigma.m, equation for J(m,x_r*x_s) for general beams
% offset X, kappa is a diagonal matrix with the same size of sigma
function [out,Xhat,Yhat]=Jmrs_sigma4D(m,mu,kappa,sigma,X)
dd=eye(size(sigma))+2i*m*sigma*kappa;
denom=sqrt(det(dd));
if phase(denom)<0,denom=-denom; end    % ensure consistent complex branch 
ddinv=inv(dd);
psi=-1i*m*X'*kappa*X-2*m^2*X'*kappa*ddinv*sigma*kappa*X;
Y=ddinv*X;
factor=exp(-1i*m*mu+psi)/denom;        % sign of psi fixed, 220424
out=factor*(ddinv*sigma+Y*Y');         % eq.34
Xhat=factor*(Y(1)+1i*Y(2));            % eq.18, if m=n
Yhat=factor*(Y(3)+1i*Y(4));
end

%............................................Bmags
function out=Bmags(sig1,sig2)
eps1=sqrt(det(sig1)); beta1=sig1(1,1)/eps1; alpha1=-sig1(1,2)/eps1;
eps2=sqrt(det(sig2)); beta2=sig2(1,1)/eps2; alpha2=-sig2(1,2)/eps2;
out=0.5*(beta1/beta2+beta2/beta1+beta1*beta2*(alpha1/beta1-alpha2/beta2)^2);
end