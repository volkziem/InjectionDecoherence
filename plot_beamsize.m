% plot_beamsize.m, V. Ziemann 220316 
clear; close all
nn=1:5000;     % number of turns
mu=2*pi*0.028; % phase advance per turn  
kappa=0.001;   % amplitude dependent tune shift
eps2=1;        % emittance of injected beam (normalizes all other parameters)
X=[0;1];       % injection steering error
beta=2;        % injection beam Twiss parameters
alpha=0.;
%-------------------------no need to change below
gamma=(1+alpha^2)/beta;
sigma=eps2*[beta,-alpha;-alpha,gamma];   % injected beam matrix

data=zeros(length(nn),9);     % allocate space
[JJ0,~]=Jmrs_sigma(0,mu,kappa,sigma,X);  % constant, taken from the loop

for n=1:length(nn) 
  [~,Xhat]=Jmrs_sigma(n,mu,kappa,sigma,X);      % eq.18
  XX1=real(Xhat); % Xhat1
  XX2=imag(Xhat); % Xhat2
  data(n,1)=XX1; 
  data(n,2)=XX2; 
  JJn=Jmrs_sigma(-2*n,mu,kappa,sigma,X); % the second moments, eq.34
  xx11=0.5*(JJ0(1,1)+real(JJn(1,1)))+imag(JJn(1,2))+0.5*(JJ0(2,2)-real(JJn(2,2)));
  xx12=-0.5*imag(JJn(1,1))+real(JJn(1,2))+0.5*imag(JJn(2,2));
  xx22=0.5*(JJ0(1,1)-real(JJn(1,1)))-imag(JJn(1,2))+0.5*(JJ0(2,2)+real(JJn(2,2)));
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
end

%..................only display below
set(gcf,'Position',[3200,100,1200,800])
subplot(3,1,1)
plot(nn,data(:,1),'k','LineWidth',2); 
xlabel('Number of turns n'); ylabel('X_1'); legend('Position X_1')
set(gca,'Fontsize',16)
subplot(3,1,2)
plot(nn,data(:,6),'k',nn,data(:,7),'r','LineWidth',2);  
ylim([1.1*min(data(:,7)),1.05*max(data(:,6))])
xlabel('Number of turns n'); ylabel('\sigma_{11},\sigma_{12}');
legend('Beam size \sigma_{11}','Correlation \sigma_{12}')
set(gca,'Fontsize',16)
subplot(3,1,3)
plot(nn,data(:,9),'k','LineWidth',2);  
ylim([0.89*min(data(:,9)),1.05*max(data(:,9))])
xlabel('Number of turns n'); ylabel('\epsilon');
legend('Emittance \epsilon','Location','SouthEast')
set(gca,'Fontsize',16)
saved_values_in_data_at_the_end=data(end,:)