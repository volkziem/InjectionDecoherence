% plot_amplitude_dependence.m
clear; close all
clc
a0=2;
theta2=0.01; % 4*kappa^2*eps^2
f=@(n,a0)a0*exp(-0.5*theta2*a0^2*n.^2./(1+theta2*n.^2))./(1+theta2*n.^2);
%g=@(n,a0)0.5*theta2*a0^2*n.^2./(1+theta2*n.^2);

n=0:50;
plot(sqrt(theta2)*n,f(n,1),'k',sqrt(theta2)*n,f(n,2),'r-.','LineWidth',2)
xlabel('\theta = 2\kappa\epsilon_0n'); ylabel('a_n')
legend('a_0=1','a_0=2')
% hold on; plot(sqrt(theta2)*n,f(n,0.1),'b','LineWidth',2)
set(gca,'Fontsize',16)

