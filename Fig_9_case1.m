% Fig.9, case 1: load 'FF_case1.mat'
% x: equilibrium of the system, obtained from MatCont. The last row
% represents the value of lambda
% y4bar_p: the average of y4, corresponding to lambda_pp. See the way of 
% caculation in line 122-127 of Fig_8.m
% max_lc, min_lc: the maximum and the minimum value of y4 around the limit 
% cycle, corresponding to lambda_p. See the way of calculation in line 
% 122-127 of Fig_8.m
% dy4bar_ana: analytical calculation of the derivative of average y4, in the
% range for limit cycles, corresponding to lambda_p, using iSRC method.
% See line 129-139 of Fig_8.m

clear;
clc;

load('FF_case1.mat')

a = 0.73; b = 0.01; c = 4.78;

eta = @(x) 1/(1+exp((c-x)/a));
dy4 = @(x) 1/(4*(1+eta(x/2)))-(x*exp((2*c-x)/(2*a))*eta(x/2)^2)/(8*a*(1+eta(x/2))^2);

% A smooth curve fitting the data of y4bar_p 
p1 = -0.0001248; p2 = 0.005466; p3 = -0.08849; p4 = 0.6324; p5 = -1.733; p6 = 1.99;
fit_y4bar_p = @(x) p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6;
fit_dy4bar_p = @(x) 5*p1*x^4 + 4*p2*x^3 + 3*p3*x^2 + 2*p4*x + p5; % the derivative of fit_y4bar_p

% direct differentiation of y4bar
sd1 = zeros(1,552);
sd2 = zeros(1,2498);
sd3 = zeros(1,length(lambda_pp));
for i = 1:552
    sd1(i) = (x(7,i+1)-x(7,i))/(x(8,i+1)-x(8,i));
end
for i = 1:2498
    sd2(i) = (x(7,i+1333)-x(7,i+1332))/(x(8,i+1333)-x(8,i+1332));
end
for i = 1:length(lambda_pp)
    sd3(i) = fit_dy4bar_p(lambda_pp(i));
end

% For the ranges of stable equilibria, analytical calculation of the derivative of y4
lambda1 = 5:0.01:7.24;
lambda2 = 12.36:0.01:15;
sa1 = zeros(1,length(lambda1));
sa2 = zeros(1,length(lambda2));
for i = 1:length(lambda1)
    sa1(i) = dy4(lambda1(i));
end
for i = 1:length(lambda2)
    sa2(i) = dy4(lambda2(i));
end

figure(4)

subplot(3,1,1) %the value of y4 at equilibrium and the amplitude of the limit cycle
plot(x(8,1:553),x(7,1:553),'-b','LineWidth',2.5)
hold on
plot(x(8,1333:end),x(7,1333:end),'-b','LineWidth',2.5)
plot(x(8,553:1333),x(7,553:1333),'-r','LineWidth',2.5)
plot(lambda_p,min_lc,'-g','LineWidth',2.5)
plot(lambda_p,max_lc,'-g','LineWidth',2.5)
plot(x(8,553),x(7,553),'.m','MarkerSize',24)
plot(x(8,1333),x(7,1333),'.m','MarkerSize',24)
hold off
ylabel('$y_4$','interpreter','latex')
axis([5 15 1 2])
title('(a) case 1')

subplot(3,1,2)  %the average of y4
plot(x(8,1:553),x(7,1:553),'-b','LineWidth',2.5)
hold on
plot(x(8,1333:end),x(7,1333:end),'-b','LineWidth',2.5)
plot(lambda_pp,y4bar_p,'-r','LineWidth',2.5)
plot(x(8,553),x(7,553),'.m','MarkerSize',24)
plot(x(8,1333),x(7,1333),'.m','MarkerSize',24)
hold off
ylabel('$\overline{y}_4$','interpreter','latex')
axis([5 15 1 2])

subplot(3,1,3)  %the derivative of y4 with respect to lambda
plot(lambda_p,dy4bar_ana,'-r','LineWidth',2.5)
hold on
plot(lambda1,sa1,'-b','LineWidth',2.5)
plot(lambda2,sa2,'-b','LineWidth',2.5)
plot(x(8,2:553),sd1,'--g','LineWidth',1.7)
plot(x(8,1334:end),sd2,'--g','LineWidth',1.7)
plot(lambda_pp,sd3,'--g','LineWidth',1.7)
plot(x(8,553),sd1(end),'.m','MarkerSize',24)
plot(x(8,1333),sd2(1),'.m','MarkerSize',24)
plot(5:0.1:15,zeros(1,length(5:0.1:15)),'--k','LineWidth',1.5)
hold off
xlabel('\lambda')
ylabel('$\partial\overline{y}_4/\partial\lambda$','interpreter','latex')
axis([5 15 -.05 .25])