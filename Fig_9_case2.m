% Fig.9, case 2: load 'FF_case2.mat'
% x: equilibrium of the system, obtained from MatCont. The last row
% represents the value of lambda
% y4bar_p(1): the average of y4, corresponding to lambda_ppp(1). See the 
% way of caculation in line 112-117 of Fig_8.m
% max_lc(1), min_lc(1): the maximum and the minimum value of y4 around the 
% limit cycle, corresponding to lambda_pp(1). See the way of calculation in 
% line 112-116 of Fig_8.m
% dy4bar_ana(1): analytical calculation of the derivative of average y4 in
% the range for limit cycles, corresponding to lambda_p(1), using iSRC method.
% See line 119-129 of Fig_8.m

clear;
clc;

load('FF_case2.mat')

a = 0.65; b = 0.01; c = 4.88;

eta = @(x) 1/(1+exp((c-x)/a));
dy4 = @(x) 1/(4*(1+eta(x/2)))-(x*exp((2*c-x)/(2*a))*eta(x/2)^2)/(8*a*(1+eta(x/2))^2);

% A smooth curve fitting the data of y4bar_p 
p1 = 0.0007896; p2 = -0.02565; p3 = 0.2663; p4 = -0.9202; p5 = 1.873;
fit_y4bar_p1 = @(x) p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5;
fit_dy4bar_p1 = @(x) 4*p1*x^3 + 3*p2*x^2 + 2*p3*x + p4;
% A smooth curve fitting the data of y4bar_p1
q1 = -0.001432; q2 = 0.06438; q3 = -1.06; q4 = 7.585; q5 = -18.26;
fit_y4bar_p2 = @(x) q1*x^4 + q2*x^3 + q3*x^2 + q4*x + q5;
fit_dy4bar_p2 = @(x) 4*q1*x^3 + 3*q2*x^2 + 2*q3*x + q4;

% direct differentiation of y4bar
sd1 = zeros(1,538);
sd2 = zeros(1,233);
sd3 = zeros(1,2433);
sd4 = zeros(1,length(lambda_ppp));
sd5 = zeros(1,length(lambda_ppp1));
for i = 1:538
    sd1(i) = (x(7,i+1)-x(7,i))/(x(8,i+1)-x(8,i));
end
for i = 1:233
    sd2(i) = (x(7,i+1054)-x(7,i+1053))/(x(8,i+1054)-x(8,i+1053));
end
for i = 1:2433
    sd3(i) = (x(7,i+1930)-x(7,i+1929))/(x(8,i+1930)-x(8,i+1929));
end
for i = 1:length(lambda_ppp)
    sd4(i) = fit_dy4bar_p1(lambda_ppp(i));
end
for i = 1:length(lambda_ppp1)
    sd5(i) = fit_dy4bar_p2(lambda_ppp1(i));
end

% For the ranges of stable equilibria, analytical calculation of the derivative of y4
lambda1 = 5:0.01:6.74;
lambda2 = 7.8:0.01:9.18;
lambda3 = 12.48:0.01:15;
sa1 = zeros(1,length(lambda1));
sa2 = zeros(1,length(lambda2));
sa3 = zeros(1,length(lambda3));
for i = 1:length(lambda1)
    sa1(i) = dy4(lambda1(i));
end
for i = 1:length(lambda2)
    sa2(i) = dy4(lambda2(i));
end
for i = 1:length(lambda3)
    sa3(i) = dy4(lambda3(i));
end

figure(4)

subplot(3,1,1)  %the value of y4 at equilibrium and the amplitude of the limit cycle
plot(x(8,1:539),x(7,1:539),'-b','LineWidth',2.5)
hold on
plot(x(8,1053:1287),x(7,1053:1287),'-b','LineWidth',2.5)
plot(x(8,1930:end),x(7,1930:end),'-b','LineWidth',2.5)
plot(x(8,539:1053),x(7,539:1053),'-r','LineWidth',2.5)
plot(x(8,1287:1930),x(7,1287:1930),'-r','LineWidth',2.5)
plot(lambda_pp,min_lc,'-g','LineWidth',2.5)
plot(lambda_pp,max_lc,'-g','LineWidth',2.5)
plot(lambda_pp1,min_lc1,'-g','LineWidth',2.5)
plot(lambda_pp1,max_lc1,'-g','LineWidth',2.5)
plot(x(8,539),x(7,539),'.m','MarkerSize',24)
plot(x(8,1053),x(7,1053),'.m','MarkerSize',24)
plot(x(8,1287),x(7,1287),'.m','MarkerSize',24)
plot(x(8,1930),x(7,1930),'.m','MarkerSize',24)
hold off
axis([5 15 1 2])
ylabel('$y_4$','interpreter','latex')
title('(b) case 2')

subplot(3,1,2)  %the average of y4
plot(x(8,1:539),x(7,1:539),'-b','LineWidth',2.5)
hold on
plot(x(8,1053:1287),x(7,1053:1287),'-b','LineWidth',2.5)
plot(x(8,1930:end),x(7,1930:end),'-b','LineWidth',2.5)
plot(lambda_ppp,y4bar_p,'-r','LineWidth',2.5)
plot(lambda_ppp1,y4bar_p1,'-r','LineWidth',2.5)
plot(x(8,539),x(7,539),'.m','MarkerSize',24)
plot(x(8,1053),x(7,1053),'.m','MarkerSize',24)
plot(x(8,1287),x(7,1287),'.m','MarkerSize',24)
plot(x(8,1930),x(7,1930),'.m','MarkerSize',24)
hold off
axis([5 15 1 2])
ylabel('$\overline{y}_4$','interpreter','latex')

subplot(3,1,3)  %the derivative of y4 with respect to lambda
plot(lambda_p,dy4bar_ana,'-r','LineWidth',2.5)
hold on
plot(lambda_p1,dy4bar_ana1,'-r','LineWidth',2.5)
plot(lambda1,sa1,'-b','LineWidth',2.5)
plot(lambda2,sa2,'-b','LineWidth',2.5)
plot(lambda3,sa3,'-b','LineWidth',2.5)
plot(x(8,2:539),sd1,'--g','LineWidth',1.7)
plot(x(8,1055:1287),sd2,'--g','LineWidth',1.7)
plot(x(8,1931:end),sd3,'--g','LineWidth',1.7)
plot(lambda_ppp,sd4,'--g','LineWidth',1.7)
plot(lambda_ppp1,sd5,'--g','LineWidth',1.7)
plot(x(8,553),sd1(end),'.m','MarkerSize',24)
plot(x(8,1053),sd2(1),'.m','MarkerSize',24)
plot(x(8,1287),sd2(end),'.m','MarkerSize',24)
plot(x(8,1930),sd3(1),'.m','MarkerSize',24)
plot(5:0.1:15,zeros(1,length(5:0.1:15)),'--k','LineWidth',1.5)
hold off
axis([5 15 -.08 0.25])
xlabel('\lambda')
ylabel('$\partial\overline{y}_4/\partial\lambda$','interpreter','latex')