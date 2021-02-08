% Fig.4: a reference figure showing V, n, m, h for two different values of
% the current, I=50 and I=60

C = 1; 
I1 = 50; I2 = 60;   %two different values of the current
gK = 36; EK = -12;
gNa = 120; ENa = 120;
gL = 0.3; EL = 10.6;

%starting point: intersection of the Poincare section {V=0, dV/dt>0} and
%the limit cycle
initials1 = [0,0.582483641907492,0.048145598184594,0.193952946600929];  
initials2 = [0,0.609153018344056,0.048401876541377,0.156817258793739];

%period
T1 = 8.490726403344505; T2 = 7.982819594507163; 
t0 = 0; dt = 0.01;

[Tu1,Pu1] = ode45(@HH_model,t0:dt:3*T1,initials1,[],C,I1,gK,EK,gNa,ENa,gL,EL);
Vu1 = Pu1(:,1); nu1 = Pu1(:,2); mu1 = Pu1(:,3); hu1 = Pu1(:,4);
[Tu2,Pu2] = ode45(@HH_model,t0:dt:3*T2,initials2,[],C,I2,gK,EK,gNa,ENa,gL,EL);
Vu2 = Pu2(:,1); nu2 = Pu2(:,2); mu2 = Pu2(:,3); hu2 = Pu2(:,4);

figure(1)
%plot for I=50
subplot(2,2,1)
plot(Tu1,Vu1,'-k','LineWidth',2.5)
axis([0 3*T1 -9 82])
ylabel('V (mV)');
title('(a) case 1')
subplot(2,2,3)
plot(Tu1,nu1,'-b','LineWidth',2.5)
hold on
plot(Tu1,mu1,'-r','LineWidth',2.5)
plot(Tu1,hu1,'-m','LineWidth',2.5)
hold off
axis([0 3*T1 -.02 1.02])
xlabel('t'); ylabel('gates')
%plot for I=60
subplot(2,2,2)
plot(Tu2,Vu2,'-k','LineWidth',2.5)
axis([0 3*T2 -9 82])
ylabel('V (mV)');
title('(b) case 2')
subplot(2,2,4)
plot(Tu2,nu2,'-b','LineWidth',2.5)
hold on
plot(Tu2,mu2,'-r','LineWidth',2.5)
plot(Tu2,hu2,'-m','LineWidth',2.5)
hold off
axis([0 3*T2 -.02 1.02])
xlabel('t'); ylabel('gates')