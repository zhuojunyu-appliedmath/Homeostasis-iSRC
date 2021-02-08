% Fig.6: Compare the difference of two distinct iSRCs with the unperturbed flow 
% Need to first run Fig.5 to find the iSRCs

[Tu,Pu] = ode45(@HH_model,0:0.01:3*T0,initials,[],C,I,gK,EK,gNa,ENa,gL,EL);
Vu = Pu(:,1); nu = Pu(:,2); mu = Pu(:,3); hu = Pu(:,4);
FF = zeros(4,length(Tu));
Delta = zeros(4,length(Tu));
for i = 1:length(Tu)
    FF(:,i) = F(Vu(i),nu(i),mu(i),hu(i));
    Delta(:,i) = eta(:,i)-xi(:,i);
end

phi = -0.0585;
figure(1)
subplot(2,2,1) %V
plot(Tu,Delta(1,:),'-k','LineWidth',2.5)
hold on
plot(Tu,phi*FF(1,:),'--r','LineWidth',2.5)
hold off
ylabel('V (mV)')
axis([0 3*T0 -7.3 3.6])
subplot(2,2,2) %n
plot(Tu,Delta(2,:),'-k','LineWidth',2.5)
hold on
plot(Tu,phi*FF(2,:),'--r','LineWidth',2.5)
hold off
ylabel('n')
axis([0 3*T0 -.0142 .0055])
subplot(2,2,3) %m
plot(Tu,Delta(3,:),'-k','LineWidth',2.5)
hold on
plot(Tu,phi*FF(3,:),'--r','LineWidth',2.5)
hold off
xlabel('t'); ylabel('m')
axis([0 3*T0 -.071 .071])
subplot(2,2,4) %h
plot(Tu,Delta(4,:),'-k','LineWidth',2.5)
hold on
plot(Tu,phi*FF(4,:),'--r','LineWidth',2.5)
hold off
xlabel('t'); ylabel('h')
axis([0 3*T0 -.0052 .01])