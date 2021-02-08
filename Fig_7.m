% Fig.7: Empirical curves with analytically derived tangent curves for V,
% sodium current and potassium current
% Need to first run Fig.5 to find iSRC eta(t)/xi(t)

I = 50;
initials = [0,0.582483641907492,0.048145598184594,0.193952946600929];
t0 = 0; tF = 100; dt = 0.01; 

%% Compute the average of the quantities for I over [45,55]
eps = -5:0.2:5;
Vbar = zeros(1,length(eps));
I_Na_bar = zeros(1,length(eps));
I_K_bar = zeros(1,length(eps));

for k = 1:length(eps)
Ip = I+eps(k);
[T1,P1] = ode45(@HH_model,t0:dt:tF,initials,[],C,Ip,gK,EK,gNa,ENa,gL,EL);
initials_p=P1(end,:);
options = odeset('Events',@HH_event1);
[~,~,~,initials_p,~] = ode45(@HH_model,[0 tF],initials_p,options,C,Ip,gK,EK,gNa,ENa,gL,EL);
[~,~,Teps,~,~] = ode45(@HH_model,[0 tF],initials_p,options,C,Ip,gK,EK,gNa,ENa,gL,EL);
Teps = Teps(end);  %period
dteps = (Teps/T0)*dt; tspan = t0:dteps:Teps;
[T2,P2] = ode45(@HH_model,tspan,initials_p,[],C,Ip,gK,EK,gNa,ENa,gL,EL);
Vp = P2(:,1); np = P2(:,2); mp = P2(:,3); hp = P2(:,4);

%average of voltage
Vbar(k) = (1/Teps)*(trapz(Vp)*dteps); 

%average of sodium current
I_Na = @(V,n,m,h) -gNa*m^3*h*(V-ENa);
I_Na_p = zeros(length(T2),1); 
for i = 1:length(T2)
    I_Na_p(i) = I_Na(Vp(i),np(i),mp(i),hp(i));
end
I_Na_bar(k) = (1/Teps)*(trapz(I_Na_p)*dteps);

%average of potassium current
I_K = @(V,n,m,h) -gK*n^4*(V-EK);
I_K_p = zeros(length(T2),1); 
for i = 1:length(T2)
    I_K_p(i) = I_K(Vp(i),np(i),mp(i),hp(i));
end
I_K_bar(k) = (1/Teps)*(trapz(I_K_p)*dteps);

end

%% Compute the derivative of the quantities at I = 50 using iSRC eta(t)
T0 = 8.490726403344505;
[Tu,Pu] = ode45(@HH_model,t0:dt:T0,initials,[],C,I,gK,EK,gNa,ENa,gL,EL);
Vu = Pu(:,1); nu = Pu(:,2); mu = Pu(:,3); hu = Pu(:,4);

% Voltage
JV = zeros(length(Tu),4);
inte_V = zeros(length(Tu),1);
for i = 1:length(Tu)
    JV(i,:) = [1 0 0 0]; %Jacobian of V 
    inte_V(i) = JV(i,:)*eta(:,i);
end
dVbar = (1/T0)*(dt*trapz(inte_V));  %derivative at I=50

Ip = 48:0.02:52;
V_tan = Vbar(26)*ones(1,length(Ip))+dVbar*(Ip-I*ones(1,length(Ip))); %tangent curve at I=50

% Sodium current
JJNa = @(V,n,m,h) [-gNa*m^3*h 0 -3*gNa*m^2*h*(V-ENa) -gNa*m^3*(V-ENa)];  %Jacobian of I_Na
JNa = zeros(length(Tu),4);
inte_Na = zeros(length(Tu),1);
for i = 1:length(Tu)
    JNa(i,:) = JJNa(Vu(i),nu(i),mu(i),hu(i));
    inte_Na(i) = JNa(i,:)*eta(:,i);
end
dI_Na_bar = (1/T0)*(dt*trapz(inte_Na)); %derivative at I=50

I_Na_tan = I_Na_bar(26)*ones(1,length(Ip))+dI_Na_bar*(Ip-I*ones(1,length(Ip))); %tangent curve at I=50

% Potassium current
JJK = @(V,n,m,h) [-gK*n^4 -4*gK*n^3*(V-EK) 0 0];  %Jacobian of I_K
JK = zeros(length(Tu),4);
inte_K = zeros(length(Tu),1);
for i = 1:length(Tu)
    JK(i,:) = JJK(Vu(i),nu(i),mu(i),hu(i));
    inte_K(i) = JK(i,:)*eta(:,i);
end
dI_K_bar = (1/T0)*(dt*trapz(inte_K));  %derivative at I=50

I_K_tan = I_K_bar(26)*ones(1,length(Ip))+dI_K_bar*(Ip-I*ones(1,length(Ip))); %tangent curve at I=50

%% Plot
figure(1)
subplot(1,3,1) %V
plot(I*ones(1,length(eps))+eps,Vbar,'k*','MarkerSize',4)
hold on
plot(I,Vbar(26),'k.','MarkerSize',23)
plot(Ip,V_tan,'-r','LineWidth',2.5)
hold off
axis([45 55 15.6 16.8])
xlabel('I'); ylabel('$\overline{V}$','Interpreter','latex');
subplot(1,3,2) %Sodium current
plot(I*ones(1,length(eps))+eps,I_Na_bar,'k*','MarkerSize',4)
hold on
plot(I,I_Na_bar(26),'k.','MarkerSize',23)
plot(Ip,I_Na_tan,'-r','LineWidth',2.5)
hold off
axis([45 55 107 109])
xlabel('I'); ylabel('$\overline{I}_{Na}$','Interpreter','latex');
subplot(1,3,3) %Potassium current
plot(I*ones(1,length(eps))+eps,I_K_bar,'k*','MarkerSize',4)
hold on
plot(Ip,I_K_tan,'-r','LineWidth',2.5)
plot(I,I_K_bar(26),'k.','MarkerSize',23)
hold off
axis([45 55 -161 -151.5])
xlabel('I'); ylabel('$\overline{I}_{K}$','Interpreter','latex');