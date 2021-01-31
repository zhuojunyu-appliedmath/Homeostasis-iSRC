% Fig.4: iSRC components for four variables of the Hodgkin-Huxley system, 
% for two different Poincare sections / initial conditions.

C = 1; 
eps = 0.5; 
I = 50; Ip = I+eps;
gK = 36; EK = -12;
gNa = 120; ENa = 120;
gL = 0.3; EL = 10.6;

alpha_n = @(V) 0.01*(10-V)/(exp((10-V)/10)-1);
alpha_m = @(V) 0.1*(25-V)/(exp((25-V)/10)-1);
alpha_h = @(V) 0.07*exp(-V/20);
beta_n = @(V) 0.125*exp(-V/80);
beta_m = @(V) 4*exp(-V/18);
beta_h = @(V) 1/(exp((30-V)/10)+1);

dalpha_n = @(V) 0.01*(1-0.1*V*exp((10-V)/10))/(exp((10-V)/10)-1)^2;
dalpha_m = @(V) 0.1*((1.5-0.1*V)*exp((25-V)/10)+1)/(exp((25-V)/10)-1)^2;
dalpha_h = @(V) -(0.07/20)*exp(-V/20);
dbeta_n = @(V) -(0.125/80)*exp(-V/80);
dbeta_m = @(V) -(4/18)*exp(-V/18);
dbeta_h = @(V) 0.1*exp((30-V)/10)/(exp((30-V)/10)+1)^2;

F = @(V,n,m,h) [(1/C)*(I-gK*n^4*(V-EK)-gNa*m^3*h*(V-ENa)-gL*(V-EL));
                alpha_n(V)*(1-n)-beta_n(V)*n;
                alpha_m(V)*(1-m)-beta_m(V)*m;
                alpha_h(V)*(1-h)-beta_h(V)*h];
            
A = @(V,n,m,h) [(1/C)*(-gK*n^4-gNa*m^3*h-gL) (1/C)*(-4*gK*n^3*(V-EK)) (1/C)*(-3*gNa*m^2*h*(V-ENa)) (1/C)*(-gNa*m^3*(V-ENa));
                (1-n)*dalpha_n(V)-n*dbeta_n(V) -alpha_n(V)-beta_n(V) 0 0;
                (1-m)*dalpha_m(V)-m*dbeta_m(V) 0 -alpha_m(V)-beta_m(V) 0;
                (1-h)*dalpha_h(V)-h*dbeta_h(V) 0 0 -alpha_h(V)-beta_h(V)];

initials = [0,0.582483641907492,0.048145598184594,0.193952946600929];
t0 = 0; dt = 0.001; tF = 100; T0 = 8.490726403344505;

%% Monodromy matrix
initials_funda = [initials, 1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1];
[~,P] = ode45(@HH_fundamental,t0:dt:T0,initials_funda,[],C,I,gK,EK,gNa,ENa,gL,EL);
U = P(end,5:20);
M = [U(1:4); U(5:8); U(9:12); U(13:16)]';  %monodromy matrix

[v d] = eigs(M);  %The monodromy matrix has an eigenvalue equal to 1
z0 = v(:,1);  %The initial condition of the iPRC: the eigevector of M associated with unit eigenvalue

%% iPRC: backward integration
tspan = t0:dt:7*T0;
[~,P] = ode45(@HH_model,tspan,initials,[],C,I,gK,EK,gNa,ENa,gL,EL);
V = P(:,1); n = P(:,2); m = P(:,3); h = P(:,4); 
FF = zeros(4,length(tspan)); DF = zeros(4,4,length(tspan));
for i = 1:length(tspan)
    FF(:,i) = F(V(i),n(i),m(i),h(i));
    DF(:,:,i) = A(V(i),n(i),m(i),h(i));
end

z = zeros(4,length(tspan));  %iPRC
k0 = FF(:,end)'*z0;  %normalization
z(:,length(tspan)) = z0/k0;
for i = length(tspan)-1:-1:1
    z(:,i) = z(:,i+1)+DF(:,:,i+1)'*z(:,i+1)*dt;
    k = FF(:,i)'*z(:,i);
    z(:,i) = z(:,i)/k;
end

dFp = [1/C; 0; 0; 0];
int = zeros(8490,1);
for i = 1:length(int)
    int(i) = z(:,5*8490+i)'*dFp;
end
T1_prc = -dt*trapz(int); %linear shift of the period

%% iSRC: inhomogeneous variational equation
nu1 = T1_prc/T0;  %timing sensitivity

% find the point of the perturbed limit cycle transverse at section {V=0,dV/dt>0}
[Tp,Pp] = ode45(@HH_model,t0:dt:100,initials,[],C,Ip,gK,EK,gNa,ENa,gL,EL);
initials_p = Pp(end,:);
options = odeset('Events',@HH_event1);
[~,~,~,initials_p1,~] = ode45(@HH_model,[0 tF],initials_p,options,C,Ip,gK,EK,gNa,ENa,gL,EL);

eta0 = (initials_p1-initials)/eps;  %initial condition for iSRC eta(t)
initials_iSRC1 = [initials,eta0];
[~,G1] = ode45(@HH_SRC,0:0.01:3*T0,initials_iSRC1,[],C,I,gK,EK,gNa,ENa,gL,EL,nu1);
eta = G1(:,5:8)'; 

% find the point of the perturbed limit cycle transverse at section {n=0.582483641907492,dn/dt<0}
options = odeset('Events',@HH_event2);
[~,~,~,initials_p2,~] = ode45(@HH_model,[0 tF],initials_p,options,C,Ip,gK,EK,gNa,ENa,gL,EL);

xi0 = (initials_p2-initials)/eps;  %initial condition for iSRC xi(t)
initials_iSRC2 = [initials,xi0];
[T,G2] = ode45(@HH_SRC,0:0.01:3*T0,initials_iSRC2,[],C,I,gK,EK,gNa,ENa,gL,EL,nu1);
xi = G2(:,5:8)'; 

figure(1)
% for section {V=0,dV/dt>0}
subplot(2,2,1)
plot(T,eta(1,:),'-k','LineWidth',2.5)
axis([0 3*T0 -4.8 2.6])
ylabel('\eta_{V}')
title('(a) \eta(t)')
subplot(2,2,3)
plot(T,eta(2,:),'-b','LineWidth',2.5)
hold on
plot(T,eta(3,:),'-r','LineWidth',2.5)
plot(T,eta(4,:),'-m','LineWidth',2.5)
hold off
xlabel('t'); ylabel('\eta_{gates}')
axis([0 3*T0 -.044 .054])
% for section {n=0.582483641907492,dn/dt<0}
subplot(2,2,2)
plot(T,xi(1,:),'-k','LineWidth',2.5)
ylabel('\xi_{V}')
axis([0 3*T0 -1.7 3])
title('(b) \xi(t)')
subplot(2,2,4)
plot(T,xi(2,:),'-b','LineWidth',2.5)
hold on
plot(T,xi(3,:),'-r','LineWidth',2.5)
plot(T,xi(4,:),'-m','LineWidth',2.5)
hold off
xlabel('t'); ylabel('\xi_{gates}')
axis([0 3*T0 -.021 .0315])