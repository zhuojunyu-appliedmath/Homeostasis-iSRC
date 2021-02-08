% Fig.9: A reference figure showing the time evolution of y4 for two casees
% and their respective iSRC in y4 direction

clear
clc

eps = 0.05;
t0 = 0; dt = 0.01; tF = 5000;

eta = @(x1,a,c) 1/(1+exp((c-x1)/a));
zeta = @(y4,b) 10/(1+y4^10)+b;

% Flow
F = @(x1,x2,x3,y1,y2,y3,y4,a,b,c,lambda) [lambda-2*x1; x1-2*x2; x2-(1+eta(x1,a,c))*x3; 
                                          x3-zeta(y4,b)*y1; zeta(y4,b)*y1-y2; y2-y3; y3-y4];

%Jacobian of F
A = @(x1,x3,y1,y4,a,b,c) [-2 0 0 0 0 0 0;
                          1 -2 0 0 0 0 0;
                          -(x3*exp((c-x1)/a))/(a*(1+exp((c-x1)/a))^2) 1 -(1+eta(x1,a,c)) 0 0 0 0;
                          0 0 1 -zeta(y4,b) 0 0 (100*y1*y4^9)/(1+y4^10)^2;
                          0 0 0 zeta(y4,b) -1 0 -(100*y1*y4^9)/(1+y4^10)^2;
                          0 0 0 0 1 -1 0;
                          0 0 0 0 0 1 -1];  

dFp = [1; 0; 0; 0; 0; 0; 0];  %the derivative of F w.r.t lambda

% tentative initial condition
initials = [3.4,1.7,1.555277820576277,12.041807913178383,1.394235483634850,1.494685201593891,1.564220169105648]; 

%% case 1:
a1 = 0.73; b1 = 0.01; c1 = 4.78;
lambda1 = 10.5;

% find the exact intersection point on the Poincare section {y4=1.6,dy4/dt>0}
% and the limit cycle as the initial point for further running
[~,PP1] = ode45(@FF_model,[t0 tF],initials,[],lambda1,a1,b1,c1);
options = odeset('Events',@FF_event);
[~,~,~,initials1,~] = ode45(@FF_model,[0 tF],PP1(end,:),options,lambda1,a1,b1,c1);  %initial1 is the initial point
[~,~,T1,~,~] = ode45(@FF_model,[0 tF],initials1,options,lambda1,a1,b1,c1);
T1 = T1(end);  %period

[Tu1,Pu1] = ode45(@FF_model,t0:dt:3*T1,initials1,[],lambda1,a1,b1,c1);  %run for three periods
y41 = Pu1(:,7); 

% find the monodromy matrix
initials_funda1 = [initials1,...
                   1 0 0 0 0 0 0, 0 1 0 0 0 0 0, 0 0 1 0 0 0 0,...
                   0 0 0 1 0 0 0, 0 0 0 0 1 0 0, 0 0 0 0 0 1 0, 0 0 0 0 0 0 1];

[~,P1] = ode45(@FF_fundamental,t0:dt:T1,initials_funda1,[],lambda1,a1,b1,c1);
U1 = P1(end,8:56);

M1 = [U1(1:7); U1(8:14); U1(15:21); U1(22:28); U1(29:35); U1(36:42); U1(43:49)]';  %monodromy matrix
[v1 d1] = eig(M1);  %The monodromy matrix has an eigenvalue equal to 1

% iPRC
tspan_z1 = t0:dt:100*T1;  %first run 100 periods
[~,Q1] = ode45(@FF_model,tspan_z1,initials1,[],lambda1,a1,b1,c1);
xq11 = Q1(:,1); xq21 = Q1(:,2); xq31 = Q1(:,3); 
yq11 = Q1(:,4); yq21 = Q1(:,5); yq31 = Q1(:,6); yq41 = Q1(:,7); 
FFq1 = zeros(7,length(tspan_z1)); DFq1 = zeros(7,7,length(tspan_z1));
for i = 1:length(tspan_z1)
    FFq1(:,i) = F(xq11(i),xq21(i),xq31(i),yq11(i),yq21(i),yq31(i),yq41(i),a1,b1,c1,lambda1);
    DFq1(:,:,i) = A(xq11(i),xq31(i),yq11(i),yq41(i),a1,b1,c1);
end

z1 = zeros(7,length(tspan_z1));
z10 = real(v1(:,4)); %The initial condition of the iPRC: the eigevector of M associated with unit eigenvalue
k0 = FFq1(:,end)'*z10;
z1(:,length(tspan_z1)) = z10/k0;
for i = length(tspan_z1)-1:-1:1  %backward integration
    z1(:,i) = z1(:,i+1)+DFq1(:,:,i+1)'*z1(:,i+1)*dt;
    k = FFq1(:,i)'*z1(:,i);
    z1(:,i) = z1(:,i)/k;
end

zz10 = z1(:,1); %take the last point of last run as the initial condition for z
tspan1 = t0:dt:T1;  %run for one period
[~,QQ1] = ode45(@FF_model,tspan1,initials1,[],lambda1,a1,b1,c1);
xqq11 = QQ1(:,1); xqq21 = QQ1(:,2); xqq31 = QQ1(:,3); 
yqq11 = QQ1(:,4); yqq21 = QQ1(:,5); yqq31 = QQ1(:,6); yqq41 = QQ1(:,7); 
FF1 = zeros(7,length(tspan1)); DF1 = zeros(7,7,length(tspan1));
for i = 1:length(tspan1)
    FF1(:,i) = F(xqq11(i),xqq21(i),xqq31(i),yqq11(i),yqq21(i),yqq31(i),yqq41(i),a1,b1,c1,lambda1);
    DF1(:,:,i) = A(xqq11(i),xqq31(i),yqq11(i),yqq41(i),a1,b1,c1);
end

zz1 = zeros(7,length(tspan1));
k0 = FF1(:,end)'*zz10;
zz1(:,length(tspan1))=zz10/k0;
for i = length(tspan1)-1:-1:1
    zz1(:,i) = zz1(:,i+1)+DF1(:,:,i+1)'*zz1(:,i+1)*dt;
    k = FF1(:,i)'*zz1(:,i);
    zz1(:,i) = zz1(:,i)/k;
end

int1 = zeros(length(tspan1),1);
for i = 1:length(tspan1)
    int1(i) = zz1(:,i)'*dFp;
end
T1_prc1 = -dt*trapz(int1);

% iSRC
nu1 = T1_prc1/T1;

% find the perturbed point on the Poincare section and the perturbed limit
% cycle
lambda_p1 = lambda1+eps;
[~,Pp1] = ode45(@FF_model,[t0 tF],initials1,[],lambda_p1,a1,b1,c1);
options = odeset('Events',@FF_event);
[~,~,~,initials_p1,~] = ode45(@FF_model,[t0 tF],Pp1(end,:),options,lambda_p1,a1,b1,c1);
[~,~,Teps1,~,~] = ode45(@FF_model,[t0 tF],initials_p1,options,lambda_p1,a1,b1,c1);
Teps1 = Teps1(end);  %perturbed period

eta10 = (initials_p1-initials1)/eps;  
initials_SRC1 = [initials1,eta10];

[~,G1] = ode45(@FF_SRC,t0:dt:3*T1,initials_SRC1,[],lambda1,a1,b1,c1,nu1);
eta1 = G1(:,8:14)';

% Caculate the average of y4
% dteps1 = (Teps1/T1)*dt; tspan2 = t0:dteps1:Teps1;
% y4p = G1(:,7);
% min_lc = min(y4p); %minimum of y4 around the limit cycle
% max_lc = max(y4p); %maximum of y4 around the limit cycle
% y4bar_p = (1/Teps1)*(trapz(y4p)*dteps1); %average of y4

% Caculate the derivative of y4 average, with respect to lambda, using the 
% iSRC method
% [Tu,Pu] = ode45(@FF_model,t0:dt:T1,initials1,[],lambda1,a1,b1,c1);
% y4 = Pu(:,7); 
% Jy4 = zeros(length(Tu),7);
% inte_y4 = zeros(length(Tu),1);
% for i = 1:length(Tu)
%     Jy4(i,:) = [0 0 0 0 0 0 1];
%     inte_y4(i) = Jy4(i,:)*eta1(:,i);
% end
% dy4bar_ana = (1/T1)*(dt*trapz(inte_y4));

%% case 2: follow the same procedures as case 1
a2 = 0.65; b2 = 0.01; c2 = 4.88;
lambda2 = 10;

% find the exact intersection point on the Poincare section {y4=1.6,dy4/dt>0}
% and the limit cycle as the initial point for further running
[~,PP2] = ode45(@FF_model,[t0 tF],initials,[],lambda2,a2,b2,c2);
options = odeset('Events',@FF_event);
[~,~,~,initials2,~] = ode45(@FF_model,[0 tF],PP2(end,:),options,lambda2,a2,b2,c2);  %initial2 is the initial point
[~,~,T2,~,~] = ode45(@FF_model,[0 tF],initials2,options,lambda2,a2,b2,c2);
T2 = T2(end);  %period

[Tu2,Pu2] = ode45(@FF_model,t0:dt:3*T2,initials2,[],lambda2,a2,b2,c2);
y42 = Pu2(:,7); 

% find the monodromy matrix
initials_funda2 = [initials2,...
                   1 0 0 0 0 0 0, 0 1 0 0 0 0 0, 0 0 1 0 0 0 0,...
                   0 0 0 1 0 0 0, 0 0 0 0 1 0 0, 0 0 0 0 0 1 0, 0 0 0 0 0 0 1];

[~,P2] = ode45(@FF_fundamental,t0:dt:T2,initials_funda2,[],lambda2,a2,b2,c2);
U2 = P2(end,8:56);

M2 = [U2(1:7); U2(8:14); U2(15:21); U2(22:28); U2(29:35); U2(36:42); U2(43:49)]';  %monodromy matrix
[v2 d2] = eig(M2);  %The monodromy matrix has an eigenvalue equal to 1

% iPRC
tspan_z2 = t0:dt:100*T2;
[~,Q2] = ode45(@FF_model,tspan_z2,initials2,[],lambda2,a2,b2,c2);
xq12 = Q2(:,1); xq22 = Q2(:,2); xq32 = Q2(:,3); 
yq12 = Q2(:,4); yq22 = Q2(:,5); yq32 = Q2(:,6); yq42 = Q2(:,7); 
FFq2 = zeros(7,length(tspan_z2)); DFq2 = zeros(7,7,length(tspan_z2));
for i = 1:length(tspan_z2)
    FFq2(:,i) = F(xq12(i),xq22(i),xq32(i),yq12(i),yq22(i),yq32(i),yq42(i),a2,b2,c2,lambda2);
    DFq2(:,:,i) = A(xq12(i),xq32(i),yq12(i),yq42(i),a2,b2,c2);
end

z2 = zeros(7,length(tspan_z2));
z20 = real(v2(:,3)); %The initial condition of the iPRC: the eigevector of M associated with unit eigenvalue
k0 = FFq2(:,end)'*z20;
z2(:,length(tspan_z2)) = z20/k0;
for i = length(tspan_z2)-1:-1:1  %backward integration
    z2(:,i) = z2(:,i+1)+DFq2(:,:,i+1)'*z2(:,i+1)*dt;
    k = FFq2(:,i)'*z2(:,i);
    z2(:,i) = z2(:,i)/k;
end

zz20 = z2(:,1); %take the last point of last run as the initial condition for z
tspan2 = t0:dt:T2;  %run for one period
[~,QQ2] = ode45(@FF_model,tspan2,initials2,[],lambda2,a2,b2,c2);
xqq12 = QQ2(:,1); xqq22 = QQ2(:,2); xqq32 = QQ2(:,3); 
yqq12 = QQ2(:,4); yqq22 = QQ2(:,5); yqq32 = QQ2(:,6); yqq42 = QQ2(:,7); 
FF2 = zeros(7,length(tspan2)); DF2 = zeros(7,7,length(tspan2));
for i = 1:length(tspan2)
    FF2(:,i) = F(xqq12(i),xqq22(i),xqq32(i),yqq12(i),yqq22(i),yqq32(i),yqq42(i),a2,b2,c2,lambda2);
    DF2(:,:,i) = A(xqq12(i),xqq32(i),yqq12(i),yqq42(i),a2,b2,c2);
end

zz2 = zeros(7,length(tspan2));
k0 = FF2(:,end)'*zz20;
zz2(:,length(tspan2)) = zz20/k0;
for i = length(tspan2)-1:-1:1
    zz2(:,i) = zz2(:,i+1)+DF2(:,:,i+1)'*zz2(:,i+1)*dt;
    k = FF2(:,i)'*zz2(:,i);
    zz2(:,i) = zz2(:,i)/k;
end

int2 = zeros(length(tspan2),1);
for i = 1:length(tspan2)
    int2(i) = zz2(:,i)'*dFp;
end
T1_prc2 = -dt*trapz(int2);

% iSRC
nu2 = T1_prc2/T2;

% find the perturbed point on the Poincare section and the perturbed limit
% cycle
lambda_p2 = lambda2+eps;
[~,Pp2] = ode45(@FF_model,[t0 tF],initials2,[],lambda_p2,a2,b2,c2);
options = odeset('Events',@FF_event);
[~,~,~,initials_p2,~] = ode45(@FF_model,[t0 tF],Pp2(end,:),options,lambda_p2,a2,b2,c2);

eta20 = (initials_p2-initials2)/eps;  
initials_SRC2 = [initials2,eta20];

[~,G2] = ode45(@FF_SRC,t0:dt:3*T2,initials_SRC2,[],lambda2,a2,b2,c2,nu2);
eta2 = G2(:,8:14)';


%% Plot
figure(1)
%for case 1
subplot(2,2,1) 
plot(Tu1,y41,'-k','LineWidth',2.5)
ylabel('y_4')
axis([0 3*T1 1.514 1.688])
title('(a) case 1')
subplot(2,2,3)
plot(Tu1,eta1(7,:),'-b','LineWidth',2.5);
hold on
plot(Tu1,zeros(1,length(Tu1)),'--r','LineWidth',2.1);
hold off
axis([0 3*T1 -.097 .097])
xlabel('t'); ylabel('\eta_{y_4}')
%for case 2
subplot(2,2,2) 
plot(Tu2,y42,'-k','LineWidth',2.5)
ylabel('y_4')
axis([0 3*T2 1.514 1.688])
title('(a) case 2')
subplot(2,2,4)
plot(Tu2,eta2(7,:),'-b','LineWidth',2.5);
hold on
plot(Tu2,zeros(1,length(Tu2)),'--r','LineWidth',2.1);
hold off
axis([0 3*T2 -.097 .097])
xlabel('t'); ylabel('\eta_{y_4}')