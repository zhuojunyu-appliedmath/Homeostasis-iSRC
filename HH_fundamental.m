function dPdt = HH_fundamental(t,P,C,I,gK,EK,gNa,ENa,gL,EL)
% This function find the fundamental matrix of the model

V=P(1); n=P(2); m=P(3); h=P(4); 
u1 = P(5:8); u2 = P(9:12); u3 = P(13:16); u4 = P(17:20);

alpha_n = 0.01*(10-V)/(exp((10-V)/10)-1);
alpha_m = 0.1*(25-V)/(exp((25-V)/10)-1);
alpha_h = 0.07*exp(-V/20);
beta_n = 0.125*exp(-V/80);
beta_m = 4*exp(-V/18);
beta_h = 1/(exp((30-V)/10)+1);

dalpha_n = 0.01*(1-0.1*V*exp((10-V)/10))/(exp((10-V)/10)-1)^2;
dalpha_m = 0.1*((1.5-0.1*V)*exp((25-V)/10)+1)/(exp((25-V)/10)-1)^2;
dalpha_h = -(0.07/20)*exp(-V/20);
dbeta_n = -(0.125/80)*exp(-V/80);
dbeta_m = -(4/18)*exp(-V/18);
dbeta_h = 0.1*exp((30-V)/10)/(exp((30-V)/10)+1)^2;

A =  @(V,n,m,h) [(1/C)*(-gK*n^4-gNa*m^3*h-gL) (1/C)*(-4*gK*n^3*(V-EK)) (1/C)*(-3*gNa*m^2*h*(V-ENa)) (1/C)*(-gNa*m^3*(V-ENa));
                 (1-n)*dalpha_n-n*dbeta_n -alpha_n-beta_n 0 0;
                 (1-m)*dalpha_m-m*dbeta_m 0 -alpha_m-beta_m 0;
                 (1-h)*dalpha_h-h*dbeta_h 0 0 -alpha_h-beta_h];

dVdt = (1/C)*(I-gK*n^4*(V-EK)-gNa*m^3*h*(V-ENa)-gL*(V-EL));
dndt = alpha_n*(1-n)-beta_n*n;
dmdt = alpha_m*(1-m)-beta_m*m;
dhdt = alpha_h*(1-h)-beta_h*h;

du1dt = A(V,n,m,h)*u1;
du2dt = A(V,n,m,h)*u2;
du3dt = A(V,n,m,h)*u3;
du4dt = A(V,n,m,h)*u4;

dPdt = [dVdt;dndt;dmdt;dhdt; du1dt; du2dt; du3dt; du4dt];

end