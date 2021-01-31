function dPdt = HH_SRC(t,P,C,I,gK,EK,gNa,ENa,gL,EL,nu1)
% This function calculates the iSRC which satisfies an inhomogeneous
% variational equation

V=P(1); n=P(2); m=P(3); h=P(4); y1=P(5); y2=P(6); y3=P(7); y4=P(8);
y=[y1;y2;y3;y4];

alpha_n = @(V) 0.01*(10-V)/(exp((10-V)/10)-1);
alpha_m = @(V) 0.1*(25-V)/(exp((25-V)/10)-1);
alpha_h = @(V) 0.07*exp(-V/20);
beta_n = @(V) 0.125*exp(-V/80);
beta_m = @(V) 4*exp(-V/18);
beta_h = @(V) 1/(exp((30-V)/10)+1);

Fu = @(V,n,m,h) [(1/C)*(I-gK*n^4*(V-EK)-gNa*m^3*h*(V-ENa)-gL*(V-EL));
                 alpha_n(V)*(1-n)-beta_n(V)*n;
                 alpha_m(V)*(1-m)-beta_m(V)*m;
                 alpha_h(V)*(1-h)-beta_h(V)*h];

dFp = [1/C; 0; 0; 0];

dalpha_n = @(V) 0.01*(1-0.1*V*exp((10-V)/10))/(exp((10-V)/10)-1)^2;
dalpha_m = @(V) 0.1*((1.5-0.1*V)*exp((25-V)/10)+1)/(exp((25-V)/10)-1)^2;
dalpha_h = @(V) -(0.07/20)*exp(-V/20);
dbeta_n = @(V) -(0.125/80)*exp(-V/80);
dbeta_m = @(V) -(4/18)*exp(-V/18);
dbeta_h = @(V) 0.1*exp((30-V)/10)/(exp((30-V)/10)+1)^2;

A = @(V,n,m,h) [(1/C)*(-gK*n^4-gNa*m^3*h-gL) (1/C)*(-4*gK*n^3*(V-EK)) (1/C)*(-3*gNa*m^2*h*(V-ENa)) (1/C)*(-gNa*m^3*(V-ENa));
                (1-n)*dalpha_n(V)-n*dbeta_n(V) -alpha_n(V)-beta_n(V) 0 0;
                (1-m)*dalpha_m(V)-m*dbeta_m(V) 0 -alpha_m(V)-beta_m(V) 0;
                (1-h)*dalpha_h(V)-h*dbeta_h(V) 0 0 -alpha_h(V)-beta_h(V)];
            
dVdt = (1/C)*(I-gK*n^4*(V-EK)-gNa*m^3*h*(V-ENa)-gL*(V-EL));
dndt = alpha_n(V)*(1-n)-beta_n(V)*n;
dmdt = alpha_m(V)*(1-m)-beta_m(V)*m;
dhdt = alpha_h(V)*(1-h)-beta_h(V)*h;
dydt = A(V,n,m,h)*y+nu1*Fu(V,n,m,h)+dFp;

dPdt = [dVdt;dndt;dmdt;dhdt;dydt];

end