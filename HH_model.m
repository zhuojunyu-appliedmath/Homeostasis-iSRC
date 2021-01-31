function dPdt = HH_model(t,P,C,I,gK,EK,gNa,ENa,gL,EL)
% 4-dim ODEs for Hodgkin-Huxley model

V=P(1); n=P(2); m=P(3); h=P(4); 

alpha_n = 0.01*(10-V)/(exp((10-V)/10)-1);
alpha_m = 0.1*(25-V)/(exp((25-V)/10)-1);
alpha_h = 0.07*exp(-V/20);
beta_n = 0.125*exp(-V/80);
beta_m = 4*exp(-V/18);
beta_h = 1/(exp((30-V)/10)+1);

dVdt = (1/C)*(I-gK*n^4*(V-EK)-gNa*m^3*h*(V-ENa)-gL*(V-EL));
dndt = alpha_n*(1-n)-beta_n*n;
dmdt = alpha_m*(1-m)-beta_m*m;
dhdt = alpha_h*(1-h)-beta_h*h;

dPdt = [dVdt;dndt;dmdt;dhdt];

end