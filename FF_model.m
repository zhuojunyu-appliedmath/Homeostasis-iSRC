function dPdt = FF_model(t,P,lambda,a,b,c)
% 7-dim ODEs for a feed-forward network

x1=P(1); x2=P(2); x3=P(3); 
y1=P(4); y2=P(5); y3=P(6); y4=P(7);

eta = @(x1) 1/(1+exp((c-x1)/a));
zeta = @(y4) 10/(1+y4^10)+b;

dx1dt = lambda-2*x1;
dx2dt = x1-2*x2;
dx3dt = x2-(1+eta(x1))*x3;
dy1dt = x3-zeta(y4)*y1;
dy2dt = zeta(y4)*y1-y2;
dy3dt = y2-y3;
dy4dt = y3-y4;

dPdt = [dx1dt;dx2dt;dx3dt;dy1dt;dy2dt;dy3dt;dy4dt];

end