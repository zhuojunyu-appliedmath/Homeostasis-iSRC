function dPdt = FF_fundamental(t,P,lambda,a,b,c)
% This function computes the fundamental matrix of the model

x1=P(1); x2=P(2); x3=P(3);
y1=P(4); y2=P(5); y3=P(6); y4=P(7);
u1 = P(8:14); u2 = P(15:21); u3 = P(22:28); 
u4 = P(29:35); u5 = P(36:42); u6 = P(43:49); u7 = P(50:56);

zeta = @(y4) 10/(1+y4^10)+b;
eta = @(x1) 1/(1+exp((c-x1)/a));

% Jacobian 
A = @(x1,x3,y1,y4) [-2 0 0 0 0 0 0;
                    1 -2 0 0 0 0 0;
                    -(x3*exp((c-x1)/a))/(a*(1+exp((c-x1)/a))^2) 1 -(1+eta(x1)) 0 0 0 0;
                    0 0 1 -zeta(y4) 0 0 (100*y1*y4^9)/(1+y4^10)^2;
                    0 0 0 zeta(y4) -1 0 -(100*y1*y4^9)/(1+y4^10)^2;
                    0 0 0 0 1 -1 0;
                    0 0 0 0 0 1 -1];

dx1dt = lambda-2*x1;
dx2dt = x1-2*x2;
dx3dt = x2-(1+eta(x1))*x3;
dy1dt = x3-zeta(y4)*y1;
dy2dt = zeta(y4)*y1-y2;
dy3dt = y2-y3;
dy4dt = y3-y4;

du1dt = A(x1,x3,y1,y4)*u1;
du2dt = A(x1,x3,y1,y4)*u2;
du3dt = A(x1,x3,y1,y4)*u3;
du4dt = A(x1,x3,y1,y4)*u4;
du5dt = A(x1,x3,y1,y4)*u5;
du6dt = A(x1,x3,y1,y4)*u6;
du7dt = A(x1,x3,y1,y4)*u7;

dPdt = [dx1dt;dx2dt;dx3dt;dy1dt;dy2dt;dy3dt;dy4dt;...
        du1dt;du2dt;du3dt;du4dt;du5dt;du6dt;du7dt];

end