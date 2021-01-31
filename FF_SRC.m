function dPdt = FF_SRC(t,P,lambda,a,b,c,nu1)
% This function finds the iSRC of the model

x1=P(1); x2=P(2); x3=P(3); y1=P(4); y2=P(5); y3=P(6); y4=P(7); 
SRC1=P(8); SRC2=P(9); SRC3=P(10); SRC4=P(11); SRC5=P(12); SRC6=P(13); SRC7=P(14);
SRC=[SRC1;SRC2;SRC3;SRC4;SRC5;SRC6;SRC7];
            
eta = @(x1) 1/(1+exp((c-x1)/a));
zeta = @(y4) 10/(1+y4^10)+b;

Fu = @(x1,x2,x3,y1,y2,y3,y4) [lambda-2*x1; x1-2*x2; x2-(1+eta(x1))*x3; 
                              x3-zeta(y4)*y1; zeta(y4)*y1-y2; y2-y3; y3-y4];
             
dFp = [1; 0; 0; 0; 0; 0; 0];

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
dSRCdt = A(x1,x3,y1,y4)*SRC+nu1*Fu(x1,x2,x3,y1,y2,y3,y4)+dFp;

dPdt = [dx1dt;dx2dt;dx3dt;dy1dt;dy2dt;dy3dt;dy4dt;dSRCdt];

end