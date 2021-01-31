function [value,isterminal,direction] = HH_event2(t,P,C,I,gK,EK,gNa,ENa,gL,EL)
% This function finds the point transverse to Poincare section 
% {n=0.582483641907492, dn/dt<0}

value=P(2)-0.582483641907492;
isterminal=1;
direction=-1;
end