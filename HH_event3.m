function [value,isterminal,direction] = HH_event3(t,P,C,I,gK,EK,gNa,ENa,gL,EL)
% This function finds the point transverse to Poincare section {V=20,dV/dt>0}

value=P(1)-20;
isterminal=1;
direction=1;

end