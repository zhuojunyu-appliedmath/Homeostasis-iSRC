function [value,isterminal,direction] = FF_event(t,P,lambda,a,b,c)
% This function finds the point on Poincare section {y4=1.6.dy4>dt>0}

value = P(7)-1.6;
isterminal = 1;
direction = 1;
end