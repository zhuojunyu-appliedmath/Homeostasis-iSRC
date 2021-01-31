This repository holds code used to generate the figures in our paper "A Homeostasis Criterion for Limit Cycle Systems Based on Infinitesimal Shape Response Curves".
Figures 3-6 are for the Hodgkin-Huxley model (HH). Figures 8-9 are for the feedforward network (FF). The rest figures are not generated via simulation.

Fig. 3 is a reference figure for the Hodgkin-Huxley model showing the evolution of variables V, n, m, h for two different values of the current, I = 50 and I = 60.
To generate figure, run Fig_3.m

Fig. 4 shows the iSRC components for the four variables of the Hodgkin-Huxley model, for two different Poincare sections / initial conditions (S1 = {V=0, dV/dt>0}, S2 = {n=0.582483641907492, dn/dt<0}. 
To generate figure, run Fig_4.m

Fig. 5 compares the difference of the distinct iSRCs obtained from Fig_4.m with the unperturbed flow.
To generate figure, run Fig_5.m

Fig.6 shows the empirical curves for I over [45,55] with analytically derived tangent curves at I = 50 using our iSRC method, for voltage, sodium current and potassium current.
To generate figure, run Fig_6.m. (One needs to first run Fig_4.m to find iSRC eta(t) or xi(t).)

Fig. 8 is a reference figure for the feedforward network showing the time evolution of variable y4 at λ = 10.5 for case 1 and at λ = 10 for case 2, and their respective iSRC in y4 component, defined by section {y4=1.6, dy4/dt>0}.
To generate figure, run Fig_8.m

Fig. 9 shows the behavior of y4 as λ varies in the chain network. 
The top panels show the value of y4 at equilibrium and the amplitudes of limit cycles. Middle panels show the average of y4. Bottom panels show the derivative of y4 with respect to λ, calculated either analytically or directly.
The data of this figure was precomputed via MatCont and MATLAB, with FF_model.m as the simulation code.
Please refer to line 112-117 of Fig_8.m for the way of calculating the results in middle panels and refer to line 119-129 of Fig_8.m for the way of calculating the results in bottom panels.
To produce the plot of case 1 from precomputed data (FF_case1.mat), run Fig_9_case1.m.
To produce the plot of case 2 from precomputed data (FF_case2.mat), run Fig_9_case2.m
