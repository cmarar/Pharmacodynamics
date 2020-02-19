%PD week 3
%receptor transducer model

close all


dose = 0.396/10; %initial dose (ug/ml)

C0 = [dose 100 0];
tspan = [0 24];

%run model
[t,y] = ode15s(@odefcn,tspan,C0,0);


Emax = 0.8;
EC50 = 3.7E-6; %ug/ml
E = Emax*y./(EC50 + y(:,1));

% figure
% plot(y,R)
% xlabel('Drug Concentration (ug/ml)')
% ylabel('Response (Osmolarity)')



function dCdt = odefcn(t,C,x)
Ke = .693/4.36; %1/hr
Kd = 0.04;

dCdt = [-Ke*C(1)+x; 0; -x];
x = ode15s(@odefcn2,dCdt);

end 

function dxdt = odefcn2(dCdt)
Kd = 0.04;
dxdt = -dCdt(1)*C(2) - Kd*dCdt(3);

end 




