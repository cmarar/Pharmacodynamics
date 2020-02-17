%direct link model week 3
%CGP 51901

close all

%initialize----------------------------------------------------------------
dose = 60/6.02E23; %initial dose (mol/L)
IgE = 100/6.02E23; %baseline IgE conc in mol/L
FcR = 0.00634; %mol/L FcR taken from stelara model

C0 = [dose, IgE, 0, 0, FcR, 0];
tspan = [0 14];

%run model
[t,y] = ode15s(@odefcn,tspan,C0);

%multiple doses------------------------------------------------------------
y_G = y(:,1);
y_E = y(:,2);
t2 = t(end);
t_cr = t;

C0 = y(end,:);
C0(1) = C0(1)+dose;

for i = 1:6
    [t1,y1] = ode15s(@odefcn,tspan,C0);

    y_G = [y_G;y1(:,1)];
    y_E = [y_E;y1(:,2)];
    t2 = t1+t2(end);
    
    C0 = y1(end,:);
    C0(1) = C0(1)+dose;
    
    t_cr = [t_cr;t2];
end

%calculate percentage reduction
y_red = y(:,2)/IgE;


%plot
figure
plot(y(:,1), y_red)
xlabel('CGP 51901 (mol/L)')
ylabel('% IgE reduction')

%model
function dCdt = odefcn(t,C)
%rate constants (1/day)-----------------
k_d = 0.046;%degradation
k_21 = 0.161; %flow rate to peripheral compartment
k_prod = 2.9E-21; %rate of IgE production mol/day

%binding constants
k_ge = 8E9; %IgG to IgE
k_fe = 8E9; %IgE to FcR

%Volume of distribution (L)
V_1 = 3; 
V_2 = 0.45;

%Equilibrium expressions-------------
C(3) = k_ge*C(1)*C(2);
C(4) = k_fe*C(2)*C(5);


%differential equations--------------
dCdt = [- k_21*C(1) + k_21*C(6) ;... %IgG free comp 1
        k_prod;...%IgE free
        0;...
        0;...
        0;...
        k_21*C(1) - k_21*C(6) - k_d*C(6)];
%1. IgG free comp 1
%2. IgE free
%3. IgG-IgE
%4. FcR-IgE
%5. FcR free
%6. IgG free comp 2

end 
