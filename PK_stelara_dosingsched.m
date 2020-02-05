%Stelara multicompartment
%SC Dosing Schedule
%2/4/20

close all

%initialize parameters--------------------------------------------
tspan = [0 120960]; %tspan is 12 weeks in minutes
dose = 6.179279E-4; %initial SC dose in M
FcR_init = 0.00634; %M FcR in cells
FcRn_init = 1660E-8; %M FcRn in cells
C0 = [dose,FcR_init,0,0,0,0,FcRn_init,0,0,0,0,0,0,0,... %SC
    0,FcR_init,0,0,0,0,FcRn_init,0,0,0,0,0,0,0,... %lymph
    0,FcR_init,0,0,0,0,FcRn_init,0,0,0,0,0,0,0,... %plasma
    0,0,0]; % surface FcRn

%run ODE---------------------------------------------------------
[t,y] = ode15s(@odefcn,tspan,C0);

%multiple doses--------------------------------------------------
A = y(end,:);%putt out the last values for everything
y_tot = y(:,:);
y_P = y(:,29);

t2 = t(end);
t_P = t;


for i = 1:3
    C0 = A;%C0 starts where last left off
    
    [t1,y1] = ode45(@odefcn,tspan,C0);

    y_P = [y_P; y1(:,29)];
    A = y1(end,:);
    t2 = t1+t2(end);
    t_P = [t_P;t2];
end

%plot------------------------------------------------------------

plot(t_P,y_P)

xlabel('time (weeks)')
ylabel('concentration (M)')
xticks([0:120960:120960*4])
xticklabels({'0','12','24','48'})
%ylim([0 5.9*10^-8])
title('Stelara Psoriasis')

function dCdt = odefcn(t,C)
%flow rates (L/min)-----------------
k_d = 0.000229;%degradation
k_pino_sc = 1.3787E-9;%pinocytosis SC
k_pino_L = 1.3787E-11;%pinocytosis lymph
k_pino_P = 6.6445E-13;%pinocytosis plasma
k_in = .0462; %rate of A-FcR entering cell-----ESTIMATE
k_r = 0.000924; %endosomal recycling rate
k_L = 0.0000417/5; %flow rate sc to lymph
k_P = 0.000056/5; %flow rate sc to plasma
k_LP = 0.000278; %flow rate lymph to plasma

k_fcr_out = 10; %Ka for fcr binding outside cell-----random number, was too high
k_fcr_in = 0.0000001; %Ka for fcr binding in cell-----ESTIMATE
k_n_out = 0.0000001; %Ka for fcrn binding outside cell-----ESTIMATE
k_n_in = 1315780; %Ka for fcrn binding in cell

R = 100; %Partition Coefficient for stelara

%Equilibrium expressions-------------
%SC
C(3) = C(1)*C(2)*k_fcr_out;
C(13) = C(14)*C(8)*k_fcr_in;
C(9) = C(1)*C(43)*k_n_out;
C(4) = C(12)*C(7)*k_n_in;
C(8) = C(14)*C(1)*k_fcr_in;
C(11) = C(1)*C(9)*k_n_out;
C(5) = C(12)*C(4)*k_n_in;

%lymph
C(17) = C(15)*C(16)*k_fcr_out;
C(23) = C(28)*C(26)*k_fcr_in;
C(27) = C(23)*C(26)*k_fcr_in;
C(18) = C(15)*C(44)*k_n_out;
C(22) = C(28)*C(21)*k_n_in;
C(20) = C(15)*C(18)*k_n_out;
C(25) = C(28)*C(22)*k_n_in;

%plasma
C(31) = C(29)*C(30)*k_fcr_out;
C(37) = C(42)*C(40)*k_fcr_in;
C(41) = C(37)*C(40)*k_fcr_in;
C(32) = C(29)*C(45)*k_n_out;
C(36) = C(42)*C(35)*k_n_in;
C(34) = C(29)*C(32)*k_n_out;
C(39) = C(42)*C(36)*k_n_in;

dCdt = [
    %Subcutaneous---------------------
    k_r*(C(4)+C(5)+C(6)) - k_L/R*C(1) - k_P/R*C(1) - k_pino_sc - k_in*C(3);...%1
    0;...%2
    -k_in*C(3);...%3
    -k_r*C(4);...%4
    -k_r*C(5);...%5
    -k_r*C(6);...%6
    0;...%7
    -k_d*C(8)+k_in*C(3);...%8
    k_r*C(4);...%9
    k_r*C(6);...%10
    k_r*C(5);...%11
    k_pino_sc*C(1);...%12
    -k_d*C(13);...%13
    0;...%14
    %Lymph-----------------------------
    k_L/R*C(1) - k_in*C(17)+ k_r*(C(22)+C(24)+C(25)) - k_LP*C(15)-k_pino_L;...%15
    0;...%16
    -k_in*C(17);...%17
    k_r*C(22);...%18
    k_r*C(24);...%19
    k_r*C(25);...%20
    0;...%21
    -k_r*C(22);...%22
    k_in*C(17);...%23
    -k_r*C(24);...%24
    -k_r*C(25);...%25
    0;...%26
    -k_d*C(27);...%27
    k_pino_L - k_d*C(28);...%28
    %Plasma----------------------------
    k_P/R*C(1) + k_LP*C(15) - k_in*C(31) + k_r*(C(36)+C(38)+C(39))- k_pino_P;...%29
    0;...%30
    -k_in*C(31);...%31
    k_r*C(36);...%32
    k_r*C(38);...%33
    k_r*C(39);...%34
    0;...%35
    -k_r*C(36);...%36
    k_in*C(31);...%37
    -k_r*C(38);...%38
    -k_r*C(39);...%39
    0;...%40
    -k_d*C(41);...%41
    k_pino_P - k_d*C(42);...%42
    %surface FcRn----------------------
    0;...%43
    0;...%44
    0]; %45

% 1. A_SC
% 2. FcR surface 
% 3. A-FcR surface
% 4. A-FcRn cell
% 5. FcRn-A-FcRn cell
% 6. FcR-A-FcRn cell
% 7. FcRn cell
% 8. A-FcR cell
% 9. A-FcRn surface
% 10.FcR-A-FcRn surface
% 11.FcRn-A-FcRn surface
% 12.A cell
% 13.FcR-A-FcR cell
% 14.FcR cell

% 15.A lymph
% 16.FcR surface
% 17.A-FcR surface
% 18.A-FcRn surface
% 19.FcR-A-FcRn surface
% 20.FcRn-A-FcRn surface
% 21.FcRn cell
% 22.A-FcRn cell
% 23.A-FcR cell
% 24.FcR-A-FcRn cell
% 25.FcRn-A-FcRn cell
% 26.FcR cell
% 27.FcR-A-FcR cell
% 28.A cell


% 29.A plasma
% 30.FcR surface
% 31.A-FcR surface
% 32.A-FcRn surface
% 33.FcR-A-FcRn surface
% 34.FcRn-A-FcRn surface
% 35.FcRn cell
% 36.A-FcRn cell
% 37.A-FcR cell
% 38.FcR-A-FcRn cell
% 39.FcRn-A-FcRn cell
% 40.FcR cell
% 41.FcR-A-FcR cell
% 42.A cell

%surface FcRn
% 43. FcRn SC surface
% 44. FcRn lymph surface
% 45. FcRn plasma surface
end 
