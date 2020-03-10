%PD week 4
%NSAID model
%aspirin clotting
%attempt 2

close all

dose = 500/0.15; %initial dose (mg/L)

Cox0 = 0.0055; %initial Cox
PT0 = 0.000013; %initial PGE2/TxA2
PLat0 = 450000000000; %initial number of platelets/L
tspan = [0 500];
C0 = [dose 0 0 0 ...
      Cox0 PT0 0 0 0 0 PLat0 Cox0 PT0 0 0 0 0 PLat0 Cox0 PT0 0 0 0 0 PLat0 Cox0 PT0 0 0 0 0 PLat0...
      Cox0 PT0 0 0 0 0 PLat0 Cox0 PT0 0 0 0 0 PLat0 Cox0 PT0 0 0 0 0 PLat0 Cox0 PT0 0 0 0 0 PLat0];

%run model
[t,y] = ode15s(@odefcn,tspan,C0);

%clot = 8*PT0 - (y(:,6)+y(:,12)+y(:,19)+y(:,26)+y(:,33)+y(:,40)+y(:,47)+y(:,54));
clot = (y(:,6)+y(:,13)+y(:,20)+y(:,27)+y(:,34)+y(:,41)+y(:,48)+y(:,55));
cox = (y(:,5)+y(:,12)+y(:,19)+y(:,26)+y(:,33)+y(:,40)+y(:,47)+y(:,54));

figure
plot(t,cox)
xlabel('time')
ylabel('Cox')

figure
plot(t,clot)

figure
plot(t,y(:,4))



function dCdt = odefcn(t,C)
Ke = 0.693/15; %1/min
Ka = 0.029; %rate of absorption
Kph = 1.13; %L/min plasma to hand
Khp = 1.13;
Kps = 0.7; %L/min plasma to stomach
Ksp = 0.971; %L/min
Kpt = 0.0024; %plasma to platelet
Ktp = 0.0024;
K1 = 0.01; %production of PGE2/TxA2 (binding to AA)
K2 = 0.01; %production of A-Cox (unbinding ASA)
Keq1 = 0.001; %eq const for ASA + Cox -> ASA-Cox
Keq2 = 0.001; %eq const for Cox + SA -> SA-Cox
%Kday = 1/(60*24)*0.0055; %rate of platelets moving to next stage
Kborn = 1E11/(24*60)*0.0055; %1/min rate of platelet production*Cox/plat
Keq2on = 1/0.001; %Cox + SA -> SA-Cox
Keq2off = 0.001; %SACox -> Cox + SA
Kina = 1/(4*24*60); %rate of platelet deactivation halflife 8 days
Kday = 1/(4*24*60); %rate of platelets moving to next stage
Cox0 = 0.0055; %initial Cox

dCdt = [%Compartments------------------------------------------------------
        -Ka*C(1);...%1. GI ASA
        Ka*C(1)-(Kph+Kps+Kpt*8+Ke)*C(2)+Khp*C(3)+Ksp*C(4)+Ktp*(C(9)+C(15)+C(21)+C(27)+C(33)+C(39)+C(45)+C(51));...%2. plasma ASA
        Kph*C(2) - Khp*C(3);...%3. hand ASA
        Kps*C(2) - Ksp*C(4);...%4. stomach ASA
        %platelet 1--------------------------------------------------------
        Kborn*Cox0 - K1*C(5) - K2*Keq1*C(9)*C(5) - Keq2on*C(5)*C(10)+Keq2off*C(8) - Kday*C(5);...%5. plat1 UB Cox
        K1*C(5);...%6. plat1 PGE2/TxA2
        K2*Keq1*C(9)*C(5);...%7. plat1 A-Cox
        Keq2*C(5)*C(9);...%8. plat1 SA-Cox
        Kpt*C(2) - Ktp*C(9);...%9. plat1 ASA
        K2*Keq1*C(5)*C(9);...%10. plat1 SA
        Kborn - Kday*C(11) - Kina*C(5)/Cox0;...%11. Active Platelets 1
        %platelet 2--------------------------------------------------------
        Kday*C(5) - K1*C(12) - K2*Keq1*C(16)*C(12) - Keq2on*C(12)*C(17)+Keq2off*C(15) - Kday*C(12);...%12.UB Cox
        K1*C(12);...%13. PGE2/TxA2
        K2*Keq1*C(16)*C(12);...%14. A-Cox
        Keq2*C(12)*C(16);...%15. SA-Cox
        Kpt*C(2) - Ktp*C(16);...%16. ASA
        K2*Keq1*C(12)*C(16);...%17. SA
        Kday*C(11) - Kday*C(18) - Kina*C(12)/Cox0;...%18. Active Platelets 2
        %platelet 3--------------------------------------------------------
        Kday*C(12) - K1*C(19) - K2*Keq1*C(23)*C(19) - Keq2on*C(19)*C(24)+Keq2off*C(22) - Kday*C(19);...%19.UB Cox
        K1*C(19);...%20. PGE2/TxA2
        K2*Keq1*C(23)*C(19);...%21. A-Cox
        Keq2*C(19)*C(23);...%22. SA-Cox
        Kpt*C(2) - Ktp*C(23);...%23. ASA
        K2*Keq1*C(19)*C(23);...%24. SA
        Kday*C(18) - Kday*C(25) - Kina*C(19)/Cox0;...%25. Active Platelets 3
        %platelet 4--------------------------------------------------------
        Kday*C(19) - K1*C(26) - K2*Keq1*C(30)*C(26) - Keq2on*C(26)*C(31)+Keq2off*C(29) - Kday*C(26);...%26.UB Cox
        K1*C(26);...%27. PGE2/TxA2
        K2*Keq1*C(30)*C(26);...%28. A-Cox
        Keq2*C(26)*C(30);...%29. SA-Cox
        Kpt*C(2) - Ktp*C(30);...%30. ASA
        K2*Keq1*C(26)*C(30);...%31. SA
        Kday*C(25) - Kday*C(32) - Kina*C(26)/Cox0;...%32. Active Platelets 4
        %platelet 5--------------------------------------------------------
        Kday*C(26) - K1*C(33) - K2*Keq1*C(37)*C(33) - Keq2on*C(33)*C(38)+Keq2off*C(36) - Kday*C(33);...%33.UB Cox
        K1*C(33);...%34. PGE2/TxA2
        K2*Keq1*C(37)*C(33);...%35. A-Cox
        Keq2*C(33)*C(37);...%36. SA-Cox
        Kpt*C(2) - Ktp*C(37);...%37. ASA
        K2*Keq1*C(33)*C(37);...%38. SA
        Kday*C(32) - Kday*C(39) - Kina*C(33)/Cox0;...%39. Active Platelets 5
        %platelet 6--------------------------------------------------------
        Kday*C(33) - K1*C(40) - K2*Keq1*C(44)*C(40) - Keq2on*C(40)*C(45)+Keq2off*C(43) - Kday*C(40);...%40.UB Cox
        K1*C(40);...%41. PGE2/TxA2
        K2*Keq1*C(44)*C(40);...%42. A-Cox
        Keq2*C(40)*C(44);...%43. SA-Cox
        Kpt*C(2) - Ktp*C(44);...%44. ASA
        K2*Keq1*C(40)*C(44);...%45. SA
        Kday*C(39) - Kday*C(46) - Kina*C(40)/Cox0;...%46. Active Platelets 6
        %platelet 7--------------------------------------------------------
        Kday*C(40) - K1*C(47) - K2*Keq1*C(51)*C(47) - Keq2on*C(47)*C(52)+Keq2off*C(50) - Kday*C(47);...%47.UB Cox
        K1*C(47);...%48. PGE2/TxA2
        K2*Keq1*C(51)*C(47);...%49. A-Cox
        Keq2*C(47)*C(51);...%50. SA-Cox
        Kpt*C(2) - Ktp*C(51);...%51. ASA
        K2*Keq1*C(47)*C(51);...%52. SA
        Kday*C(46) - Kday*C(53) - Kina*C(47)/Cox0;...%53. Active Platelets 7
        %platelet 8--------------------------------------------------------
        Kday*C(47) - K1*C(54) - K2*Keq1*C(58)*C(54) - Keq2on*C(54)*C(59)+Keq2off*C(57) - Kday*C(54);...%54.UB Cox
        K1*C(54);...%55. PGE2/TxA2
        K2*Keq1*C(58)*C(54);...%56. A-Cox
        Keq2*C(54)*C(58);...57. SA-Cox
        Kpt*C(2) - Ktp*C(58);...%58. ASA
        K2*Keq1*C(54)*C(58);...%59. SA
        Kday*C(53) - Kday*C(60) - Kina*C(54)/Cox0;...%60. Active Platelets 8
        ];
end 
