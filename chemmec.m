clear
clc
close all

%Difine the Components 
COH    =C(1);
CO     =C(2);
CH2    =C(3);
CO2    =C(4);
CH2O   =C(5);
CH2O2  =C(6);
CCO    =C(7);
CCH4   =C(8);
CCH3   =C(9);
CCH2O  =C(10);
CHO    =C(11);
CH     =C(12);
CCH2   =C(13);

% Rate Constant K 

% Rate Law r = K()*CR1*CR2*....
r1  = k1*COH*CO;
r2  = k2*CO*CH2;
r3  = k3*COH*CH2;
r4  = k4*COH*COH;
r5  = k5*CH*CH;
r6  = k6*CH*COH;
r7  = k7*CO*CO;
r8  = k8*CH*CO2;
r9  = k9*CH*CHO2;
r10 = k10*CH*CHO2;
r11 = k11*CO*CHO2;
r12 = k12*COH*CHO2;
r13 = k13*CHO2*CHO2;
r14 = k14*COH*COH;
r15 = k15*CH*CH2O2;
r16 = k16*CH*CH2O2;
r17 = k17*CO*CH2O2;
r18 = k18*COH*CH2O2;
r19 = k19*CCO*CH;
r20 = k20*CCO*CO;
r21 = k21*CCO*COH;
r22 = k22*CCO*CHO2;
r23 = k23*CCO*CO2;
r24 = k24*CCH4*CH;
r25 = k25*CCH4*CO;
r26 = k26*CCH4*COH;
r27 = k27*CCH4;
r28 = k28*CCH3*CO;
r29 = k29*CCH3*CO2;
r30 = k30*CCH2O*CH;
r31 = k31*CCH2O*CO;
r32 = k32*CCH2O*COH;
r33 = k33*CCH2O*CCH2;
r34 = k34*CCH2O;
r35 = k35*CCH2O;
r36 = k36*CCHO*CH;
r37 = k37*CCHO*CO;
r38 = k38*CCHO*CO;
r39 = k39*CCHO*COH;
r40 = k40*CCHO*CO2;

% Mass Balance dR1 = r1 + r2 +... 
dOH   = -r1 + r2 - r3 - 2*r4 - r6 + 2*r9 + r11 - r12 - 2*r14 + r16 + r17 - r18 - r21 + r25 - r26 + r31 - r32 + r37 - r39;
dO    = -r1 - r2 + r4 - 2*r7 - r11 - r17 - r20 + r22 + r23 - r25 - r28 + r29 - r31 - r37 - r38;
dH2   = -r2 - r3 + r4 + r5 + r10 + r15 + r24 + r30 + r35 + r36;
dH    =  r1 + r2 + r3 -2*r5 - r6 - r8 - r9 - r10 - r15 - r16 - r19 + r21 - r24 + r27 + r28 - r30 + r34 - r36 + r38;
dO2   =  r1 + r7 - r8 + r10 + r11 + r12 + r13 - r23 - r29 - r40;
dHO2  =  r8 - r9 - r10 - r11 - r12 - 2*r13 + r15 + r17 + r18 - r22 + r40;
dH2O2 =  r13 + r14 - r15 - r16 - r17 - r18; 
dCO   = -r19 - r20 - r21 - r22 - r23 + r35 + r36 + r37 + r39 + r40;
dCH4  = -r24 - r25 - r26 - r27 + r33; 
dCH3  =  r24 + r25 + r26 + r27 - r28 - r29;
dCH2O =  r28 + r29 - r30 - r31 - r32 - r33 - r34 - r35;
dCHO  =  r19 + r30 + r31 + r32 + r33 + r34 - r36 - r37 - r38 - r39 - r40;
dCH2  = -r33;
dH2O  =  r3 + r6 + r12 + r16 + r18 + r26 + r32 +r39;
% Assign The Output Variables 
dC(1,:)  = dOH;
dC(2,:)  = dO;
dC(3,:)  = dH2;
dC(4,:)  = dH;
dC(5,:)  = dO2;
dC(6,:)  = dHO2;
dC(7,:)  = dH2O2;
dC(8,:)  = dCO;
dC(9,:)  = dCH4;
dC(10,:) = dCH3;
dC(11,:) = dCH2O;
dC(12,:) = dCHO;
dC(13,:) = dCH2;
dC(14,:) = dH2O;