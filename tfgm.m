C=0.00102;
Em=400*sqrt(2)/sqrt(3);
R=0.5;
L=0.0054;
k=1.36;
t=0.0045;
wn=sqrt(k*326.6/t)
e=sqrt(t*k*326.6)/2
Q_=-5000;
P_=-7000;
Tc=0.001;
Kpc=L/Tc;
Kic=R/Tc;
E_=800;
edc=(sqrt(2))/2;
wdc=418.825;
Kpdc=C*edc*wdc
Kidc=(C*wdc^2)/2
we=314.16;

%PV Tesla SR25T3
Temp=298; %Temperature K (free parameter)
G0=700; %Irradiance w/m2 (free parameter)
Rs=0.043068; %Rseries ohm 
Rsh=20.818; %Rshunt ohm
n=1.04; %diode quality factor
q=1.6e-19; %Electron charge coulombs
Kboltz=1.38e-23; %Boltzmann's constant m2 kg s-2 K-1
QAK=q/(n*Kboltz); %Constant including Qcharge,n and Kboltz
Ns=190*6; %number of arrays connected in series (6 cells/module)
Np=12; %number of arrays connected in parallel


T1= 273+25; % T(K)
Voc_T1=215/(50*6); % (V)one cell
Isc_T1=7.7; % (A)
T2=273+45; % T(K)
Voc_T2=202.5/(50*6); % (V)one cell
Isc_T2=7.735;%(A)
Id0=Isc_T1/(exp(QAK*Voc_T1/T1)-1);
K0=(Isc_T2-Isc_T1)/(T2-T1);
ILConst= Isc_T1/G0;

Vg=1.12; % bandgap voltage from Walker
IL=7.7293; 
Voc=4.3;
k=0.814;
Vmp=Voc*k;

% %Matlab script for the PV model SunPower E19/245
% Temp=298; %Temperature K (free parameter)
% G0=1000; %Irradiance w/m2 (free parameter)
% Rs=0.4804/72; %Rseries ohm 
% Rsh=370.7525; %Rshunt ohm
% n=0.92671; %diode quality factor 
% Kboltz=1.38e-23; %Boltzmann's constant m2 kg s-2 K-1
% Qcharge=1.6e-19; %Electron charge coulombs
% QAK=Qcharge/(n*Kboltz); %Constant including Qcharge,n and Kboltz
% Ns=72*18; %number of arrays connected in series (72 cells/array)
% Np=7; %number of arrays connected in parallel
% 
% %SunPower E19/245 module data 
% T1= 273+25; % T(K)
% Voc_T1=48.8/72; % (V)one cell
% Isc_T1=6.43; % (A)
% T2=273+45; % T(K)
% Voc_T2=46.32/72; % (V)one cell
% Isc_T2=6.48;%(A)
% 
% Id0=Isc_T1/(exp(QAK*Voc_T1/T1)-1);
% K0=(Isc_T2-Isc_T1)/(T2-T1);
% Vg=1.12; % bandgap voltage from Walker
% ILConst= Isc_T1/G0; %proportionality constant I/G- [A/(wm-2)]

