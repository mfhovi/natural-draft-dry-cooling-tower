clear all

%% condenser two
Dc=895.8836111;%kg/s

L1=0.5;%wall thickness(mm)
mw=49127;% kg/s
tw=18;
n=11500*4;
dic=0.024;
doc=0.025;
rhoaw=(1.49343e-3-3.7164e-6*(tw+273.15)+7.09782e-9*(tw+273.15)^2-1.90321e-20*(tw+273.15)^6)^-1;%kg/m3
vw=mw/(rhoaw*n*pi*dic^2/4);
K0= 12.87*vw^5 - 145.1*vw^4 + 638.4*vw^3 - 1525*vw^2 + 3007*vw + 714.7;% velocty unit m/s,K0 unit W/m^2*k;
Fw= 1.151*exp(-((tw-55.95)/75.15)^2) + 0.07026 *exp(-((tw-20.01)/14.15)^2);%water temp 
Fm= (0.002362*L1^2 - 0.2093*L1 + 4.056) / (L1 + 3.863);%wall thickness-mm
Fc=0.9;
Kc=0.001*K0*Fw*Fm*Fc;%KJ/(m^2*K*s）
cpw=4.1868;%KJ/(Kg*K)
A=24769.5*4;
L=A/(pi*doc*n);%  
ts=28.9;%input
H1 = py.CoolProp.CoolProp.PropsSI('H','T',ts+273.15,'Q',1,'water')*.001;
H2 = py.CoolProp.CoolProp.PropsSI('H','T',ts+273.15,'Q',0,'water')*.001;
%H1=refpropm('H','T',ts+273.15,'Q',1,'water')*0.001;%kJ/kg
%H2=refpropm('H', 'T',ts+273.15,'Q',0, 'water')*0.001;
Q=Dc*(H1-H2);%flow rate-kg/s，KJ/(Kg)，KJ/(m^2*K*s），
tw1=Dc*(H1-H2)/(cpw*mw)+tw;
Tm=(tw1-tw)/log((ts-tw)/(ts-tw1));
Qh=Kc*A*Tm;%KJ/(m^2*K*s） m^2 


while abs((Q-Qh)/Q)>0.000001
    if ((Q-Qh)/Q)>0
ts=ts+0.05*abs((Q-Qh)/Q);
    else
        ts=ts-0.05*abs((Q-Qh)/Q);
    end
H1 = py.CoolProp.CoolProp.PropsSI('H','T',ts+273.15,'Q',1,'water')*.001;
H2 = py.CoolProp.CoolProp.PropsSI('H','T',ts+273.15,'Q',0,'water')*.001;
%H1=refpropm('H','T',ts+273.15,'Q',1,'water')*0.001;%
%H2=refpropm('H', 'T',ts+273.15,'Q',0, 'water')*0.001;
Q=Dc*(H1-H2);%flow -kg/s，KJ/(Kg*K)，KJ/(m^2*K*s），
tw1=Dc*(H1-H2)/(cpw*mw)+tw;

Tm=(tw1-tw)/log((ts-tw)/(ts-tw1));
Qh=Kc*A*Tm;%KJ/(m^2*K*s） m^2 

%fprintf("ts=%e\n",ts);
end 
pc = py.CoolProp.CoolProp.PropsSI('P', 'T',ts+273.15, 'Q',0, 'water');%Kpa
%pc=refpropm('P', 'T',ts+273.15, 'Q',0, 'water');%Kpa


