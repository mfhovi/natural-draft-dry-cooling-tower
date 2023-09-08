tic
clear all
%% ambient conditions
Ta1=-2.7+273.15;
%Ta1= 5+273.15; 
%Ta1= 10+273.15;
%Ta1= 15+273.15;
%Ta1= 20+273.15;
%Ta1= 25+273.15;
%Ta1= 30+273.15;

Twb1= -10+273.15;
pa1=102841;
%% tower parameter
h4=26; %inlet height of dry section
d4=120.53;

h3=h4/2;
d3=d4;

h5=150; %tower height
d5=84.53;
A5=pi*d5^2/4;

nts=80;% number of tower supports
Lts=24;%length of tower support,m
dts=0.5;%diameter of tower support,m
Cdts=2.0;%drag coeffient of support
ts=0.8;% thickness of the square-edged(90°)shell at the inlet to the tower,m
ri_d3=0.02;%the tower shell rounded inlet 
%% dry section
d=0.025;%hydraulic diameter of tube,m
relativeroughness=5.24*10^(-4);
Ats=3.664*10^(-4);%inside cross-sectional flow area,m^2
Lt=24;%length of finned tube,m
Lte=22;%effective length of finned tube,ineffective length due to obstructions on air side,m
Ati=pi*d;%inside area of the tube per unit length,m
nr=4;%number of tube rows
ntb=320;%number of tubes per bundle
nwp=2;%number of water passes
nb=176;%number of deltas
Afrd=Lte*2.4*2*nb; %total effective frontal area of bundles,m^2
aj=49.08/2; %2aj=61.5°,apex angle of a-frame
ratiomf=0.433;%ratio of minimum to free stream flow area
Kci=0.05;%inlet contraction loss coefficient

a11=-0.605;
a12=4.34;
a13=-9.72;
a14=7.54;
a21=0.0231;
a22=0.0059;
a23=-0.248;
a24=0.287;
a31=0.294;
a32=-1.99;
a33=4.32;
a34=-3;
a41=0.0198;
a42=-0.305;
a43=0.897;
a44=-0.731;


%% operation parameter
R=8.3144598/28.966*1e3;
g=9.81;
Rv=461.52;

cpa1=1.045356e3-3.161783e-1*(Ta1+273.15)/2+7.083814e-4*((Ta1+273.15)/2)^2-2.705209e-7*((Ta1+273.15)/2)^3; % the specific heat of air inlet tower
cpv1=1.3605e3+2.31334*(Ta1+273.15)/2-2.46784e-10*((Ta1+273.15)/2)^5+5.91332e-13*((Ta1+273.15)/2)^6; % the specific heat of vapor after drift elimination
cpw1=8.15599e3-2.80627e1*(Ta1+273.15)/2+5.11283e-2*((Ta1+273.15)/2)^2-2.17582e-13*((Ta1+273.15)/2)^6; % the specific heat of water after drift elimination

Ta2=Ta1-0.00975*h3;
pa2=pa1*(1-(0.00975*h3/Ta1))^3.5;
rhoa2=pa2/(R*Ta2);
cpa2=1.045356e3-3.161783e-1*(Ta2+273.15)/2+7.083814e-4*((Ta2+273.15)/2)^2-2.705209e-7*((Ta2+273.15)/2)^3; % the specific heat of air inlet tower
cpv2=1.3605e3+2.31334*(Ta2+273.15)/2-2.46784e-10*((Ta2+273.15)/2)^5+5.91332e-13*((Ta2+273.15)/2)^6; % the specific heat of vapor after drift elimination
cpw2=8.15599e3-2.80627e1*(Ta2+273.15)/2+5.11283e-2*((Ta2+273.15)/2)^2-2.17582e-13*((Ta2+273.15)/2)^6; % the specific heat of water after drift elimination

pa6=pa1*(1-(0.00975*h5/Ta1))^3.5;
Ta6=Ta1-0.00975*h5;

cpa6=1.045356e3-3.161783e-1*(Ta6+273.15)/2+7.083814e-4*((Ta6+273.15)/2)^2-2.705209e-7*((Ta6+273.15)/2)^3; % the specific heat of air inlet tower
cpv6=1.3605e3+2.31334*(Ta6+273.15)/2-2.46784e-10*((Ta6+273.15)/2)^5+5.91332e-13*((Ta6+273.15)/2)^6; % the specific heat of vapor after drift elimination
cpw6=8.15599e3-2.80627e1*(Ta6+273.15)/2+5.11283e-2*((Ta6+273.15)/2)^2-2.17582e-13*((Ta6+273.15)/2)^6; % the specific heat of water after drift elimination

rhoa6=pa1/(R*Ta6);


%% mass flow rate distribution

%mw= (0.9)*64000*1000/3600;%kg/s ; 10% less than base condition 
%mw= (0.95)*64000*1000/3600;%kg/s ; 5% less than base condition
mw=64000*1000/3600;%kg/s ;base condition
%mw= (1.05)*64000*1000/3600;%kg/s ; 5% more than base condition
%mw= (1.10)*64000*1000/3600;%kg/s ; 10% more than base condition

Tin=313.5475;
%% integration initial
pa3=pa1*(1-(0.00975*h3/Ta1))^3.5-15;
v02=1.55;
Twod=(Tin*3+Ta6)/4; 
i=1;
p=1;
m=1;
n=1;
q=1;
delta_pf=0.51;
while abs(delta_pf)>=0.5 % Outer Loop
if delta_pf>0
    pa3=pa3+min(0.5,0.5*delta_pf);
else
    pa3=pa3-min(0.5,0.5*abs(delta_pf));
end

delta_fl_bad=0.21;
while abs(delta_fl_bad)>=0.2 % Dry cooling section
if delta_fl_bad>0
    v02=v02+min(0.002,abs(delta_fl_bad));
else
    v02=v02-min(0.002,abs(delta_fl_bad));
end

ma23=v02*Afrd*rhoa2;

delta_en_bad=0.001;
while abs(delta_en_bad)>=1e-5

if delta_en_bad>0
    Twod=Twod+min(0.01,0.01*abs(delta_en_bad));
else
    Twod=Twod-min(0.01,0.01*abs(delta_en_bad));
end
Twdm=(Tin+Twod)/2;
cpwdm=8.15599e3-2.80627e1*(Twdm+273.15)/2+5.11283e-2*((Twdm+273.15)/2)^2-2.17582e-13*((Twdm+273.15)/2)^6; % the specific heat of water after drift elimination
Qwd=mw*cpwdm*(Tin-Twod);
Ta3=Qwd/(ma23*cpa2)+Ta2;

niuwd=2.414*1e-5*10^(247.8/((Tin+Twod)/2-140));%water kinematic viscosity
Lamdawd=-6.14255e-1+6.9963e-3*((Tin+Twod)/2)-1.01075e-5*((Tin+Twod)/2)^2+4.74737e-12*((Tin+Twod)/2)^4;%water thermal conductivity
Prwd=niuwd*cpwdm/Lamdawd;
Rhoaw=1/(1.49343*1e-3-3.7164*1e-6*((Tin+Twod)/2)+7.09782*1e-9*((Tin+Twod)/2)^2-1.90321*1e-20*((Tin+Twod)/2)^6);
Rewd=mw*nwp*d/Ats/ntb/nb/niuwd;
f=power(1.82*log10(Rewd)-1.64,-2);%Water side drag coefficient
Nuwd=(f/8)*(Rewd-1000)*Prwd*(1+power(d/Lte,2/3))/(1+12.7*power(f/8,0.5)*(power(Prwd,2/3)-1));%water side reynolds no
hw=Nuwd*Lamdawd/d;%Heat transfer coefficient at current water side
ha=45.448*power(v02,0.327);
haeAa=ha*Afrd*31.93/0.3044;
Aw=Ati*Lte*ntb*nb;
HA=1/(1/(haeAa)+1/(hw*Aw));

 delta_t=((Tin-Ta3)-(Twod-Ta2))/log(abs((Tin-Ta3)/(Twod-Ta2)));%the logarithmic mean temperature difference
 s1=(Tin-Twod)/(Tin-Ta2);
 s2=(Ta3-Ta2)/(Tin-Ta2);
 s12=s1/s2;
 s3=abs((s1-s2)/log(abs((1-s2)/(1-s1))));
 arctan=atan(s12);
 Ft= 1-a11*(1-s3)*sin(2*arctan)-a12*(1-s3)^2*sin(2*arctan)-a13*(1-s3)^3*sin(2*arctan)-...
     a14*(1-s3)^4*sin(2*arctan)-a21*(1-s3)*sin(4*arctan)-a22*(1-s3)^2*sin(4*arctan)-...
     a23*(1-s3)^3*sin(4*arctan)-a24*(1-s3)^4*sin(4*arctan)-a31*(1-s3)*sin(6*arctan)-...
     a32*(1-s3)^2*sin(6*arctan)-a33*(1-s3)^3*sin(6*arctan)-a34*(1-s3)^4*sin(6*arctan)-...
     a41*(1-s3)*sin(8*arctan)-a42*(1-s3)^2*sin(8*arctan)-a43*(1-s3)^3*sin(8*arctan)-...
     a44*(1-s3)^4*sin(8*arctan);%temperature correction factor for a two-pass.crossflow heat exchanger 

Qhd=Ft*HA*delta_t;
delta_q=Qwd-Qhd;
delta_en_bad=(Qwd-Qhd)/Qwd;
end
Ta23=(Ta2+Ta3)/2;
rhoa3=pa3/(R*Ta3);
rhoa23=2/(1/rhoa2+1/rhoa3);
va3=ma23/rhoa3/Afrd;
Hd=pa2-pa3;

 niana34=2.287973*10^(-6)+6.259793*10^(-8)*Ta23-3.131956*10^(-11)*Ta23^2+8.15038*10^(-15)*Ta23^3;%viscosity of Ta34,kg/sm
 
 Khe=51.756-32.614*v02+11.682*v02^2-2.0986*v02^3+0.1808*v02^4-0.0059*v02^5;
 ajm=0.0019*aj^2+0.9133*aj-3.1558;% mean inlet flow angle
 Kd=exp(5.488405-0.2131209*aj+3.533265*10^(-3)*aj^2-0.2901016*10^(-4)*aj^3);% the downstream loss coefficient
 Kheaj=Khe+2*rhoa3*(1/sin(ajm/180*pi)-1)/(rhoa2+rhoa3)*((1/sin(ajm/180*pi)-1)+2*Kci^0.5)+2*rhoa2*Kd/(rhoa2+rhoa3);%for the bundle loss coefficient
 Zd=Kheaj*(ma23/Afrd)^2/(2*rhoa23)+rhoa3*va3^2/2;

delta_fl_bad=Hd-Zd
end

Ta5=Ta3-0.00975*(h5-h3);

delta_pa5=0.2;
pa5=pa3-30;
while abs(delta_pa5)>=0.1
    if delta_pa5>0
         pa5=pa5-abs(delta_pa5);
    else pa5=pa5+abs(delta_pa5);
    end

rhoa5=pa5/(R*Ta5);

Frd=(ma23/A5)^2/rhoa5/(rhoa6-rhoa5)/g/d5;
Kto=-0.129*(Frd*d5/d3)^(-1)+0.0144*(Frd*d5/d3)^(-1.5);
ae5=1.004+5.8*(d5/d3)^9+(0.007+0.043*(d5/d3)^2.5)*Frd^(-1.5);
pa5c=pa6+(Kto+ae5)*(ma23/A5)^2/rhoa5/2;
delta_pa5=pa5-pa5c;
end
va5=ma23/(rhoa5*A5);
 
delta_p35=pa3-pa3*(1-(0.00975*(h5-h3)/Ta3))^3.5;
pa3m=pa5+(rhoa5*va5^2)/2-(rhoa3*va3^2)/2+delta_p35;

delta_pf=pa3m-pa3
pf(i)=delta_pf;
i=i+1;
end
Qa=ma23*cpa2*(Ta3-Ta2);
Qw=mw*cpwdm*(Tin-Twod);

toc

Ans=[ma23 Qa Qw Twod]
Ans1=[Hd Ta3-Ta2]
[Ans Ans1]
