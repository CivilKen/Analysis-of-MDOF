function[RA,RV,RD,SA,SV,SD]=MDOFR(acc,Wk,d,dt,pf,fi,M,u0,ud0)

%%set structure properties%%
m=1;
Wn=Wk;
Wd=Wn*sqrt(1-d^2);
f=-1*pf*acc;  %forcing function

%%initial condition%%
z0=(fi'*M*u0)/(fi'*M*fi)
zd0=(fi'*M*ud0)/(fi'*M*fi)
f0=0;

%%results matrix%%
RA=zeros(1,15000);
RV=zeros(1,15000);
RD=zeros(1,15000);

%%solution of ODE at step 1%%
   g0=f0; %ai=p(0)
   h0=(f(1)-f0)/dt; %b0
   
   A0=(g0/Wn^2)-(2*d*h0/Wn^3);
   A1=h0/Wn^2;
   A2=z0-A0;
   A3=(zd0+d*Wn*A2-A1)/Wd;
   z(1)=A0+A1*dt+A2*exp(-d*Wn*dt)*cos(Wd*dt)+A3*exp(-d*Wn*dt)*sin(Wd*dt);
   zd(1)=A1+exp(-d*Wn*dt)*((Wd*A3-d*Wn*A2)*cos(Wd*dt)-(Wd*A2+d*Wn*A3)*...
       sin(Wd*dt));
   zdd(1)=exp(-d*Wn*dt)*sin(Wd*dt)*(d*Wn*(Wd*A2+d*Wn*A3)-...
       Wd*(Wd*A3-d*Wn*A2))-exp(-d*Wn*dt)*cos(Wd*dt)*(d*Wn*(Wd*A3-d*Wn*A2)...
       +Wd*(Wd*A2+d*Wn*A3));
   RD(1,1)=z(1);
   RV(1,1)=zd(1);
   RA(1,1)=zdd(1);
%solutions of ODE for the folloing steps%
for i=1:(length(acc)-1)  %causion:'acc'start with 1
   %external force%
   gi=f(i); %ai=p(0)
   hi=(f(i+1)-f(i))/dt; %b0
   
   A0=(gi/Wn^2)-(2*d*hi/Wn^3);
   A1=hi/Wn^2;
   A2=z(i)-A0;
   A3=(zd(i)+d*Wn*A2-A1)/Wd;
   z(i+1)=A0+A1*dt+A2*exp(-d*Wn*dt)*cos(Wd*dt)+A3*exp(-d*Wn*dt)*sin(Wd*dt);
   zd(i+1)=A1+exp(-d*Wn*dt)*((Wd*A3-d*Wn*A2)*cos(Wd*dt)-(Wd*A2+d*Wn*A3)*...
       sin(Wd*dt));
   zdd(i+1)=exp(-d*Wn*dt)*sin(Wd*dt)*(d*Wn*(Wd*A2+d*Wn*A3)-...
       Wd*(Wd*A3-d*Wn*A2))-exp(-d*Wn*dt)*cos(Wd*dt)*(d*Wn*(Wd*A3-d*Wn*A2)...
       +Wd*(Wd*A2+d*Wn*A3));;
    
   RD(1,i+1)=z(i+1);
   RV(1,i+1)=zd(i+1);
   RA(1,i+1)=zdd(i+1);
   
end
  SA=max(abs(RA))
  SV=max(abs(RV))
  SD=max(abs(RD))
 
end