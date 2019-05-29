function[RA,RV,RD,SA,SV,SD,PSA,PSV]=SDOFR(acc,Tn,d,dt)

%%set structure properties%%
%m=1;
Wn=2*pi/Tn;
Wd=Wn*(1-d^2)^0.5;
f=-1*acc;  %pi:forcing function

%%initial condition%%
u0=0;
ud0=0;
udd0=0;
f0=0;

%%results matrix%%
RA=zeros(1,15000);
RV=zeros(1,15000);
RD=zeros(1,15000);

%%solution of ODE at step 1%%
   g=f0; %ai=p(ti)
   h=(f(1)-f0)/dt; %bi
   
   a0=(g/Wn^2)-(2*d*h/Wn^3);
   a1=h/Wn^2;
   a2=u0-a0;
   a3=(ud0+d*Wn*a2-a1)/Wd;
   u(1)=a0+a1*dt+a2*exp(-d*Wn*dt)*cos(Wd*dt)+a3*exp...
    (-d*Wn*dt)*sin(Wd*dt);
   ud(1)=a1+exp(-d*Wn*dt)*((Wd*a3-d*Wn*a2)*cos(Wd*dt)...
       -(Wd*a2+d*Wn*a3)*sin(Wd*dt));
   udd(1)=exp(-d*Wn*dt)*sin(Wd*dt)*(d*Wn*(Wd*a2+d*Wn*a3)...
       -Wd*(Wd*a3-d*Wn*a2))-exp(-d*Wn*dt)*cos(Wd*dt)...
       *(d*Wn*(Wd*a3-d*Wn*a2)+Wd*(Wd*a2+d*Wn*a3));
   RD(1,1)=u(1);
   RV(1,1)=ud(1);
   RA(1,1)=udd(1);
%solutions of ODE for the folloing steps%
for i=1:(length(acc)-1)  %causion:'acc'start with 1
   %external force%
   g=f(i); %ai=p(ti)
   h=(f(i+1)-f(i))/dt; %bi
   
   a0=(g/Wn^2)-(2*d*h/Wn^2);
   a1=h/Wn^2;
   a2=u(i)-a0;
   a3=(ud(i)+d*Wn*a2-a1)/Wd;
   u(i+1)=a0+a1*dt+a2*exp(-d*Wn*dt)*cos(Wd*dt)+a3*exp...
    (-d*Wn*dt)*sin(Wd*dt);
   ud(i+1)=a1+exp(-d*Wn*dt)*((Wd*a3-d*Wn*a2)*cos(Wd*dt)...
       -(Wd*a2+d*Wn*a3)*sin(Wd*dt));
   udd(i+1)=exp(-d*Wn*dt)*sin(Wd*dt)*(d*Wn*(Wd*a2+d*Wn*a3)...
       -Wd*(Wd*a3-d*Wn*a2))-exp(-d*Wn*dt)*cos(Wd*dt)...
       *(d*Wn*(Wd*a3-d*Wn*a2)+Wd*(Wd*a2+d*Wn*a3));
    
   RD(1,i+1)=u(i+1);
   RV(1,i+1)=ud(i+1);
   RA(1,i+1)=udd(i+1);
   
end
  SA=max(abs(RA))
  SV=max(abs(RV))
  SD=max(abs(RD))
  PSV=SD*Wn
  PSA=SD*Wn^2
    

end