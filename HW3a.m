clc
clear all
close all

% read earthquake data 
fid = fopen('TCU076.V2','r');
tline = fgetl(fid);
col = 1; % stand for which row to fill in
allData = zeros(9,15000); %所有資料存到all data裡面(使用cell指令)
while ischar(tline)
   
   matchesa1 = strfind(tline,'points of accel'); 
   numa1 = length(matchesa1); %判斷有無符合accl之字串
   if numa1 > 0 %有讀到accl之字串
      %fprintf(1,'%d:%s\n',numa,tline);
      
      A=fscanf(fid,'%f');
      allData(col,1:15001) = A;
      col=col+1;
   end
      
   matchesv1 = strfind(tline,'points of veloc');
   numv1 = length(matchesv1); %判斷有無符合veloc之字串
   if numv1 > 0 %有讀到veloc之字串
      %fprintf(1,'%d:%s\n',numv,tline);
      
      V=fscanf(fid,'%f');
      allData(col,1:15001) = V;
      col=col+1;
   end
      
   matchesd1 = strfind(tline,'points of displ');
   numd1 = length(matchesd1); %判斷有無符合displ之字串
   if numd1 > 0 %有讀到displ之字串
      %fprintf(1,'%d:%s\n',numd,tline);
      
      D=fscanf(fid,'%f');
      allData(col,1:15000) = D;
      col=col+1;
   end   
   
   
   tline = fgetl(fid); %再讀下一行  
   
end
fclose(fid);

t = 0:0.01:149.99;
A1 = allData(1,1:15000);
V1 = allData(2,1:15000);
D1 = allData(3,1:15000);
A2 = allData(4,1:15000);
V2 = allData(5,1:15000);
D2 = allData(6,1:15000);
A3 = allData(7,1:15000);
V3 = allData(8,1:15000);
D3 = allData(9,1:15000);

% [第一題]SDOF system properies
mass=[12000 15000 20000];
M=diag(mass); % mass matrix
K=[2700000 -1200000 0;-1200000 2200000 -1000000;0 -1000000 1000000]*2; %stiffness of the system

% [第二題]
%set modal shape fi
%set modal frequency W
% v:eigen vector ;d:eigenvalue
[v,d]=eig(K,M); %Generalized Eigenvalues (K-W^2*M)=0
W1=d(1,1)^0.5;
W2=d(2,2)^0.5;
W3=d(3,3)^0.5;
W=[W1 W2 W3]
%3X1 matrix %3 different modes
%every mode has 3 floors
fi1=v(:,1); % mode 1
fi2=v(:,2); % mode 2
fi3=v(:,3); % mode 3
fi=v; %set modal shape matrix fi (3X3 matrix)
figure %graphics of fi-value if different modes
subplot(1,3,1)
plot([0,fi1(1,1),fi1(2,1),fi1(3,1)],[0,1,2,3])
title('mode1')
subplot(1,3,2)
plot([0,fi2(1,1),fi2(2,1),fi2(3,1)],[0,1,2,3])
title('mode2')
subplot(1,3,3)
plot([0,fi3(1,1),fi3(2,1),fi3(3,1)],[0,1,2,3])
title('mode3')


% [第三題]cauculate participation factor for every mode
one=[1;1;1];
Lk1=fi1'*M*one;
Mk1=fi1'*M*fi1;
PF1=Lk1/Mk1;
Lk2=fi2'*M*one;
Mk2=fi2'*M*fi2;
PF2=Lk2/Mk2;
Lk3=fi3'*M*one;
Mk3=fi3'*M*fi3;
PF3=Lk3/Mk3;
PF=[PF1 PF2 PF3]; %store PF value into a matrix

% [第四題]solution of each degree of freedom
RAM=zeros(3,15000);%store the data of every degree of freedom in the matrix
RVM=zeros(3,15000);%this value is z-value, not real response
RDM=zeros(3,15000);%to obtain real response, multiply z-value with fi
%initial condition of the system
u0=zeros(3,1); 
ud0=zeros(3,1);
for k=1:3 % 'k' represents the number of modes
    [RA,RV,RD,SA,SV,SD]=MDOFR(A3,W(k),0.05,0.01,PF(k),fi(:,k),M,u0,ud0)
    
    RAM(k,:)=RA; %zdd for every mode
    RVM(k,:)=RV; %zd for every mode
    RDM(k,:)=RD; %z for every mode
end

%response of each floor
% 1st floor
udd1=fi1(1)*RAM(1,:)+fi2(1)*RAM(2,:)+fi3(1)*RAM(3,:);
ud1=fi1(1)*RVM(1,:)+fi2(1)*RVM(2,:)+fi3(1)*RVM(3,:);
u1=fi1(1)*RDM(1,:)+fi2(1)*RDM(2,:)+fi3(1)*RDM(3,:);
% 2nd floor
udd2=fi1(2)*RAM(1,:)+fi2(2)*RAM(2,:)+fi3(2)*RAM(3,:);
ud2=fi1(2)*RVM(1,:)+fi2(2)*RVM(2,:)+fi3(2)*RVM(3,:);
u2=fi1(2)*RDM(1,:)+fi2(2)*RDM(2,:)+fi3(2)*RDM(3,:);
% 3rd floor
udd3=fi1(3)*RAM(1,:)+fi2(3)*RAM(2,:)+fi3(3)*RAM(3,:);
ud3=fi1(3)*RVM(1,:)+fi2(3)*RVM(2,:)+fi3(3)*RVM(3,:);
u3=fi1(3)*RDM(1,:)+fi2(3)*RDM(2,:)+fi3(3)*RDM(3,:);
% all response
U=[u1;u2;u3];
UD=[ud1;ud2;ud3];
UDD=[udd1;udd2;udd3];

t = 0:0.01:149.99;
figure
subplot(3,1,1)
plot(t,udd1)
title('1st floor')
xlabel('time')
ylabel('acceleration')
subplot(3,1,2)
plot(t,ud1)
xlabel('time')
ylabel('velocity')
subplot(3,1,3)
plot(t,u1)
xlabel('time')
ylabel('displacement')

figure
subplot(3,1,1)
plot(t,udd2)
title('2nd floor')
xlabel('time')
ylabel('acceleration')
subplot(3,1,2)
plot(t,ud2)
xlabel('time')
ylabel('velocity')
subplot(3,1,3)
plot(t,u2)
xlabel('time')
ylabel('displacement')

figure
subplot(3,1,1)
plot(t,udd3)
title('3rd floor')
xlabel('time')
ylabel('acceleration')
subplot(3,1,2)
plot(t,ud3)
xlabel('time')
ylabel('velocity')
subplot(3,1,3)
plot(t,u3)
xlabel('time')
ylabel('displacement')

%[第五題]equivalent lateral force
  F=K*U % lateral force of each floor
  f1=F(1,:); % for 1st floor
  F1=max(f1); % maximum lateral force for 1st floor
  f2=F(2,:); % for 2nd floor
  F2=max(f2); % maximum lateral force for 2nd floor
  f3=F(3,:); % for 3rd floor
  F3=max(f3); % maximum lateral force for 3rd floor
  figure
  plot(t,f1)
  title('equivalent lateral force at 1st floor')
  xlabel('time')
  ylabel('force')
  figure
  plot(t,f2)
  title('equivalent lateral force at 2nd floor')
  xlabel('time')
  ylabel('force')
  figure
  plot(t,f3)
  title('equivalent lateral force at 3rd floor')
  xlabel('time')
  ylabel('force')
  

% [第六題]
% Calculate response spectrum from SDOF
% [RAS,RVS,RDS,SAS,SVS,SDS,PSAS,PSVS]=SDOFR(A3,Wk,d,0.01) 
  % for mode 1 (W1)
  [RAS1,RVS1,RDS1,SAS1,SVS1,SDS1,PSAS1,PSVS1]=SDOFR(A3,W1,0.05,0.01);
  % Sa,1=SAS1
   u11=fi1(1)*PF1*SAS1/W1^2; %floor 1
   u12=fi1(2)*PF1*SAS1/W1^2; %floor 2
   u13=fi1(3)*PF1*SAS1/W1^2; %floor 3
  % for mode 2(W2)
  [RAS2,RVS2,RDS2,SAS2,SVS2,SDS2,PSAS2,PSVS2]=SDOFR(A3,W2,0.05,0.01);
   u21=fi2(1)*PF2*SAS2/W2^2; %floor 1
   u22=fi2(2)*PF2*SAS2/W2^2; %floor 2
   u23=fi2(3)*PF2*SAS2/W2^2; %floor 3
  % for mode 3(W3)
  [RAS3,RVS3,RDS3,SAS3,SVS3,SDS3,PSAS3,PSVS3]=SDOFR(A3,W3,0.05,0.01);
   u31=fi3(1)*PF3*SAS3/W3^2; %floor 1
   u32=fi3(2)*PF3*SAS3/W3^2; %floor 2
   u33=fi3(3)*PF3*SAS3/W3^2; %floor 3
   
   U1=[u11,u21,u31]; % 1st floor
   U2=[u12,u22,u32]; % 2nd floor
   U3=[u13,u23,u33]; % 3rd floor
   
   u1m=max(U1); % maximum displacement for 1st floor
   f1m=u1m*12000; % maximum lateral force for 1st floor
   u2m=max(U2); % maximum displacement for 1st floor
   f2m=u2m*15000; % maximum lateral force for 1st floor
   u3m=max(U3); % maximum displacement for 3rd floor
   f3m=u3m*20000; % maximum lateral force for 3rd floor