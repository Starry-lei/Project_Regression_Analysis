clc;
%data of C60.m in the sequence like Temp,Phi,PhiDot,kf.
%-----------------------------------------CHOICE 03-----------------------------------------------------
mdata=load('C60.m');
%mdata=load('C15.m');
%mdata=load('100Cr6.m');
Temp=mdata(:,1);
Phi=mdata(:,2);
PhiDot=mdata(:,3);
kf=mdata(:,4);
con={'Phi','PhiDot','Temp'};
%-----------------------------------------CHOICE 04-----------------------------------------------------
% Function expressions of C60:
funcC6001=@(Temp,Phi,PhiDot)2336.6042*exp(-0.0021374*Temp)*PhiDot^(0.011358+(8.7041e-05*Temp))*Phi^(0.29178)*exp(-0.66393*Phi);
funcC6002=@(Temp,Phi,PhiDot)2203.6714*exp(-0.0022399*Temp)*PhiDot^(-0.045741+0.00010135*Temp+0.0084669*Phi)*Phi^(-0.045741+-0.00017449*Temp+2.1399e-05*PhiDot)*exp(-0.55383*Phi+0.0035972*PhiDot);
 
% Function expressions of C15:
%funcC6001=@(Temp,Phi,PhiDot)767.0981*exp(-0.00094523*Temp)*PhiDot^(0.052005+(-9.0079e-05*Temp))*Phi^(0.098706)*exp(-0.047643*Phi);
%funcC6002=@(Temp,Phi,PhiDot)965.2317*exp(-0.0014811*Temp)*PhiDot^(0.067929+-8.5815e-05*Temp+-0.011947*Phi)*Phi^(0.3055+-0.00042756*Temp+0.00050364*PhiDot)*exp(0.065144*Phi+-0.00031097*PhiDot);

% Function expressions of 100Cr6:
%funcC6001=@(Temp,Phi,PhiDot)1275.2557*exp(-0.0012855*Temp)*PhiDot^(0.062255+(-0.00010549*Temp))*Phi^(0.15201)*exp(-0.16814*Phi);
%funcC6001=@(Temp,Phi,PhiDot)1377.2266*exp(-0.0015931*Temp)*PhiDot^(0.082002+-0.00010518*Temp+-0.052111*Phi)*Phi^(0.25864+-0.00025809*Temp+0.00082974*PhiDot)*exp(-0.029658*Phi+0.00095732*PhiDot);

%We set Phi amd PhiDot as constants to get 2-di graph

minT=min(Temp);
maxT=max(Temp);
Tempnew=(minT:1:maxT)';
datalen2=length(Tempnew);
Phinew=(ones(size(Tempnew)))';
phicons=0.5;
Phinew=phicons*Phinew';
PhiDotnew=(ones(size(Tempnew)))';
phidotcons=1.5;
PhiDotnew=phidotcons*PhiDotnew';
kf_func2d01=ones(datalen2,1);
kf_func2d02=ones(datalen2,1);
for i=1:length(Tempnew)
kf_func2d01(i)=funcC6001(Tempnew(i,1),Phinew(i,1),PhiDotnew(i,1));
end
for i=1:length(Tempnew)
kf_func2d02(i)=funcC6002(Tempnew(i,1),Phinew(i,1),PhiDotnew(i,1));
end

%Scatter the points on the graphy

mdatan =[Temp,Phi,PhiDot,kf];
ind1=mdatan(:,2)==0.5;
datapoint1=mdatan(ind1,:);
ind2=datapoint1(:,3)==1.5;
datapoint2=datapoint1(ind2,:);

figure(1);
plot(Tempnew,kf_func2d01);
hold on
scatter(datapoint2(:,1),datapoint2(:,4));
grid on
xlabel('Temp'); 
ylabel('kf'); 
s1=['kf-values of model01 for ',char(con(1)),'=',num2str(phicons),';',char(con(2)),'=',num2str(phidotcons)];
title(s1);
hold off

figure(2);
plot(Tempnew,kf_func2d02);
hold on
scatter(datapoint2(:,1),datapoint2(:,4));
grid on
xlabel('Temp'); 
ylabel('kf'); 
s12=['kf-values of model02 for ',char(con(1)),'=',num2str(phicons),';',char(con(2)),'=',num2str(phidotcons)];
title(s12);
hold off

% we set only Phi as constants to get 3-di graph

minPhid=min(PhiDot);
maxPhid=max(PhiDot);
PhiDotnew2=(minPhid:10:maxPhid)';
q=datalen2;
x=linspace(minT,maxT,q);
y=linspace(minPhid,maxPhid,q);
[X,Y]=meshgrid(x,y);
Z1 = zeros(q);
Z2 = zeros(q);
for j=1:q 
    for i=1:q 
        Z1(j,i) =funcC6001(x(i),Phinew(i,1),y(j));
    end
end
for j=1:q 
    for i=1:q 
        Z2(j,i) =funcC6002(x(i),Phinew(i,1),y(j));
    end
end

figure(3);
surf(X,Y,Z1);
hold on
shading interp;
colormap default;
xlabel('Temp'); 
ylabel('PhiDot'); 
zlabel('kf');
s2 = ['kf-values of model 01 for ',char(con(1)),' = ',num2str(phicons)];
title(s2);
grid on
scatter3(datapoint1(:,1),datapoint1(:,3),datapoint1(:,4));
hold off
figure(4);
surf(X,Y,Z2);
hold on
scatter3(datapoint1(:,1),datapoint1(:,3),datapoint1(:,4));
shading interp;
colormap default;
xlabel('Temp'); 
ylabel('PhiDot'); 
zlabel('kf');
s22 = ['kf-values of model 02 for ',char(con(1)),' = ',num2str(phicons)];
title(s22);
grid on
hold off
%-----------------------------------------Hope you have a nice day!:)-----------------------------------------------------



