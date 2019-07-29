clear
clc
%-------------------------------------------------CHOICE 01----------------------------------------------------- 
% Import data
%FilePosition=('C:\Users\ASUS\Desktop\CBEM\Project_Regression_Analysis\C60.csv');
%FilePosition=('C:\Users\ASUS\Desktop\CBEM\Project_Regression_Analysis\C15.csv');
FilePosition=('C:\Users\ASUS\Desktop\CBEM\Project_Regression_Analysis\100Cr6.csv');
fid=fopen(FilePosition);
data=textscan(fid,'%f %f %f %f','headerlines',1,'delimiter',';');
fclose(fid);
PhiDot=data{1};
Phi=data{2};
Temp=data{3};
kf=data{4};
DM=[Temp,Phi,PhiDot,kf];
dlmwrite('C60.m',DM);
%-----------------------------------------CHOICE 02-----------------------------------------------------
%dlmwrite('C15.m',DM);
%dlmwrite('100Cr6.m',DM);

%Linearization of model function 1
m=log(kf);
n_1=Temp;
n_2=log(PhiDot);
n_3=log(Phi);
n_4=Phi;
n_5=n_1.*n_2;
%Linearization of model function 2
w_1=Temp;
w_2=log(PhiDot);
w_3=log(Phi);
w_4=Phi;
w_5=w_1.*w_2;
w_6=w_4.*w_2;
w_7=w_1.*w_3;
w_8=PhiDot.*w_3;
w_9=PhiDot;
datalen=length(m);
one=ones(datalen,1);
%A is the measured data from Moodle;
A1=[n_1,n_2,n_3,n_4,n_5,one];
A2=[w_1,w_2,w_3,w_4,w_5,w_6,w_7,w_8,w_9,one];
X=(A1'*A1)\A1'*m;
X2=(A2'*A2)\A2'*m;
%X=inv(A'*A)*A'*m;
%disp(X);
%vn=variable name
n_6=ones(datalen,1);
w_10=ones(datalen,1);
vn={'n_1','n_2','n_3','n_4','n_5','1'};
vn2={'w_1','w_2','w_3','w_4','w_5','w_6','w_7','w_8','w_9','1'};
bfc=@(n_1, n_2 ,n_3, n_4, n_5,n_6)X(1)*n_1+X(2)*n_2+X(3)*n_3+X(4)*n_4+X(5)*n_5+X(6)*n_6;
bfc2=@(w_1, w_2 ,w_3, w_4, w_5,w_6,w_7,w_8, w_9, w_10)X2(1)*w_1+X2(2)*w_2+X2(3)*w_3+X2(4)*w_4+X2(5)*w_5+X2(6)*w_6+X2(7)*w_7+X2(8)*w_8+X2(9)*w_9+X2(10)*w_10;
fdata=ones(datalen,1);
   for i=1:datalen
       fdata(i)=bfc(A1(i,1),A1(i,2),A1(i,3),A1(i,4),A1(i,5),A1(i,6));
   end
   fdata2=ones(1188,1);
   for i=1:datalen
       fdata2(i)=bfc2(A2(i,1),A2(i,2),A2(i,3),A2(i,4),A2(i,5),A2(i,6),A2(i,7),A2(i,8),A2(i,9),A2(i,10));
   end
    A1=[A1,m];
    A2=[A2,m];
    r201=RSquared(A1,fdata);
    r202=RSquared(A2,fdata2);
% Format the output for function01

C='';
F='(';  
for i=1:length(X)
    v1=num2str(i);
    v2=num2str(X(i));
    if i==1
        C=strcat('x',v1,'=',v2);
        F=strcat(F,v2,')*',char(vn(i)));
    else
        
        C=strcat(C,';x',v1,'=',v2);
        F=strcat(F,'+(',v2,')*',char(vn(i)));
    end
end
% Format the output for function02
C22='';
F22='(';
for i=1:length(X2)
    v12=num2str(i);
    v22=num2str(X2(i));
    if i==1
        C22=strcat('x',v12,'=',v22);
        F22=strcat(F22,v22,')*',char(vn2(i)));
    else
        
        C22=strcat(C22,';x',v12,'=',v22);
        F22=strcat(F22,'+(',v22,')*',char(vn2(i)));
    end
end
%Represent the best fitting curve as inline function
disp(['Curve Fitting-Start->',char(datestr(now))]);
fprintf('\n');
    disp(['Parameter:',C]);
    fprintf('\n');
    disp(['Best fitting curve:m(n1,n2,n3,n4)=',F]);
   %output Rsquared
    v=num2str(r201);
    fprintf('\n');
    disp(['Rsquared:',v]);
   fprintf('\n');
    X(6)=exp(X(6));
    disp('Transformed parameters for the non linear model function:');
    v1=num2str(X(1));
    v2= num2str(X(2));
    v3 = num2str(X(3));
    v4 = num2str(X(4));
    v5 = num2str(X(5));
    v6 = num2str(X(6));
fprintf('x1=%s\n',v1);
fprintf('x2=%s\n',v2);
fprintf('x3=%s\n',v3);
fprintf('x4=%s\n',v4);
fprintf('x5=%s\n',v5);
fprintf('x6=%s\n',v6);
fprintf('\n');
    disp('Non linear function 01:');
    C2= 'kf = g(Temp,Phi,PhiDot) = ';
    F2=strcat(v6,'*exp(',v1,'*Temp)*PhiDot^(',v2,'+(',v5,'*Temp))*Phi^(',v3,')*exp(',v4,'*Phi)');
    disp(strcat(C2,F2));

   %OUTPUT for function model 02
   
   fprintf('\n');
    disp(['Parameter:',C22]);
    fprintf('\n');
    disp(['Best fitting curve:m(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10)=',F22]);
   %output Rsquared
    v2=num2str(r202);
    fprintf('\n');
    disp(['Rsquared:',v2]);
   fprintf('\n');
    X2(10)=exp(X2(10));
    disp('Transformed parameters for the non linear model function 02:');
    v12=num2str(X2(1));
    v22= num2str(X2(2));
    v32 = num2str(X2(3));
    v42 = num2str(X2(4));
    v52 = num2str(X2(5));
    v62 = num2str(X2(6));
    v72 = num2str(X2(7));
    v82 = num2str(X2(8));
    v92 = num2str(X2(9));
    v102 = num2str(X2(10));
fprintf('x1=%s\n',v12);
fprintf('x2=%s\n',v22);
fprintf('x3=%s\n',v32);
fprintf('x4=%s\n',v42);
fprintf('x5=%s\n',v52);
fprintf('x6=%s\n',v62);
fprintf('x7=%s\n',v72);
fprintf('x8=%s\n',v82);
fprintf('x9=%s\n',v92);
fprintf('x10=%s\n',v102);
fprintf('\n');
    disp('Non linear function 02:');
    C202= 'kf = g(Temp,Phi,PhiDot) = ';
    F202=strcat(v102,'*exp(',v12,'*Temp)*PhiDot^(',v22,'+',v52,'*Temp','+',v62,'*Phi',')','*Phi^(',v32,'+',v72,'*Temp','+',v82,'*PhiDot)*exp(',v42,'*Phi+',v92,'*PhiDot)');
    disp(strcat(C202,F202));