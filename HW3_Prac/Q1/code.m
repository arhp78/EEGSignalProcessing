%Amirreza Hatamipour
% 97101507
%% read data
clc;clear;
EEG = load("Ex1_data.mat");
Fs=100;
X_org=EEG.X_org;
X1=EEG.X1;
X2=EEG.X2;
X3=EEG.X3;
X4=EEG.X4;
T1=EEG.T1;
T2=EEG.T2;
plotEEG(X_org,'original')
%% GEVD
% section a
clc; close all;
X_org_shift=X_org(:,401:end);
X_org_1=X_org(:,1:end-400);
Px_n=X_org_1*X_org_shift.'/length(X_org_1);
Px=(Px_n+Px_n.');

Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S1_est=U(:,1).'*X_org;
S_est=[S1_est; zeros(7,length(S1_est))];
X1_est=inv(U.')*S_est;
plotEEG(X1_est,'estimated X1')
plotEEG(X1,'original X1')
RRMSE_X1 = RRMSE(X1,X1_est)

figure
t=0:0.01:100-0.01;
plot(t,S1_est)
hold on
plot(t,U(:,1).'*X1)
hold on
xlabel('time(s)')
legend('estimates','original')
title('S1')
%% section b
clc; close all;
min_RRMSE=1;
T=0;
for i=300:10:700
X_org_shift=X_org(:,i+1:end);
X_org_1=X_org(:,1:end-i);
Px_n=X_org_1*X_org_shift.'/length(X_org_1);
Px=(Px_n+Px_n.')/2;

Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S1_est=U(:,1).'*X_org;
S_est=[S1_est; zeros(7,length(S1_est))];
X1_est=inv(U.')*S_est;
RRMSE_X1 = RRMSE(X1,X1_est);
if RRMSE_X1 < min_RRMSE
    min_RRMSE=RRMSE_X1;
    T=i;
end
end
for i=T-15:T+15
        X_org_shift=X_org(:,i+1:end);
X_org_1=X_org(:,1:end-i);
Px_n=X_org_1*X_org_shift.'/length(X_org_1);
Px=(Px_n+Px_n.')/2;

Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S1_est=U(:,1).'*X_org;
S_est=[S1_est; zeros(7,length(S1_est))];
X1_est=inv(U.')*S_est;
RRMSE_X1 = RRMSE(X1,X1_est);
if RRMSE_X1 < min_RRMSE
    min_RRMSE=RRMSE_X1;
    T=i;
end


end
% now found T
        X_org_shift=X_org(:,T+1:end);
X_org_1=X_org(:,1:end-T);
Px_n=X_org_1*X_org_shift.'/length(X_org_1);
Px=(Px_n+Px_n.')/2;

Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S1_est=U(:,1).'*X_org;
S_est=[S1_est; zeros(7,length(S1_est))];
X1_est=inv(U.')*S_est;
RRMSE_X1 = RRMSE(X1,X1_est)
plotEEG(X1_est,'estimated X1 with T=492')
plotEEG(X1,'original X1')


figure
t=0:0.01:100-0.01;
plot(t,S1_est)
hold on
plot(t,U(:,1).'*X1)
hold on
xlabel('time(s)')
legend('estimates','original')
title('S1')
%%
max_D(1,1)=0;
for i=300:1:700
    X_org_shift=X_org(:,i+1:end);
X_org_1=X_org(:,1:end-i);
Px_n=X_org_1*X_org_shift.'/length(X_org_1);
Px=(Px_n+Px_n.')/2;

Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
if D(1,1) > max_D(1,1)
    max_D=D;
    T=i;
end
end
X_org_shift=X_org(:,T+1:end);
X_org_1=X_org(:,1:end-T);
Px_n=X_org_1*X_org_shift.'/length(X_org_1);
Px=(Px_n+Px_n.')/2;
Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S1_est=U(:,1).'*X_org;
S_est=[S1_est; zeros(7,length(S1_est))];
X1_est=inv(U.')*S_est;
RRMSE_X1 = RRMSE(X1,X1_est)
plotEEG(X1_est,'estimated X1 with T=401')
plotEEG(X1,'original X1')
%% section c
clc;close all;
Ton=find(T1);
X_org_Ton=X_org(:,Ton(1,:));
Px=X_org_Ton*X_org_Ton.'/length(X_org_Ton);

Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S2_est=U(:,1).'*X_org;
S_est=[S2_est; zeros(7,length(S2_est))];
X2_est=inv(U.')*S_est;
RRMSE_X2 = RRMSE(X2,X2_est)
plotEEG(X2_est,'estimated X2')
plotEEG(X2,'original X2')

figure
t=0:0.01:100-0.01;
plot(t,S2_est)
hold on
plot(t,U(:,1).'*X2)
hold on
xlabel('time(s)')
legend('estimates','original')
title('S2')
%% section d
clc;close all;
Ton=find(T2);
X_org_Ton=X_org(:,Ton(1,:));
Px=X_org_Ton*X_org_Ton.'/length(X_org_Ton);
Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S2_est=U(:,1).'*X_org;
S_est=[S2_est; zeros(7,length(S2_est))];
X2_est=inv(U.')*S_est;
RRMSE_X2_before = RRMSE(X2,X2_est)
plotEEG(X2_est,'estimated X2 with T2')
plotEEG(X2,'original X2')
figure
t=0:0.01:100-0.01;
plot(t,S2_est)
hold on
plot(t,U(:,1).'*X2)
xlabel('time(s)')
legend('estimates','original')
title('S2 before iteration')

T2_new=abs(S2_est)>0.05;
Ton_new=find(T2_new);
for i=1:100

X_org_Ton_new=X_org(:,Ton_new(1,:));
Px=X_org_Ton_new*X_org_Ton_new.'/length(X_org_Ton_new);

Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S2_est=U(:,1).'*X_org;
T2_new=abs(S2_est)>0.05;
Ton_new=find(T2_new);

end
X_org_Ton_new=X_org(:,Ton_new(1,:));
Px_new=X_org_Ton_new*X_org_Ton_new.'/length(X_org_Ton_new);
Cx_new=X_org*X_org.'/length(X_org);
[U_new,D_new]=eig(Px_new,Cx_new);
D_new= diag(D_new);
[D, index] = sort(D_new, 'descend');
U_new = (U_new(:,index));
S2_est_new=U_new(:,1).'*X_org;
S_est_new=[S2_est_new; zeros(7,length(S2_est_new))];
X2_est_new=inv(U.')*S_est_new;


t=0:0.01:100-0.01;
plot(t,S2_est_new)
hold on
plot(t,U(:,1).'*X2)
xlabel('time(s)')
legend('estimates','original')
title('S2 after iteration')

RRMSE_X2_after = RRMSE(X2,X2_est_new)
plotEEG(X2_est_new,'estimated X2 ')
plotEEG(X2,'original X2')
%% section e
clc;close all;
X_fft=fft(X_org.');
X_fft1=X_fft.';
% freq 10-15 Hz
theta1=1000:1:1500;
theta2=8500:1:9000;
theta=[theta1 theta2];

a=X_fft1(:,theta(1,:));
b=conj(X_fft1(:,theta(1,:)));
Sx=a*b.'/length(theta);

Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Sx,Cx);
D = diag(D);
D=real(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S3_est=U(:,1).'*X_org;
S3_est=real(S3_est);
S_est=[S3_est; zeros(7,length(S3_est))];
X3_est=inv(U.')*S_est;
X3_est=real(X3_est);
plotEEG(X3_est,'estimated X3')
plotEEG(X3,'original X3')
RRMSE_X3 = RRMSE(X3,X3_est)

figure
subplot(1,2,1)
t=0:0.01:100-0.01;
plot(t,S3_est)
hold on
plot(t,U(:,1).'*X3)
hold on
xlabel('time(s)')
legend('estimates','original')
title('S3')
subplot(1,2,2)
f=0:0.01:100-0.01;
plot(t,abs(fft(S3_est)))
hold on
plot(t,abs(fft(U(:,1).'*X3)))
hold on
xlabel('freq(Hz)')
legend('estimates','original')
title('S3')


%% section f
clc;close all;
X_fft=fft(X_org.');
X_fft1=X_fft.';
% freq 5-25 Hz
theta1=500:1:2500;
theta2=7500:1:9500;
theta=[theta1 theta2];

a=X_fft1(:,theta(1,:));
b=conj(X_fft1(:,theta(1,:)));
Sx=a*b.'/length(theta);

Cx=X_org*X_org.'/length(X_org);
[U,D]=eig(Sx,Cx);
D = diag(D);
D=real(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
S3_est=U(:,1).'*X_org;
S3_est=real(S3_est);
S_est=[S3_est; zeros(7,length(S3_est))];
X3_est=inv(U.')*S_est;
X3_est=real(X3_est);
plotEEG(X3_est,'estimated X3')
plotEEG(X3,'original X3')
RRMSE_X3 = RRMSE(X3,X3_est)

figure
subplot(1,2,1)
t=0:0.01:100-0.01;
plot(t,S3_est)
hold on
plot(t,U(:,1).'*X3)
hold on
xlabel('time(s)')
legend('estimates','original')
title('S3')
subplot(1,2,2)
f=0:0.01:100-0.01;
plot(t,abs(fft(S3_est)))
hold on
plot(t,abs(fft(U(:,1).'*X3)))
hold on
xlabel('freq(Hz)')
legend('estimates','original')
title('S3')

%% #####################################
% DSS
% section a
clc; close all;
[Z,B]=whitening(X_org);
w=rand(8,8);
for i=1:8
    
for j=1:20
    %step 1
    rp=w(:,i).'*Z;
    %step 2
    L=10000/400;
    s=zeros(1,400);
    for k=1:L
       s=s+rp(1,1+(k-1)*400:400*k); 
    end
    s=s/L;
    rp_plus=repmat(s, [1, L]);
    
    %step 3
    w_plus=Z*rp_plus.';
    %step 4
    W_orth=(eye(8)-w(:,1:i-1)*w(:,1:i-1).')*w_plus;
    W_orth=W_orth./norm(W_orth,1);
    w(:,i)=W_orth;
    
end
end
k=1:3;
s1_est=w(:,k).'*Z;

X1_est=w(:,k)*s1_est;
X1_est=inv(B)*X1_est;
plotEEG(X1_est,'estimated X1')
plotEEG(X1,'original X1')
RRMSE_X1 = RRMSE(X1,X1_est)

%% section b
clc; close all;
[Z,B]=whitening(X_org);
min_RRMSE=4;
for l=300:700
w=rand(8,8);
for i=1:8
    
for j=1:20
    %step 1
    rp=w(:,i).'*Z;
    %step 2
    L=floor(10000/l);
    s=zeros(1,l);
    for k=1:L
       s=s+rp(1,1+(k-1)*l:l*k); 
    end
    s=s/L;
    rp_plus=repmat(s, [1, L]);
    rp_pluse1=zeros(1,10000);
    rp_pluse1(1,1:length(rp_plus))=rp_plus;
    %step 3
    w_plus=Z*rp_pluse1.';
    %step 4
    W_orth=(eye(8)-w(:,1:i-1)*w(:,1:i-1).')*w_plus;
    W_orth=W_orth./norm(W_orth,1);
    w(:,i)=W_orth;
    
end
end
s1_est=w.'*Z;

X1_est=w*s1_est;
X1_est=inv(B)*X1_est;
RRMSE_X1 = RRMSE(X1,X1_est);
if RRMSE_X1<min_RRMSE
    min_RRMSE=RRMSE_X1;
    t=l;
    s1_est_min=s1_est;
    X1_est_min=X1_est;
end
end

plotEEG(X1_est_min,'estimated X1')
plotEEG(X1,'original X1')
RRMSE_X1 = RRMSE(X1,X1_est_min)
%% section c
clc; close all;
[Z,B]=whitening(X_org);
w=rand(8,8);
T=zeros(1,length(X_org));
Ton=find(T1);
T(1,Ton(1,:))=1;
for i=1:8
    
for j=1:20
    %step 1
    rp=w(:,i).'*Z;
    %step 2
    
    rp_plus=rp.*T;
    
    %step 3
    w_plus=Z*rp_plus.';
    %step 4
    W_orth=(eye(8)-w(:,1:i-1)*w(:,1:i-1).')*w_plus;
    W_orth=W_orth./norm(W_orth,1);
    w(:,i)=W_orth;
    
end
end
k=1:3;
s2_est=w(:,k).'*Z;

X2_est=w(:,k)*s2_est;
X2_est=inv(B)*X2_est;
plotEEG(X2_est,'estimated X1')
plotEEG(X2,'original X1')
RRMSE_X2 = RRMSE(X2,X2_est)
%% section d
clc; close all;
[Z,B]=whitening(X_org);
w=rand(8,8);
T=zeros(1,length(X_org));
Ton=find(T2);
T(1,Ton(1,:))=1;
for i=1:8
    
for j=1:20
    %step 1
    rp=w(:,i).'*Z;
    %step 2
    
    rp_plus=rp.*T;
    
    %step 3
    w_plus=Z*rp_plus.';
    %step 4
    W_orth=(eye(8)-w(:,1:i-1)*w(:,1:i-1).')*w_plus;
    W_orth=W_orth./norm(W_orth,1);
    w(:,i)=W_orth;
    
end
end
k=1:3;
s2_est=w(:,k).'*Z;

X2_est=w(:,k)*s2_est;
X2_est=inv(B)*X2_est;
plotEEG(X2_est,'estimated X1')
plotEEG(X2,'original X1')
RRMSE_X2 = RRMSE(X2,X2_est)
%% section e
clc; close all;
[Z,B]=whitening(X_org);
w=rand(8,8);
f=zeros(1,length(X_org));
f(1,1000:1500)=1;
f(1,8500:9000)=1;
for i=1:8
    
for j=1:20
    %step 1
    rp=w(:,i).'*Z;
    %step 2
    
    rp_plus=rp.*f;
    
    %step 3
    w_plus=Z*rp_plus.';
    %step 4
    W_orth=(eye(8)-w(:,1:i-1)*w(:,1:i-1).')*w_plus;
    W_orth=W_orth./norm(W_orth,1);
    w(:,i)=W_orth;
    
end
end
k=[3 4 6];
s3_est=w(:,k).'*Z;

X3_est=w(:,k)*s3_est;
X3_est=inv(B)*X3_est;
plotEEG(X3,'original X3')
plotEEG(X3_est,'estimated X3')
RRMSE_X3 = RRMSE(X3,X3_est)
%% section f
clc; close all;
[Z,B]=whitening(X_org);
w=rand(8,8);
f=zeros(1,length(X_org));
f(1,500:2500)=1;
f(1,7500:9500)=1;
for i=1:8
    
for j=1:20
    %step 1
    rp=w(:,i).'*Z;
    %step 2
    
    rp_plus=rp.*f;
    
    %step 3
    w_plus=Z*rp_plus.';
    %step 4
    W_orth=(eye(8)-w(:,1:i-1)*w(:,1:i-1).')*w_plus;
    W_orth=W_orth./norm(W_orth,1);
    w(:,i)=W_orth;
    
end
end
k=[3 4 6];
s3_est=w(:,k).'*Z;

X3_est=w(:,k)*s3_est;
X3_est=inv(B)*X3_est;
plotEEG(X3,'original X3')
plotEEG(X3_est,'estimated X3')
RRMSE_X3 = RRMSE(X3,X3_est)
%% Functions
function output = RRMSE(X_org,X_den)

output=sqrt(sum((X_org-X_den).^2,'all'))/sqrt(sum(X_org.^2,'all'));

end
function [Z,B]=whitening(data)
Cx = cov(data.');
[U,Diagonal] = eig(Cx);
%calculate D
d=diag(Diagonal);
Diagonal_new=diag(d.^(-0.5));
B=Diagonal_new*U.';
Z=B*data;
end
