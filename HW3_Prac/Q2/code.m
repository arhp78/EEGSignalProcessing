%Amirreza Hatamipour
% 97101507
%% read data
clc;clear; close all;
noise= load("contaminated.mat");
signal= load("pure.mat");
pure=signal.pure;
X=noise.contaminated;
Fs=200;
% plot
disp_eeg(pure, '',Fs,'','signal');
disp_eeg(X, '',Fs,'','signal+noise');
disp_eeg(X-pure, '',Fs,'','noise');


RRMSE(pure,X)
Ton=[1.56*Fs:1:3.06*Fs 16.95*Fs:1:16.42*Fs 19.955*Fs:1:21.965*Fs 22.06*Fs:1:22.73*Fs];
%% GEVD
clc;close all;
Ton=int64(Ton);
T=zeros(1,length(X));
T(1,Ton(1,:))=1;
%X_Ton=X.*T;
Px=X(:,Ton(1,:))*X(:,Ton(1,:)).'/length(Ton);

Cx=X*X.'/length(X);
[U,D]=eig(Px,Cx);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
sEOG_est1=U(:,1).'*X;
sEOG_est2=U(:,2).'*X;

S_est=[sEOG_est1; sEOG_est2;zeros(17,length(sEOG_est1))];
EOG_est=inv(U.')*S_est;
RRMSE_X2 = RRMSE(X-pure,EOG_est)
disp_eeg(EOG_est, '',Fs,'','estimated EOG');
disp_eeg(X-pure, '',Fs,'','original EOG');

X_den=X-EOG_est;
disp_eeg(X_den, '',Fs,'','estimated denoise');
disp_eeg(pure, '',Fs,'','original EEG without noise');
RRMSE(pure,X_den)
%% DSS
clc; close all;
[Z,B]=whitening(X);
a=cov(Z.');
w=rand(19,10);
T=zeros(1,length(X));
T(1,Ton(1,:))=1;
for i=1:10
    
for j=1:20
    %step 1
    rp=w(:,i).'*Z;
    %step 2
    
    rp_plus=rp.*T;
    
    %step 3
    w_plus=Z*rp_plus.';
    %step 4
    W_orth=(eye(19)-w(:,1:i-1)*w(:,1:i-1).')*w_plus;
    W_orth=W_orth./norm(W_orth,1);
    w(:,i)=W_orth;
    
end
end
sEOG_est=w(:,1:8).'*Z;
EOG_est=w(:,1:8)*sEOG_est;
EOG_est=inv(B)*EOG_est;
RRMSE_X2 = RRMSE(X-pure,EOG_est)
disp_eeg(EOG_est, '',Fs,'','estimated EOG');
disp_eeg(X-pure, '',Fs,'','original EOG');

X_den=X-EOG_est;
disp_eeg(X_den, '',Fs,'','estimated denoise');
disp_eeg(pure, '',Fs,'','original EEG without noise');
RRMSE(pure,X_den)
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
