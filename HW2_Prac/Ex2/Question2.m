% Amirreza Hatamipour
% 97101507
%% Question 2
% section a
clc;clear;
data=load('Ex2.mat');

% plot EEG signal
X=data.X_org;  % 32-channel data
plotEEG(X,'original')
%% plot noise_1
X=data.X_noise_1;  % 32-channel noise
plotEEG(X,'noise1')
%% plot noise_2
X=data.X_noise_2;  % 32-channel noise
plotEEG(X,'noise2')
%% section a 
% noise 1, SNR=10db

X_org=data.X_org;
noise_1=data.X_noise_1;
SNR=-10;
%generate new signal with noise
X_noise1_10db=add_noise(X_org,noise_1,SNR);
%check
snr_10db=snr(X_org,X_noise1_10db-X_org)



%plot
plotEEG(X_noise1_10db,'X with noise -10db')


%% section b
% Source separation with PCA for noise_1 & SNR=-10 
[coeff,score,latent]  = pca(X_noise1_10db.');
D_noise1_10db=diag(latent.^(-0.5))*coeff.';
Y_noise1_10db_PCA =D_noise1_10db*X_noise1_10db;
plotEEG(Y_noise1_10db_PCA,'source separation with PCA for noise1 & -10db ')
%%
% Source separation with ICA for noise_1 & SNR=-10 
%close all
Pest=32;
[F,W,K]=COM2R(X_noise1_10db,Pest);
Y_noise1_10db_ICA=W*X_noise1_10db;
Cy=cov(Y_noise1_10db_ICA');
plotEEG(Y_noise1_10db_ICA,'source separation with ICA for noise1 & -10db ')
%% section c - noise reduction
clc;
Y_noise1_10db_PCA(24:32,:)=0;
Y_noise1_10db_PCA(1:21,:)=0;
plotEEG(Y_noise1_10db_PCA,'source separation with PCA for noise1 & -10db after zero some channel')


Y_noise1_10db_ICA(4:end,:)=0;
Y_noise1_10db_ICA(1:2,:)=0;

plotEEG(Y_noise1_10db_ICA,'source separation with ICA for noise1 & -10db after zero some channel')
%% section d

X_noise1_10db_DEN_PCA=D_noise1_10db'*Y_noise1_10db_PCA;
plotEEG(X_noise1_10db_DEN_PCA,'final source with -10db PCA')

X_noise1_10db_DEN_ICA=F*Y_noise1_10db_ICA;
plotEEG(X_noise1_10db_DEN_ICA,'final source with -10db ICA')

%% section e

%channel 13
figure()
suptitle('channel 13') 
subplot(2,2,1)
plot(X_noise1_10db_DEN_PCA(13,:))
title('denoise with PCA')
subplot(2,2,2)
plot(X_noise1_10db_DEN_ICA(13,:))
title('denoise with ICA')
subplot(2,2,3)
plot(X_noise1_10db(13,:))
title('noisy signal')
subplot(2,2,4)
plot(X_org(13,:))
title('original signal')
% channel 24
figure()
suptitle('channel 24') 
subplot(2,2,1)
plot(X_noise1_10db_DEN_PCA(24,:))
title('denoise with PCA')
subplot(2,2,2)
plot(X_noise1_10db_DEN_ICA(24,:))
title('denoise with ICA')
subplot(2,2,3)
plot(X_noise1_10db(24,:))
title('noisy signal')
subplot(2,2,4)
plot(X_org(24,:))
title('original signal')
%plot(X_noise1_10db_DEN_PCA(13,:))
%plotEEG([X_noise1_10db_DEN_PCA([13,24],:);X_org([13,24],:)])
%plotEEG(X_noise1_10db_DEN_PCA([13,24],:))
%plotEEG(X_org([13,24],:))

%% section f - RRMSE
RRMSE_noise1_10db_PCA=RRMSE(X_org,X_noise1_10db_DEN_PCA)
RRMSE_noise1_10db_ICA=RRMSE(X_org,X_noise1_10db_DEN_ICA)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% section a 
% noise 1, SNR=-20db

X_org=data.X_org;
noise_1=data.X_noise_1;
SNR=-20;
%generate new signal with noise
X_noise1_20db=add_noise(X_org,noise_1,SNR);
%check
snr_20db=snr(X_org,X_noise1_20db-X_org)



%plot
plotEEG(X_noise1_20db,'X with noise -20db')


%% section b
% Source separation with PCA for noise_1 & SNR=-10 
[coeff,score,latent]  = pca(X_noise1_20db.');
D_noise1_20db=diag(latent.^(-0.5))*coeff.';
Y_noise1_20db_PCA =D_noise1_20db*X_noise1_20db;
plotEEG(Y_noise1_20db_PCA,'source separation with PCA for noise1 & -20db ')
%%
% Source separation with ICA for noise_1 & SNR=-10 
%close all
Pest=32;
[F,W,K]=COM2R(X_noise1_20db,Pest);
Y_noise1_20db_ICA=W*X_noise1_20db;
Cy=cov(Y_noise1_20db_ICA');
plotEEG(Y_noise1_20db_ICA,'source separation with ICA for noise1 & -20db ')
%% section c - noise reduction
clc;
Y_noise1_20db_PCA(24:32,:)=0;
Y_noise1_20db_PCA(1:21,:)=0;
plotEEG(Y_noise1_20db_PCA,'source separation with PCA for noise1 & -20db after zero some channel')


Y_noise1_20db_ICA(17:end,:)=0;
Y_noise1_20db_ICA(1:15,:)=0;

plotEEG(Y_noise1_20db_ICA,'source separation with ICA for noise1 & -20db after zero some channel')
%% section d

X_noise1_20db_DEN_PCA=D_noise1_20db'*Y_noise1_20db_PCA;
plotEEG(X_noise1_20db_DEN_PCA,'final source with -20db PCA')

X_noise1_20db_DEN_ICA=F*Y_noise1_20db_ICA;
plotEEG(X_noise1_20db_DEN_ICA,'final source with -20db ICA')

%% section e

%channel 13
figure()
suptitle('channel 13') 
subplot(2,2,1)
plot(X_noise1_20db_DEN_PCA(13,:))
title('denoise with PCA')
subplot(2,2,2)
plot(X_noise1_20db_DEN_ICA(13,:))
title('denoise with ICA')
subplot(2,2,3)
plot(X_noise1_20db(13,:))
title('noisy signal')
subplot(2,2,4)
plot(X_org(13,:))
title('original signal')
% channel 24
figure()
suptitle('channel 24') 
subplot(2,2,1)
plot(X_noise1_20db_DEN_PCA(24,:))
title('denoise with PCA')
subplot(2,2,2)
plot(X_noise1_20db_DEN_ICA(24,:))
title('denoise with ICA')
subplot(2,2,3)
plot(X_noise1_20db(24,:))
title('noisy signal')
subplot(2,2,4)
plot(X_org(24,:))
title('original signal')

%% section f - RRMSE
RRMSE_noise1_20db_PCA=RRMSE(X_org,X_noise1_20db_DEN_PCA)
RRMSE_noise1_20db_ICA=RRMSE(X_org,X_noise1_20db_DEN_ICA)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% section a 
% noise 2, SNR=-10db

X_org=data.X_org;
noise_2=data.X_noise_2;
SNR=-10;
%generate new signal with noise
X_noise2_10db=add_noise(X_org,noise_2,SNR);
%check
snr_10db=snr(X_org,X_noise2_10db-X_org)



%plot
plotEEG(X_noise2_10db,'X with noise2 -10db')


%% section b
% Source separation with PCA for noise_2 & SNR=-10 
[coeff,score,latent]  = pca(X_noise2_10db.');
D_noise2_10db=diag(latent.^(-0.5))*coeff.';
Y_noise2_10db_PCA =D_noise2_10db*X_noise2_10db;
plotEEG(Y_noise2_10db_PCA,'source separation with PCA for noise2 & -10db ')
%% Source separation with ICA for noise_2 & SNR=-10 
%close all
Pest=32;
[F,W,K]=COM2R(X_noise2_10db,Pest);
Y_noise2_10db_ICA=W*X_noise2_10db;
Cy=cov(Y_noise2_10db_ICA');
plotEEG(Y_noise2_10db_ICA,'source separation with ICA for noise2 & -10db ')
%% section c - noise reduction
clc;
Y_noise2_10db_PCA(5:32,:)=0;
Y_noise2_10db_PCA(1,:)=0;
plotEEG(Y_noise2_10db_PCA,'source separation with PCA for noise2 & -10db after zero some channel')

Y_noise2_10db_ICA(13:end,:)=0;
Y_noise2_10db_ICA(7:11,:)=0;
Y_noise2_10db_ICA(1:5,:)=0;

plotEEG(Y_noise2_10db_ICA,'source separation with ICA for noise2 & -10db after zero some channel')
%% section d

X_noise2_10db_DEN_PCA=D_noise2_10db'*Y_noise2_10db_PCA;
plotEEG(X_noise2_10db_DEN_PCA,'final source with -10db PCA')

X_noise2_10db_DEN_ICA=F*Y_noise2_10db_ICA;
plotEEG(X_noise2_10db_DEN_ICA,'final source with -10db ICA')


%% section e

%channel 13
figure()
suptitle('channel 13') 
subplot(2,2,1)
plot(X_noise2_10db_DEN_PCA(13,:))
title('denoise with PCA')
subplot(2,2,2)
plot(X_noise2_10db_DEN_ICA(13,:))
title('denoise with ICA')
subplot(2,2,3)
plot(X_noise2_10db(13,:))
title('noisy signal')
subplot(2,2,4)
plot(X_org(13,:))
title('original signal')
% channel 24
figure()
suptitle('channel 24') 
subplot(2,2,1)
plot(X_noise2_10db_DEN_PCA(24,:))
title('denoise with PCA')
subplot(2,2,2)
plot(X_noise2_10db_DEN_ICA(24,:))
title('denoise with ICA')
subplot(2,2,3)
plot(X_noise2_10db(24,:))
title('noisy signal')
subplot(2,2,4)
plot(X_org(24,:))
title('original signal')

%% section f - RRMSE
RRMSE_noise2_10db_PCA=RRMSE(X_org,X_noise2_10db_DEN_PCA)
RRMSE_noise2_10db_ICA=RRMSE(X_org,X_noise2_10db_DEN_ICA)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% section a 
% noise 2, SNR=-20db

X_org=data.X_org;
noise_2=data.X_noise_2;
SNR=-20;
%generate new signal with noise
X_noise2_20db=add_noise(X_org,noise_2,SNR);
%check
snr_10db=snr(X_org,X_noise2_20db-X_org)



%plot
plotEEG(X_noise2_20db,'X with noise2 -20db')


%% section b
% Source separation with PCA for noise_1 & SNR=-10 
[coeff,score,latent]  = pca(X_noise2_20db.');
D_noise2_20db=diag(latent.^(-0.5))*coeff.';
Y_noise2_20db_PCA =D_noise2_20db*X_noise2_20db;
plotEEG(Y_noise2_20db_PCA,'source separation with PCA for noise2 & -20db ')
%%
% Source separation with ICA for noise_1 & SNR=-10 
%close all
Pest=32;
[F,W,K]=COM2R(X_noise2_20db,Pest);
Y_noise2_20db_ICA=W*X_noise2_20db;
Cy=cov(Y_noise2_20db_ICA');
plotEEG(Y_noise2_20db_ICA,'source separation with ICA for noise2 & -20db ')
%% section c - noise reduction
clc;
Y_noise2_20db_PCA(17:32,:)=0;
Y_noise2_20db_PCA(1:14,:)=0;
plotEEG(Y_noise2_20db_PCA,'source separation with PCA for noise2 & -20db after zero some channel')


Y_noise2_20db_ICA(16:end,:)=0;
Y_noise2_20db_ICA(1:14,:)=0;

plotEEG(Y_noise2_20db_ICA,'source separation with ICA for noise2 & -20db after zero some channel')
%% section d

X_noise2_20db_DEN_PCA=D_noise2_20db'*Y_noise2_20db_PCA;
plotEEG(X_noise2_20db_DEN_PCA,'final source with -20db PCA')


X_noise2_20db_DEN_ICA=F*Y_noise2_20db_ICA;
plotEEG(X_noise2_20db_DEN_ICA,'final source with -20db ICA')


%% section e

%channel 13
figure()
suptitle('channel 13') 
subplot(2,2,1)
plot(X_noise2_20db_DEN_PCA(13,:))
title('denoise with PCA')
subplot(2,2,2)
plot(X_noise2_20db_DEN_ICA(13,:))
title('denoise with ICA')
subplot(2,2,3)
plot(X_noise2_20db(13,:))
title('noisy signal')
subplot(2,2,4)
plot(X_org(13,:))
title('original signal')
% channel 24
figure()
suptitle('channel 24') 
subplot(2,2,1)
plot(X_noise2_20db_DEN_PCA(24,:))
title('denoise with PCA')
subplot(2,2,2)
plot(X_noise2_20db_DEN_ICA(24,:))
title('denoise with ICA')
subplot(2,2,3)
plot(X_noise2_20db(24,:))
title('noisy signal')
subplot(2,2,4)
plot(X_org(24,:))
title('original signal')

%% section f - RRMSE
RRMSE_noise2_20db_PCA=RRMSE(X_org,X_noise2_20db_DEN_PCA)
RRMSE_noise2_20db_ICA=RRMSE(X_org,X_noise2_20db_DEN_ICA)
%% Functions
function output = RRMSE(X_org,X_den)

output=sqrt(sum((X_org-X_den).^2,'all'))/sqrt(sum(X_org.^2,'all'));

end
function x_noise=add_noise(X_org,noise,SNR)
% first calculate P_S and P_N
P_S=sum(X_org.^2,'all');
P_N=sum(noise.^2,'all');
sigma_noise=(P_S/P_N)*10^(-SNR/10);

%generate new signal with noise
x_noise=X_org+sqrt(sigma_noise)*noise;
end