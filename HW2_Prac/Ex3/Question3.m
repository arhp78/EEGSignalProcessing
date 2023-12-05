% Amirreza Hatamipour
% 97101507
%% Question 3
% section a
clc;clear;
data1=load('NewData1.mat');
data2=load('NewData2.mat');
data3=load('NewData3.mat');
data4=load('NewData4.mat');
EEG_Sig_1=data1.EEG_Sig;
EEG_Sig_2=data2.EEG_Sig;
EEG_Sig_3=data3.EEG_Sig;
EEG_Sig_4=data4.EEG_Sig;
plotEEG(EEG_Sig_1,'EEG Signal 1')
plotEEG(EEG_Sig_2,'EEG Signal 2')
plotEEG(EEG_Sig_3,'EEG Signal 3')
plotEEG(EEG_Sig_4,'EEG Signal 4')
%% section c
% for EEG signal 1 
Pest=32;
[F_noise1,W_noise1,K_noise1]=COM2R(EEG_Sig_1,Pest);
EEG_Sig_1_ICA=W_noise1*EEG_Sig_1;
plotEEG(EEG_Sig_1_ICA,'EEG Signal 1 ICA')

%% section d
clc
Electrodes=load('Electrodes.mat');
i=1;
% time domain
plotEEG(EEG_Sig_1_ICA,'EEG Signal 1 ICA')
figure()
for i=1:21
subplot(7,3,i)
%pwelch
pxx = pwelch(EEG_Sig_1_ICA(i,:));
pwelch(EEG_Sig_1_ICA(i,:))
title(i+" -th source")
ylabel("pow/freq")

end
figure()
for i=1:21
subplot(4,6,i)
plottopomap(Electrodes.Electrodes.X(:,1),Electrodes.Electrodes.Y(:,1),Electrodes.Electrodes.labels(1,:),F_noise1(:,i))
title(i+" -th source")
end

%%

SelSources=[2,3,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21];
EEG_Sig_1_DEN=F_noise1(:,SelSources)*EEG_Sig_1_ICA(SelSources,:);
plotEEG(EEG_Sig_1_DEN,'finall')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for second signal-EEG_sig_3
Pest=32;
[F_noise3,W_noise3,K_noise3]=COM2R(EEG_Sig_3,Pest);
EEG_Sig_3_ICA=W_noise3*EEG_Sig_3;
plotEEG(EEG_Sig_3_ICA,'EEG Signal 3 ICA')

%%
Electrodes=load('Electrodes.mat');
plotEEG(EEG_Sig_3_ICA,'EEG Signal 3 ICA')
figure()
for i=1:21
subplot(7,3,i)
%pwelch
pxx = pwelch(EEG_Sig_3_ICA(i,:));
pwelch(EEG_Sig_3_ICA(i,:))
title(i+" -th source")
ylabel("pow/freq")

end
figure()
for i=1:21
subplot(4,6,i)
plottopomap(Electrodes.Electrodes.X(:,1),Electrodes.Electrodes.Y(:,1),Electrodes.Electrodes.labels(1,:),F_noise3(:,i))
title(i+" -th source")
end
%%

SelSources1=[2,6,10,17,18,19,20,21];

EEG_Sig_3_DEN=F_noise3(:,SelSources1)*EEG_Sig_3_ICA(SelSources1,:);
plotEEG(EEG_Sig_3_DEN,'Final X Denoise')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for second signal-EEG_sig_2
Pest=32;
[F_noise2,W_noise2,K_noise2]=COM2R(EEG_Sig_2,Pest);
EEG_Sig_2_ICA=W_noise2*EEG_Sig_2;
plotEEG(EEG_Sig_2_ICA,'EEG Signal 2 ICA')

%%
Electrodes=load('Electrodes.mat');
plotEEG(EEG_Sig_2_ICA,'EEG Signal 2 ICA')
figure()
for i=1:21
subplot(7,3,i)
%pwelch
pxx = pwelch(EEG_Sig_2_ICA(i,:));
pwelch(EEG_Sig_2_ICA(i,:))
title(i+" -th source")
ylabel("pow/freq")

end
figure()
for i=1:21
subplot(4,6,i)
plottopomap(Electrodes.Electrodes.X(:,1),Electrodes.Electrodes.Y(:,1),Electrodes.Electrodes.labels(1,:),F_noise2(:,i))
title(i+" -th source")
end
%%

SelSources1=[2,3,13,17,18,19];

EEG_Sig_2_DEN=F_noise2(:,SelSources1)*EEG_Sig_2_ICA(SelSources1,:);
plotEEG(EEG_Sig_2_DEN,'Final X Denoise')
