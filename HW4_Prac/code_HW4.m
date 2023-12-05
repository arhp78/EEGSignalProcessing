% Amirreza Hatamipour
% 97101507
%% Question 1
clc;clear; close all;
EEG_ERP=load('ERP_EEG.mat');
EEG=EEG_ERP.ERP_EEG;
EEG=EEG.';
%% section a
close all;
t=0:1/240:1-1/240;
for N=100:100:2500
    EEG_mean=mean(EEG(1:N,:));
    plot(t,EEG_mean)
    xlabel('time')
    title('mean of N sample')
    grid on
    hold on
    legendInfo{N/100} = ['N= ' num2str(N)]; 
    
end
legend(legendInfo)

%% section b
clc; close all;
EEG_max=zeros(1,2550);
EEG_max(1,1)=max(abs(EEG(1,:)));
for N=2:2550
    EEG_mean=mean(EEG(1:N,:));
    EEG_max(1,N)=max(abs(EEG_mean));
  
end
figure()
plot(EEG_max)
xlabel('N')
title('max of mean of N sample')
grid on
%% section c
clc; close all;
EEG_rms=zeros(1,2550);
EEG_mean_i_1=EEG(1,:);
for N=2:2550
    EEG_mean=mean(EEG(1:N,:));
    EEG_rms_i=rms(EEG_mean-EEG_mean_i_1);
    EEG_rms(1,N)=EEG_rms_i;
    EEG_mean_i_1=EEG_mean;
end
figure()
plot(EEG_rms)
xlabel('N')
grid on
title('difference between rms')

x=500:2550;
figure()
plot(x,EEG_rms(1,500:end))
xlabel('N')
grid on
%% section e
close all;
EEG_mean=mean(EEG(1:1000,:));
EEG_mean1=mean(EEG(1:2550,:));
EEG_mean2=mean(EEG(1:333,:));
EEG_mean3= mean(EEG((randi(length(EEG),1,1000)),:));
EEG_mean4= mean(EEG((randi(length(EEG),1,333)),:));
figure()
plot(EEG_mean)
hold on
plot(EEG_mean1)
hold on
plot(EEG_mean2)
hold on
plot(EEG_mean3)
hold on
plot(EEG_mean4)
hold on
xlabel('N')
grid on
legend('N=1000','N=2550','N=N_0 /3','N=N_0 random','N=N_0/3 random')
%% Question 2
clc;clear; close all;
EEG_SSVEP=load('SSVEP_EEG.mat');
EEG=EEG_SSVEP.SSVEP_Signal;
Events=EEG_SSVEP.Events;
Event_samples=EEG_SSVEP.Event_samples;
Fs=250;
%% section a1
EEG_filtered = zeros(6,117917);
for i=1:6
    EEG_filtered(i,:) = bandpass( EEG(i,:) , [ 1 40] , Fs);
end
%% section a2
clc;
EEG_ss=zeros(15,2);
for i=1:15
    EEG_ss(i,1)=Event_samples(1,i);
     EEG_ss(i,2)=250*5+Event_samples(1,i);
end
%% section a3
clc; close all;
figure()
for i=1:15
subplot(3,5,i)
%pwelch
for j=1:6
pxx = pwelch(EEG_filtered(j,EEG_ss(i,1):EEG_ss(i,2)));
%pwelch(EEG_filtered(j,EEG_ss(i,1):EEG_ss(i,2)))
plot(pow2db(pxx))
%plot(pxx)
hold on
end
title(i+" -th trial")
ylabel("pow/freq")
xlabel("frequency")
grid on
legend(subplot(3,5,1),{'Pz','Qz','P7','P8','O2','O1'})

end
%% section 2b
clc;
f=zeros(1,5);
y=zeros(1,5);
for i=1:5
   f(1,i)= Events(1,3*(i-1)+1);
end
t=5*Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc=length(m)/15
%% section 2c
clc;
f=zeros(1,5);
y=zeros(1,5);
for i=1:5
   f(1,i)= Events(1,3*(i-1)+1);
end
t=5*Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_Pz=length(m)/15

for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_P7=length(m)/15

for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_O2=length(m)/15
%% section 2d
clc;
f=zeros(1,5);
y=zeros(1,5);
for i=1:5
   f(1,i)= Events(1,3*(i-1)+1);
end
t=3*Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_3T=length(m)/15

t=2*Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_2T=length(m)/15

t=Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_T=length(m)/15

%% Question 3
clc;clear; close all;
EEG=load('Ex3.mat');
TestData=EEG.TestData;
TrainData=EEG.TrainData;
TrainLabel=EEG.TrainLabel;

%% CSP
clc; close all;
n=find(TrainLabel);
Cx1=0;
for i=1:n
Cx1=Cx1+TrainData(:,:,n(i))*TrainData(:,:,n(i)).'/trace(TrainData(:,:,n(i))*TrainData(:,:,n(i)).');
end
Cx1=Cx1/length(n);

m=find(TrainLabel==0);
Cx2=0;
for i=1:m
Cx2=Cx2+TrainData(:,:,m(i))*TrainData(:,:,m(i)).'./trace(TrainData(:,:,m(i))*TrainData(:,:,m(i)).');
end
Cx2=Cx2/length(m);

% GEVD
[U,D]=eig(Cx1,Cx2);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
W=[U(:,1) U(:,end)];
class0=W.'*TrainData(:,:,m(1,1));

class1=W.'*TrainData(:,:,n(1,1));
figure()
subplot(2,1,1)
plot(class0(1,:))
hold on
plot(class1(1,:))
title('first filter')
legend('class0','class1')
grid on

subplot(2,1,2)
plot(class0(2,:))
hold on
plot(class1(2,:))
title('last filter')
legend('class0','class1')
grid on
%% section b
Electrodes=load('AllElectrodes.mat');
clc; close all;
elecX=zeros(30,1);
elecY=zeros(30,1);
label=[];

i=[37 7 5 38 40 42 10 47 45  15 13 48 50 52 18 32 55 23 22 21 20 31 57  58 59 60 26 25 27 64];
for j=1:30
    elecX(j)=Electrodes.AllElectrodes(i(j)).X;
    elecY(j)=Electrodes.AllElectrodes(i(j)).Y;
    label{j}=Electrodes.AllElectrodes(i(j)).labels;
    
end
subplot(1,2,1)
plottopomap(elecX,elecY, label,abs(U(:,end)))
title('last filter')
subplot(1,2,2)
plottopomap(elecX,elecY, label,abs(U(:,1)))
title('first filter')
%% section 3 c
clc;
a=[1 56 110 165];

for filter_num=1:15
    acc=0;
    for i=1:3
    test_3fold=TrainData(:,:,a(i):a(i+1));
    test_3fold_label=TrainLabel(1,a(i):a(i+1));
    train_3fold=TrainData;
    train_3fold(:,:,a(i):a(i+1))=[];
    train_3fold_label=TrainLabel;
    train_3fold_label(:,a(i):a(i+1))=[];
    n=find(train_3fold_label);
    Cx1=0;
    for j=1:n
        Cx1=Cx1+train_3fold(:,:,n(j))*train_3fold(:,:,n(j)).'/trace(train_3fold(:,:,n(j))*train_3fold(:,:,n(j)).');
    end
    Cx1=Cx1/length(n);

    m=find(train_3fold_label==0);
    Cx2=0;
    for j=1:m
        Cx2=Cx2+train_3fold(:,:,m(j))*train_3fold(:,:,m(j)).'./trace(train_3fold(:,:,m(j))*train_3fold(:,:,m(j)).');
    end
    Cx2=Cx2/length(m);

    % GEVD
    [U,D]=eig(Cx1,Cx2);
    D = diag(D);
    [D, index] = sort(D, 'descend');
    U = (U(:,index));
    W=[U(:,1:filter_num) U(:,end-filter_num+1:end)];
    
    var_class0=zeros(length(m),2*filter_num);
    for j=1: length(m)
        class0=W.'*train_3fold(:,:,m(1,j));
        var_class0(j,:)=var(class0');
    end
    
    var_class1=zeros(length(n),2*filter_num);
    for j=1: length(n)
        class1=W.'*train_3fold(:,:,n(1,j));
        var_class1(j,:)=var(class1');
    end
    
    var_test=zeros(length(test_3fold(1,1,:)),2*filter_num);
    for j=1: length(test_3fold(1,1,:))
        class_test=W.'*test_3fold(:,:,j);
        var_test(j,:)=var(class_test');
    end
    
    feature=[var_class0 ; var_class1];
    label=[zeros(length(m),1); ones(length(n),1)];
    Mdl = fitcknn(feature,label);
    label_test_predict = predict(Mdl,var_test);
    

    acc0=length(find(label_test_predict'==test_3fold_label))/length(test_3fold_label);
    acc=acc+acc0;
    end
  accuracy(filter_num)=acc/3;
end
x=1:15;
plot(2*x,accuracy)
xlabel('number of filter')
ylabel('acc')
grid on
%% section d
%best filter num
filter_num=5;
acc=0;
TestData=EEG.TestData;
TrainData=EEG.TrainData;
TrainLabel=EEG.TrainLabel;
    n=find(TrainLabel);
    Cx1=0;
    for j=1:n
        Cx1=Cx1+TrainData(:,:,n(j))*TrainData(:,:,n(j)).'/trace(TrainData(:,:,n(j))*TrainData(:,:,n(j)).');
    end
    Cx1=Cx1/length(n);

    m=find(TrainLabel==0);
    Cx2=0;
    for j=1:m
        Cx2=Cx2+TrainData(:,:,m(j))*TrainData(:,:,m(j)).'./trace(TrainData(:,:,m(j))*TrainData(:,:,m(j)).');
    end
    Cx2=Cx2/length(m);

    % GEVD
    [U,D]=eig(Cx1,Cx2);
    D = diag(D);
    [D, index] = sort(D, 'descend');
    U = (U(:,index));
    W=[U(:,1:filter_num) U(:,end-filter_num+1:end)];
    
    var_class0=zeros(length(m),2*filter_num);
    for j=1: length(m)
        class0=W.'*TrainData(:,:,m(1,j));
        var_class0(j,:)=var(class0');
    end
    
    var_class1=zeros(length(n),2*filter_num);
    for j=1: length(n)
        class1=W.'*TrainData(:,:,n(1,j));
        var_class1(j,:)=var(class1');
    end
    
    var_test=zeros(length(TestData(1,1,:)),2*filter_num);
    for j=1: length(TestData(1,1,:))
        class_test=W.'*TestData(:,:,j);
        var_test(j,:)=var(class_test');
    end
    
    feature=[var_class0 ; var_class1];
    label=[zeros(length(m),1); ones(length(n),1)];
    Mdl = fitcknn(feature,label);
    label_test_predict = predict(Mdl,var_test);
 %% save   
save('test_label.mat',label_test_predict')


%% Function
function y=generator(f,t)
    num=floor(40/f);
    y=zeros(2*num,t);
    n=t/250;
    T=[0:1/250:n-1/250];
    for i=1:num
       y(2*i-1,:)= sin(2*pi*i*f*T);
       y(2*i,:)= cos(2*pi*i*f*T);
    end
end

function plottopomap(elocsX,elocsY,elabels,data)

% define XY points for interpolation
interp_detail = 100;
interpX = linspace(min(elocsX)-.2,max(elocsX)+.25,interp_detail);
interpY = linspace(min(elocsY),max(elocsY),interp_detail);

% meshgrid is a function that creates 2D grid locations based on 1D inputs
[gridX,gridY] = meshgrid(interpX,interpY);
% Interpolate the data on a 2D grid
interpFunction = TriScatteredInterp(elocsY,elocsX,data);
topodata = interpFunction(gridX,gridY);

% plot map
contourf(interpY,interpX,topodata);
hold on
scatter(elocsY,elocsX,10,'ro','filled');
for i=1:length(elocsX)
    text(elocsY(i),elocsX(i),elabels(i))
end
set(gca,'xtick',[])
set(gca,'ytick',[])
end