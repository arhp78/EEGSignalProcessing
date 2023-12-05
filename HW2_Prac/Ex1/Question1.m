% Amirreza Hatamipour
% 97101507
%% Question 1
%% section a
clc;clear;
data=load('Ex1.mat');
EEG_data=data.X;

%csvwrite('file.csv',EEG_data)

% plot
figure()
subplot(2,2,1)
scatter3(EEG_data(1,:),EEG_data(2,:),EEG_data(3,:),'r','.');
xlabel('x')
ylabel('y')
zlabel('z')
subplot(2,2,2)
scatter(EEG_data(2,:),EEG_data(3,:),'r','.');
xlabel('y')
ylabel('z')
grid on
subplot(2,2,3)
scatter(EEG_data(1,:),EEG_data(3,:),'r','.');
xlabel('x')
ylabel('z')
grid on
subplot(2,2,4)
scatter(EEG_data(1,:),EEG_data(1,:),'r','.');
xlabel('x')
ylabel('y')
grid on

%% section b
% PCA
clc;
Cx = cov(EEG_data.');
[U,Diagonal] = eig(Cx);
%calculate D
d=diag(Diagonal);
Diagonal_new=diag(d.^(-0.5));
D=Diagonal_new*U.';
y=D*EEG_data;

Cy=D*Cx*D'
y_matalab=D.'*EEG_data;
%plot

scatter3(y(1,:),y(2,:),y(3,:),'r','.');
xlabel('x')
ylabel('y')
zlabel('z')
title('Whitening  data')




figure
subplot(2,2,1)
biplot(D(:,1:3),'varlabels',{'v_1','v_2','v_3'})
title('PCA component of data')

subplot(2,2,2)
biplot(D(:,1:3),'varlabels',{'v_1','v_2','v_3'})
hold on
scatter3(EEG_data(1,:),EEG_data(2,:),EEG_data(3,:),'g','.');
title('data on component axis')

subplot(2,2,3)
scatter3(EEG_data(1,:),EEG_data(2,:),EEG_data(3,:),'g','.');
xlabel('x')
ylabel('y')
zlabel('z')
hold on
biplot(D(:,1:3),'varlabels',{'v_1','v_2','v_3'})
title('data with PCA coponents')

subplot(2,2,4)
scatter3(EEG_data(1,:),EEG_data(2,:),EEG_data(3,:),'g','.');
xlabel('x')
ylabel('y')
zlabel('z')
title('Original Data')

%% section c
% matlab PCA
[coeff,score,latent] = pca(EEG_data.');
D_matlab=diag(latent.^(-0.5))*coeff.';

y_matalb=D_matlab*EEG_data;
Cy_matlab1 = cov(y_matalb.')
Cy_matlab=D_matlab*Cx*D_matlab.'
%Cy(:,1)=Cy(:,1)./sum(coeff(:,1).^2);
%plot
scatter3(y_matalb(1,:),y_matalb(2,:),y_matalb(3,:),'r','.');
xlabel('x')
ylabel('y')
zlabel('z')
title('Whitening  data')

figure
subplot(2,2,1)
biplot(D_matlab(:,1:3),'varlabels',{'v_1','v_2','v_3'})
title('PCA component of data')

subplot(2,2,2)
biplot(D_matlab(:,1:3),'varlabels',{'v_1','v_2','v_3'})
hold on
scatter3(EEG_data(1,:),EEG_data(2,:),EEG_data(3,:),'g','.');
title('data on component axis')

subplot(2,2,3)
scatter3(EEG_data(1,:),EEG_data(2,:),EEG_data(3,:),'g','.');
xlabel('x')
ylabel('y')
zlabel('z')
hold on
biplot(D_matlab(:,1:3),'varlabels',{'v_1','v_2','v_3'})
title('data with PCA coponents')

subplot(2,2,4)
scatter3(EEG_data(1,:),EEG_data(2,:),EEG_data(3,:),'g','.');
xlabel('x')
ylabel('y')
zlabel('z')
title('Original Data')

%% section d
[U,S,V] = svd(EEG_data.','econ');
figure()

scatter3(U(:,1) , U(:,2) , U(:,3),'b')
title("Whitening data")

Cy_SVD = cov(U)
figure()
scatter3(EEG_data(1,:) , EEG_data(2,:) ,EEG_data(3,:),'r','.')
hold on
biplot(V(:,1:3),'varlabels',{'v_1','v_2','v_3'})



