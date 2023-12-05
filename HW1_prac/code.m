%Amirreza Hatamipour
%97101507
%% Qustion 1
%section a
%x(t)=cos(2pif(t)t)
%f(t)=100+100t^2
clc;
t=0:1/1e3:2;
fo=100;
f1=200;
y = chirp(t,fo,1,f1,'quadratic');

figure;
plot(t,y)
title('signal in time domain x(t)=cos(2*pi*(100+100t^2)*t)')
xlabel('time(s)')
%% section b
wvtool(rectwin(128))
wvtool(triang(128))
wvtool(gausswin(128))
wvtool(hamming(128))
wvtool(rectwin(128),triang(128),gausswin(128),hamming(128))

%% section c
clc;
fs=1000;
nfft=128;overlap=0;
figure;
subplot(2,2,1)
spectrogram(y,rectwin(128),overlap,nfft,fs,'yaxis');
title('STFT of Rectangle')
subplot(2,2,2)
spectrogram(y,triang(128),overlap,nfft,fs,'yaxis');
title('STFT of Tringle')
subplot(2,2,3)
spectrogram(y,gausswin(128),overlap,nfft,fs,'yaxis');
title('STFT of Gauss')
subplot(2,2,4)
spectrogram(y,hamming(128),overlap,nfft,fs,'yaxis');
title('STFT of Hamming')

%% section d
overlap1=0;
overlap2=64;
overlap3=127;
figure;
subplot(1,3,1)
spectrogram(y,triang(128),overlap1,nfft,fs,'yaxis');
title('STFT of Tringle, overlap=0')
subplot(1,3,2)
spectrogram(y,triang(128),overlap2,nfft,fs,'yaxis');
title('STFT of Tringle, overlap=64')
subplot(1,3,3)
spectrogram(y,triang(128),overlap3,nfft,fs,'yaxis');
title('STFT of Tringle, overlap=127')
%% section e
L1=32;
L2=128;
L3=512;
figure;
subplot(1,3,1)
spectrogram(y,triang(L1),L1-1,L1,fs,'yaxis');
title('STFT of Tringle, length=32, overlap=31')
subplot(1,3,2)
spectrogram(y,triang(L2),L2-1,L2,fs,'yaxis');
title('STFT of Tringle, length=128, overlap=127')
subplot(1,3,3)
spectrogram(y,triang(L3),L3-1,L3,fs,'yaxis');
title('STFT of Tringle, length=512, overlap=511')
%% section f
L=128; overlap=L/2;
figure;
subplot(1,3,1)
spectrogram(y,triang(L),overlap,L,fs,'yaxis');
title('STFT of Tringle,nfft=L')
subplot(1,3,2)
spectrogram(y,triang(L),overlap,2*L,fs,'yaxis');
title('STFT of Tringle,nfft=2L')
subplot(1,3,3)
spectrogram(y,triang(L),overlap,4*L,fs,'yaxis');
title('STFT of Tringle,nfft=4L')
%% section g

L=128;nfft=L;
overlap=64;
spectrogram_fft(y,L,overlap,nfft)

%% Question 2
clc;clear;
Fs=256;
t=0:1/Fs:2;
data=load("NewEEGSignal.mat");
mat=data.NewEEGSignal;


subplot(1,3,1)
plot(t(1,1:end-1),mat)
title('EEG signal in time domain')
xlabel('tmie(s)')
ylabel('amplitude')

fft_EEG=fftshift(fft(mat));
w= linspace(-Fs/2,Fs/2,512);
subplot(1,3,2)
plot(w,abs(fft_EEG))
title('EEG signal in ferquency domain')
xlabel('F(Hz)')
ylabel('abs(FFT)')

subplot(1,3,3)
spectrogram(mat,gausswin(Fs),Fs/2,Fs,Fs,'yaxis');
title('EEG signal in time domain')

% after edition
figure
subplot(1,3,1)
plot(t(1,1:end-1),mat)
title('EEG signal in time domain')
xlabel('tmie(s)')
ylabel('amplitude')

fft_EEG=fftshift(fft(mat));
w= linspace(-Fs/2,Fs/2,512);
subplot(1,3,2)
plot(w,abs(fft_EEG))
title('EEG signal in ferquency domain')
xlabel('F(Hz)')
ylabel('abs(FFT)')
xlim([-64 64])
subplot(1,3,3)
spectrogram(mat,gausswin(Fs),Fs/2,Fs,Fs,'yaxis');
title('EEG signal in time domain')
ylim([0 64])
%% section b
Fs=256;
t=0:1/Fs:2;
fpass=64;
mat_LP = lowpass(mat,fpass,Fs);
figure
subplot(1,3,1)
plot(t(1,1:end-1),mat_LP)
title('EEG signal in time domain')
xlabel('tmie(s)')
ylabel('amplitude')

fft_EEG=fftshift(fft(mat_LP));
w= linspace(-Fs/2,Fs/2,512);
subplot(1,3,2)
plot(w,abs(fft_EEG))
title('EEG signal in ferquency domain')
xlabel('F(Hz)')
ylabel('abs(FFT)')

subplot(1,3,3)
spectrogram(mat_LP,gausswin(Fs),Fs/2,Fs,Fs,'yaxis');
title('EEG signal in time domain')

% now downsample signal
Fs=128;
t=0:1/Fs:2;

EEG_ds = downsample(mat,2);
figure
subplot(1,3,1)
plot(t(1,1:end-1),EEG_ds)
title('EEG signal in time domain')
xlabel('tmie(s)')
ylabel('amplitude')

fft_EEG_ds=fftshift(fft(EEG_ds));
w= linspace(-Fs/2,Fs/2,256);
subplot(1,3,2)
plot(w,abs(fft_EEG_ds))
title('EEG signal in ferquency domain')
xlabel('F(Hz)')
ylabel('abs(FFT)')

subplot(1,3,3)
spectrogram(mat_LP,gausswin(Fs),Fs/2,Fs,Fs,'yaxis');
title('EEG signal in time domain')
%% section c
N=256;
w3= linspace(-Fs/2,Fs/2,32);
w2= linspace(-Fs/2,Fs/2,64);
w1= linspace(-Fs/2,Fs/2,128);

fft_EEG_ds_N8=fftshift(fft(EEG_ds,N/8));
fft_EEG_ds_N4=fftshift(fft(EEG_ds,N/4));
fft_EEG_ds_N2=fftshift(fft(EEG_ds,N/2));

subplot(1,3,1)
plot(w3,abs(fft_EEG_ds_N8))
title('FFT N/8 points')
xlabel('F(Hz)')
ylabel('abs(FFT)')


subplot(1,3,2)
plot(w2,abs(fft_EEG_ds_N4))
title('FFT N/4 points')
xlabel('F(Hz)')
ylabel('abs(FFT)')

subplot(1,3,3)
plot(w1,abs(fft_EEG_ds_N2))
xlabel('F(Hz)')
ylabel('abs(FFT)')
title('FFT N/2 points')
%% section d

EEG_ds3_1=[EEG_ds(1:N/8) zeros(1,224)];
EEG_ds3_2=[EEG_ds(N/8+1:2*N/8) zeros(1,224)];
EEG_ds3_3=[EEG_ds(2*N/8+1:3*N/8) zeros(1,224)];
EEG_ds3_4=[EEG_ds(3*N/8+1:4*N/8) zeros(1,224)];
EEG_ds3_5=[EEG_ds(4*N/8+1:5*N/8) zeros(1,224)];
EEG_ds3_6=[EEG_ds(5*N/8+1:6*N/8) zeros(1,224)];
EEG_ds3_7=[EEG_ds(6*N/8+1:7*N/8) zeros(1,224)];
EEG_ds3_8=[EEG_ds(7*N/8+1:N) zeros(1,224)];


subplot(1,3,1)
plot(w,abs(fftshift(fft(EEG_ds3_1+EEG_ds3_2+EEG_ds3_3+EEG_ds3_4+EEG_ds3_5+EEG_ds3_6+EEG_ds3_7+EEG_ds3_8,N))))
title('FFT N/8 points with zero padding')
xlabel('F(Hz)')
ylabel('abs(FFT)')

EEG_ds2_1=[EEG_ds(1:N/4) zeros(1,192)];
EEG_ds2_2=[EEG_ds(N/4+1:2*N/4) zeros(1,192)];
EEG_ds2_3=[EEG_ds(2*N/4+1:3*N/4) zeros(1,192)];
EEG_ds2_4=[EEG_ds(3*N/4+1:4*N/4) zeros(1,192)];
subplot(1,3,2)
plot(w,abs(fftshift(fft(EEG_ds2,N))))
title('FFT N/4 points with zero padding')
xlabel('F(Hz)')
ylabel('abs(FFT)')

EEG_ds1_1=[EEG_ds(1:N/2) zeros(1,128)];
EEG_ds1_2=[EEG_ds(N/2+1:N) zeros(1,128)];

subplot(1,3,3)
plot(w,abs(fftshift(fft(EEG_ds1,N))))
xlabel('F(Hz)')
ylabel('abs(FFT)')
title('FFT N/2 points with zero padding')

%% Functions
function out=spectrogram_fft(y,L,overlap,nfft)
col=length(y)/(L-overlap);
spect=zeros(L/2,floor(col));
index=L/2;
for i=1:(floor(col)-1)
    signal=y(index-L/2+1:index+L/2).*rectwin(L).';
    fft_signal=fft(signal,nfft);
    spect(:,i)=fft_signal(1:nfft/2).';
    index=index+(L/2);
end
surf((20*log(abs(spect))/nfft))
%imagesc(flipud((20*log(abs(spect))/nfft)))
xlabel('time(s)')
ylabel('frequency')
title('spectrogram using fft')
end