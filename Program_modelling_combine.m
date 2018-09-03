%Program to derive MFCC coefficients when a .wav file is input
clearvars
clc

%Input
filename1 = 'D:\Acads\IDP-sem7\data_solo\data1\solo.tar\solo\data\solo\irf02\irf02_f01_solo.wav'; %.wav solo speech file path
filename2 = 'D:\Acads\IDP-sem7\data_whisper\whsp.tar\whsp\data\whsp\irf02\irf02_f03_whsp.wav'; %.wav whisper speech file path
filename3 = 'D:\Acads\IDP-sem7\data_fast\fast.tar\fast\data\fast\irf02\irf02_f04_fast.wav'; %.wav fast speech file path
[signal1,Fs1] = audioread(filename1);%Read the solo speech signal into matrix
[signal2,Fs2] = audioread(filename2);%Read the whisper speech signal into matrix
[signal3,Fs3] = audioread(filename3);%Read the whisper speech signal into matrix
N1 = size(signal1,1);%Size of the matrix
N2 = size(signal2,1);
N3 = size(signal3,1);
%Input Parameters
preemph_c = 0.97; %Preemphasis fiter coefficicent
frame_size = 0.010; %Frame Length in seconds
frame_shift = 0.006; %Distance between left edges of successive windows in seconds
nfft = 512; %Length of fft
f_min = 0; %lower frequency for mel filter
f_max = Fs1/2; %upper frequency for mel filter
mel_n = 40; %number of mel filters
coeff_n  = 13; %number of coefficients

%Code parameters
size_n = round(Fs1*frame_size); %Number of samples in a window
shift_n = round(Fs1*frame_shift); %Number of samples between left edges of sucesssive windows

%Preemphasis
signal1 = filter([1 -preemph_c], 1, signal1); %Applying filter to emphasise the high frequencies
signal2 = filter([1 -preemph_c], 1, signal2); 
signal3 = filter([1 -preemph_c], 1, signal3); 
%Framing
frames1 =  buffer(signal1, size_n, size_n-shift_n,'nodelay');%Divide the speech signal into frames 
frames1 = frames1(:,1:end-1);%Remove the last frame as it has zeros padded
frames2 =  buffer(signal2, size_n, size_n-shift_n,'nodelay');%Divide the speech signal into frames 
frames2 = frames2(:,1:end-1);%Remove the last frame as it has zeros padded
frames3 =  buffer(signal3, size_n, size_n-shift_n,'nodelay');%Divide the speech signal into frames 
frames3 = frames3(:,1:end-1);%Remove the last frame as it has zeros padded
%Hamming Window
w = hamming(size_n);%Finding Hamming window of required size
frames1 = bsxfun(@times,frames1,w); %Apply hamming window to make smooth frame edges
frames_n1 = size(frames1);%Number of frames
frames2 = bsxfun(@times,frames2,w); %Apply hamming window to make smooth frame edges
frames_n2 = size(frames2);%Number of frames
frames3 = bsxfun(@times,frames3,w); %Apply hamming window to make smooth frame edges
frames_n3 = size(frames3);%Number of frames
%Finding FFT of each frame
frames_fft1 = fft(frames1,nfft,1);%Finding fft of each frame
frames_power1 = ((abs(frames_fft1)).^2)/(size_n);%Finding power of each frame at each frequency
frames_fft2 = fft(frames2,nfft,1);%Finding fft of each frame
frames_power2 = ((abs(frames_fft2)).^2)/(size_n);
frames_power2 = frames_power2(:,1:floor(size(frames_power1,2)/5));
frames_fft3 = fft(frames3,nfft,1);%Finding fft of each frame
frames_power3 = ((abs(frames_fft3)).^2)/(size_n);
frames_power3 = frames_power3(:,1:floor(size(frames_power1,2)/5));
%Finding Mel Filters
[mel_filters,mel_n] = MelFilter(f_min,f_max,mel_n,nfft,Fs1);%Generating the mel filter banks

%Finding the coefficients
coeff1 = (mel_filters.')*frames_power1;%Passing the Power into the mel filters
coeff2 = (mel_filters.')*frames_power2;
coeff3 = (mel_filters.')*frames_power3;

%Finding DCT of log
dctm = @(N,M)(sqrt(2.0/M)*cos(repmat([0:N-1]',1,M).*repmat(pi*([1:M]-0.5)/M,N,1)));%Definfing DCT
coeff1 = dctm(coeff_n,mel_n)*log(coeff1);%Appying Dct on the coeffiients
coeff1(isnan(coeff1)) = 0;%Make any undefined value to zero
coeff1 = coeff1(2:end,:);%Remove the first coefficient
coeff1 = coeff1(:,1:floor(size(coeff1,2)/2));
coeff2 = dctm(coeff_n,mel_n)*log(coeff2);%Appying Dct on the coeffiients
coeff2(isnan(coeff2)) = 0;%Make any undefined value to zero
coeff2 = coeff2(2:end,:);
coeff3 = dctm(coeff_n,mel_n)*log(coeff3);%Appying Dct on the coeffiients
coeff3(isnan(coeff3)) = 0;%Make any undefined value to zero
coeff3 = coeff3(2:end,:);

%save('D:\Acads\IDP-sem7\codes\testing_solo2\coeff20','coeff')
[~,model1] = EM_gmm(coeff1,12);
[~,model2] = EM_gmm(coeff2,12);
[~,model3] = EM_gmm(coeff3,12);
%[z2,model2,llh2] = mixGaussEm(coeff,8);
save('D:\Acads\IDP-sem7\codes\model_combine_12_new\model10','model1','model2','model3');
%save('D:\Acads\IDP-sem7\codes\model_8\model20','model2');
