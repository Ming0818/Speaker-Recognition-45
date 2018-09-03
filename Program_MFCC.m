%Program to derive MFCC coefficients when a .wav file is input
clearvars
clc

%Input
filename = 'D:\Acads\IDP-sem7\data_fast\fast.tar\fast\data\fast\irf02\irf02_f01_fast.wav'; %.wav file path
[signal,Fs] = audioread(filename);%Read the speech signal into matrix
N = size(signal);%Size of the matrix

%Input Parameters
preemph_c = 0.97; %Preemphasis fiter coefficicent
frame_size = 0.010; %Frame Length in seconds
frame_shift = 0.006; %Distance between left edges of successive windows in seconds
nfft = 512; %Length of fft
f_min = 0; %lower frequency for mel filter
f_max = Fs/2; %upper frequency for mel filter
mel_n = 40; %number of mel filters
coeff_n  = 13; %number of coefficients

%Code parameters
N = N(1);%Number of samples in the input speech signal
size_n = round(Fs*frame_size); %Number of samples in a window
shift_n = round(Fs*frame_shift); %Number of samples between left edges of sucesssive windows

%Preemphasis
signal = filter([1 -preemph_c], 1, signal); %Applying filter to emphasise the high frequencies

%Framing
frames =  buffer(signal, size_n, size_n-shift_n,'nodelay');%Divide the speech signal into frames 
frames = frames(:,1:end-1);%Remove the last frame as it has zeros padded

%Hamming Window
w = hamming(size_n);%Finding Hamming window of required size
frames = bsxfun(@times,frames,w); %Apply hamming window to make smooth frame edges
frames_n = size(frames);%Number of frames

%Finding FFT of each frame
frames_fft = fft(frames,nfft,1);%Finding fft of each frame
frames_power = ((abs(frames_fft)).^2)/(size_n);%Finding power of each frame at each frequency

%Finding Mel Filters
[mel_filters,mel_n] = MelFilter(f_min,f_max,mel_n,nfft,Fs);%Generating the mel filter banks

%Finding the coefficients
coeff = (mel_filters.')*frames_power;%Passing the Power into the mel filters

%Finding DCT of log
dctm = @(N,M)(sqrt(2.0/M)*cos(repmat([0:N-1]',1,M).*repmat(pi*([1:M]-0.5)/M,N,1)));%Definfing DCT
coeff = dctm(coeff_n,mel_n)*log(coeff);%Appying Dct on the coeffiients
coeff(isnan(coeff)) = 0;%Make any undefined value to zero
coeff = coeff(2:end,:);%Remove the first coefficient
%coeff = coeff(:,1:floor(size(coeff,2)/5));
%filename_model_save = 'D:\Acads\IDP-sem7\codes\testing_3_MFCC\coeff01';
%save(filename_model_save,'coeff')
k = 8;
[~,model] = EM_gmm(coeff,k);
%filename_save = 'D:\Acads\IDP-sem7\codes\model_fast_8\model10';
%save(filename_save,'model');

