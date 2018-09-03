%Program to derive LPCC coefficients when a .wav file is input
clearvars
clc

%Input
filename = 'D:\Acads\IDP-sem7\data_solo\data1\solo.tar\solo\data\solo\irf02\irf02_f01_solo.wav'; %.wav file path
[signal,Fs] = audioread(filename);%Read the speech signal into matrix
N = size(signal);%Size of the matrix

%Input Parameters
frame_size = 0.010; %Frame Length in seconds
frame_shift = 0.006; %Distance between left edges of successive windows in seconds
k = 10;%Order of the LPC filter to derive
coeff_n  = 15; %number of coefficients required per frame

%Code parameters
N = N(1);%Number of samples in the input speech signal
size_n = round(Fs*frame_size); %Number of samples in a window
shift_n = round(Fs*frame_shift); %Number of samples between left edges of sucesssive windows

%Framing
frames =  buffer(signal, size_n, size_n-shift_n,'nodelay');%Divide the speech signal into frames 
frames = frames(:,1:end-1);%Remove the last frame as it has zeros padded

%Finding the coefficients
coeff_LPC = zeros(k,size(frames,2));%Initialising the LPC matrix
for i = 1:size(frames,2)
    coeff_LPC(:,i) = levinson_recursion( frames(:,i), k );%Finding each frame's LPC by levinson recursion
end

%Finding LPCC matrix from LPC matrix
coeff = zeros(coeff_n,size(frames,2));%Initialising the LPCC matrix
help1 = [1:coeff_n];
for i = 1:size(frames,2)
    for j = 1:coeff_n
        if(j<=k)
            temp = (help1(1:j-1))'/j;
            help2 = (coeff(1:j-1,i)).*flipud(coeff_LPC(1:j-1,i));
            help2 = temp.*help2;
            coeff(j,i) = -coeff_LPC(j,i) - sum(help2);
        else
            temp = (help1(j-k:j-1))'/j;
            help2 = (coeff(j-k:j-1,i)).*flipud(coeff_LPC(1:k,i));
            help2 = temp.*help2;
            coeff(j,i) = -sum(help2);
        end
    end
end
coeff = coeff(:,1:floor(size(coeff,2)/5));
%size(coeff)
%save('D:\Acads\IDP-sem7\codes\testing_solo_LPCC\coeff10','coeff')
[~,model1] = EM_gmm(coeff,8);
[~,model2] = EM_gmm(coeff,12);
[~,model3] = EM_gmm(coeff,16);
%[z2,model2,llh2] = mixGaussVb(coeff,20);
save('D:\Acads\IDP-sem7\codes\model_8_LPCC_15\model10','model1');
save('D:\Acads\IDP-sem7\codes\model_12_LPCC_15\model10','model2');
save('D:\Acads\IDP-sem7\codes\model_16_LPCC_15\model10','model3');
