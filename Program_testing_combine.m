%Code for speaker identification
clearvars;
clc;
%Input
N  = 10; %Number of speakers
folder1 = dir('D:\Acads\IDP-sem7\codes\model_combine_12_new'); %Directory having the speaker models
%Loading the speaker models
for i =1:N
    fname=folder1(i+2).name;
    fpath=strcat('D:\Acads\IDP-sem7\codes\model_combine_12_new\',fname);
    load(fpath);
    S1(i) = model1;%Holds all the speaker models
    S2(i) = model2;
    S3(i) = model3;
end
%Testing files
folder2 =dir('D:\Acads\IDP-sem7\codes\testing_solo_MFCC');%Directory having the test speech samples
Fs = 44100; %Sampling rate of the speech samples
T = [0.1 0.3 0.5 0.8 1 2 5 7 10];%Test sample length in seconds
correct_count = zeros(N,numel(T));
vector_rate = 100; %Number of feature vector per second
M = 1000; %Number of samples to test for
for i = 1:N
    fname=folder2(i+2).name;
    fpath=strcat('D:\Acads\IDP-sem7\codes\testing_solo_MFCC\',fname);
    load(fpath);%Loading the feature vectors of testing speech samples one by one
    loglikelihood = zeros(N,M,numel(T));%Matrix storing the log likelihood for one speaker for all the speaker models and all the T
    for j = 1:numel(T)
        n = T(j)*vector_rate;%Sample length
        for k = 1: N
            loglikelihood(k,:,j) = loglikelihood_cal_combine(coeff(:,1:M+n-1),S1(k),S2(k),S3(k),M,n);%Function to find the log likelihood
        end
        [Maximum,Index] = max(loglikelihood(:,:,j),[],1);%Finding arg max
        correct_count(i,j) = nnz(Index==i);%Counting number of correct answers
    end     
end
correct_per = (correct_count*100)/M;
mean_correct_percent = mean(correct_per,1);
