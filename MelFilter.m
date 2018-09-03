function [ mel_filters, mel_n ] = MelFilter( f_min,f_max,mel_n,nfft,Fs )
%Function to find the mel filters. 
%input-----------------------------------------------
%f_min: lower frequency for mel filter
%f_max: upper frequency for mel filter
%mel_n: number of mel filters
%nfft: length of fft of the frames 
%Fs: Sampling rate of speech signal
%output----------------------------------------------
%melfilters: matrix f size (nfft x mel_n). Each column is a mel filter.
%Code------------------------------------------------
mel_fmin = 1125*log(1+(f_min/700)); %Min freq in mel scale
mel_fmax = 1125*log(1+(f_max/700)); %Max freq in mel scale
mel_freq = linspace(mel_fmin, mel_fmax, mel_n+2); %Divide equally in Mel Scale
freq = 700*(exp(mel_freq/1125)-1); %Convert back to freq from mel scale 
bin_fft = floor((nfft+1)*freq/Fs) + 1;%Find fft bin for each frequency
bin_fft = unique(bin_fft);%Take only the unique fft bins
temp = size(bin_fft);%Size of the unique fft bins
mel_n = temp(2)-2;%Update number of filters according to number of unique fft bins
mel_filters = zeros(nfft,mel_n);%Initialise
for i = 2:mel_n+1
    for j = bin_fft(i-1):bin_fft(i)
       mel_filters(j,i-1) = (j - bin_fft(i-1))/(bin_fft(i) - bin_fft(i-1));
    end
    for j = bin_fft(i)+1:bin_fft(i+1)
        mel_filters(j,i-1) = (bin_fft(i+1)-j)/(bin_fft(i+1) - bin_fft(i));
    end
    plot(mel_filters(:,i-1));
    hold on;
end

