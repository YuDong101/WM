
function [f,mag] = fun_FFT(input)

XX_FFT=input';

[N,m]=size(XX_FFT);

fs=1000;
n=0:N-1;   
y=fft(XX_FFT,N);    %对信号进行快速Fourier变换
mag=abs(y').^2/N^2;     %求得Fourier变换后的振幅
f=n*fs/N;

save date.mat