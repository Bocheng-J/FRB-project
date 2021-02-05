fsamp=500e6;        %fsamp=250MHz 输入采样频率为250MHz
Ts=1/fsamp;        %Ts为fsamp的倒数 即输入采样间隔Ts
band=200e6;         %预设的采样带宽为200MHz
Tp=60e-6;          %预设的采样时间周期Tp为60us
N=Tp*fsamp;        %N为输入采样频率与采样时间周期之积。表示在采样时间周期Tp内，以fsamp的采样率采样可以得到的采样点数   N = 15000
u=band/Tp;         %u为带宽除以时宽。表示在单位时间间隔内的频带宽度。也即这30M的带宽分布在Tp=60us的时间周期上，单位时间的频带宽度
t=-Tp/2:Tp/N:Tp/2-Tp/N;       %t取点从-Tp/ b b2开始以Tp/N为步进值增加到Tp/2-Tp/N。
f0=125e6;           %输入的已调频信号载波频率为125MHz
xs=cos(2*pi*(f0*t+0.5*u*t.^2));        %输入的已调频信号经fsamp=500MHz带通采样后的输出。相当于A/D转换后的数字信号序列
round(xs*2^11-1);%
S0=fft(xs,N);      %S0是对A/D转换后的数字信号序列进行N点fft的结果；
S1=abs(S0);      %S1是对s0求模的结果；
S2=(S1);
S2=awgn(S2,10);    %在S2的频谱中加入10db的高斯白噪声
f=0:fsamp/N:fsamp-fsamp/N;       %f的取点由0开始以fsamp/N为步进值直到fsamp-fsamp/N结束
figure(1);         %画出载频为125MHz有用信号带宽为200MHz的带通信号序列频谱
plot(f/1e6,20*log10(S2/max(S2)));           %横轴以MHz为单位，纵轴以dB形式，其中S2/max(S2)表示输出该带通信号序列相对幅度大小，对它取对数后的结果就是dB的形式了。
title('载频为125MHz有用信号带宽为200MHz的带通信号序列频谱');
xlabel('frequency(MHz)');                   %从频谱图中可以看出信号序列频谱有（11MHz，41MHz）和（55MHz，85MHz）两部分。
ylabel('Magnitude(dB)');                    %这是因为经fsamp=94MHz带通采样时，在fsamp/2=48MHz处发生了频谱折叠，原来信号序列频谱（55M，85M）折叠到（11MHz，41MHz）了。
grid on                                     %这两部分频谱形状一致，没有发生频谱混叠。这里的带通采样速率fsamp=96MHz是通过计算得出来的。
                                            %具体计算式如下：fsamp>=4f0/(2n+1) 且fsamp>=2B。这里f0指信号序列中心频率，B指信号序列带宽 n取正整数。

%NCO数控振荡器模块%
for t=1:N
    t1=(t-1)*Ts;
    ncoi_c(t)=cos(2*pi*f0*t1);             %产生频率为f0的cos数控本振（I路），这里产生数控本振的时间间隔与A/D采样间隔相同，便于序列后续相乘。
end

for t=1:N
    t1=(t-1)*Ts;
    ncoq_c(t)=sin(2*pi*f0*t1);            %产生频率为f0的sin数控本振（Q路），这里产生数控本振的时间间隔与A/D采样间隔相同，便于序列后续相乘。
end
ncoi=awgn(ncoi_c,80);         %对产生的I路本振序列加入80dB的高斯白噪声
ncoq=awgn(ncoq_c,80);         %对产生的Q路本振序列加入80dB的高斯白噪声

f=0:fsamp/N:fsamp-fsamp/N;    %f取点从0开始以fsamp/N为步进值直到fsamp-fsamp/N结束
u1=abs(fft(ncoi));            %对加入高斯白噪声的I路本振信号序列进行FFT后取模
u2=abs(fft(ncoq));            %对加入高斯白噪声的Q路本振信号序列进行FFT后取模
figure(2);                    %画出加入高斯白噪声后的数控振荡器I路信号频谱
plot(f/1e6,20*log10(u1/max(u1)));%横轴以MHz为单位，纵轴是dB形式
title('加入高斯白噪声的数控振荡器I路信号频谱');
xlabel('frequency(MHz)');
ylabel('Magnitude(dB)');
grid on
figure(3);                     %画出加入高斯白噪声后的数控振荡器Q路信号频谱
plot(f/1e6,20*log10(u2/max(u2)));%横轴以MHz为单位，纵轴是dB形式
title('加入高斯白噪声的数控振荡器Q路信号频谱');
xlabel('frequency(MHz)');
ylabel('Magnitude(dB)');
grid on

for n=1:1:N
    ysi(n)=xs(n)*ncoi(n);%A/D带通采样后信号序列与数字本振I路信号序列混频相乘（下变频）
    ysq(n)=xs(n)*ncoq(n);%A/D带通采样后信号序列与数字本振Q路信号序列混频相乘（下变频）
end