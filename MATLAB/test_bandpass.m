fs = 4000;
t = 0:1/fs:1-1/fs;
y = sin(2*pi*200*t);
figure;
plot(y);

bandpassY = bandpass(y,[100 300],fs);
figure;
plot(bandpassY);