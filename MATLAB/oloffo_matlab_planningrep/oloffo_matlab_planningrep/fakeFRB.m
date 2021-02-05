%% Simulated FRB. Only dispersion
clf, clc
%Pulse profile
pulse_width = 1; %[ms], top-hat width of pulse
DM = 1500; %[cm^-3 pc], dispersion measure

%Calculate maximum time delay
deltaT_max = 4.148808*(1.2^-2 - 1.7^-2)*DM; %[ms]

%Create window grid coordinates
dt = 40*10^-3; %[ms], grid time resolution
t_max = 1.1*deltaT_max; %[ms], window size
t_vec = linspace(0, t_max, t_max/dt); %[ms]

f_min = 1.2; %[GHz], window lower freq
f_max = 1.7; %[GHz], window upper freq
df = 125*10^-6; %[GHz], grid frequency resolution
%Nf = (f_max-f_min)/df; %#frequency channels
Nf = 2000;
f_vec = linspace(f_min, f_max, Nf); %[GHz]

[t_mat, f_mat] = meshgrid(t_vec, f_vec); %grid coordinates


%Quadratically dispersed wideband pulse
dispersion = 4.148808*(f_min.^-2 - f_mat.^-2)*DM; %[ms]
pulse_offset = max(max(dispersion)); %[ms]
quad_pulse = abs(t_mat - pulse_offset + dispersion) < pulse_width/sqrt(2);

figure(1)
imagesc(t_vec, f_vec, quad_pulse)
set(gca, 'YDir', 'normal')
title(['Dispersed FRB dynamic spectrum. DM = ', num2str(DM), ' cm^{-3}pc, W = ', num2str(pulse_width), ' ms'])
xlabel('[ms]');
ylabel('[GHz]');

figure(2)
deltaT = @(f) deltaT_max - 4.148808.*(f_min.^-2 - f.^-2).*DM;
f = linspace(f_min, f_max, Nf*10);
plot(deltaT(f), f)
title(['Dispersed FRB profile. DM = ', num2str(DM), ' cm^{-3}pc'])
xlabel('[ms]');
ylabel('[GHz]');
grid on


disp('Done!')