clc;
close all;
clear all;
d1=importdata('typeA_4.txt');  %h6
%%
x=d1.data(:,2);       %acceleration
y=d1.data(:,3);       %acceleration (big)
z=d1.data(:,4);       %acceleration new 
p=d1.data(:,5);       %point acceleration
time = d1.data(:,1);   %time
a= 20;
b= 30;
fs=4096;
x=x(a*fs:b*fs);
y=y(a*fs:b*fs);
z=z(a*fs:b*fs);
p=p(a*fs:b*fs);
time= time(a*fs:b*fs);
%
sen_acc1 = 100*10^(-3);     %acc1
sen_acc2 = 1;     %acc2 
sen_acc3 = 100*10^(-3);         %acc3
sen_acc4 = 1;   
%acc4
acc1=x/sen_acc1; 
acc2= y/sen_acc2;           %acceleration=output data/sensitivity of accelerator
acc3= z/sen_acc3;
acc4= p/sen_acc4;
figure(1)
plot(time,acc1,'k', 'linewidth',1.5);
xlabel('Time (ms)','fontsize',12,'FontWeight','bold');
ylabel('Amplitude(N)','fontsize',12,'FontWeight','bold');
title('Acceleration signal on the foot 4','fontsize',18,'FontWeight','bold');
grid minor;
figure(2)
plot(time,acc2,'k', 'linewidth',1.5);
xlabel('Time (ms)','fontsize',12,'FontWeight','bold');
ylabel('Amplitude(m/s^2)','fontsize',12,'FontWeight','bold');
title('Acceleration signal on the floor near foot 4','fontsize',18,'FontWeight','bold');
grid minor
figure(3)
plot(time,acc3,'k', 'linewidth',1.5);
xlabel('Time (ms)','fontsize',12,'FontWeight','bold');
ylabel('Amplitude(m/s^2)','fontsize',12,'FontWeight','bold');
title('Acceleration signal on the foot 7','fontsize',18,'FontWeight','bold');
grid minor;
figure(4)
plot(time,acc4,'k', 'linewidth',1.5);
xlabel('Time (ms)','fontsize',12,'FontWeight','bold');
ylabel('Amplitude(m/s^2)','fontsize',12,'FontWeight','bold');
title('Acceleration signal on the floor near foot 7','fontsize',18,'FontWeight','bold');
grid minor;
%% Power spectral density
win_len = 2^13;
win1 = rectwin(win_len);
Fs1 = 4096;
t = win_len/Fs1;
%%
[psd_x f_x]= cpsd(acc1,acc1,win1,0.5*win_len,win_len,Fs1);
[psd_y f_y]= cpsd(acc2,acc2,win1,0.5*win_len,win_len,Fs1);
[psd_z f_z]= cpsd(acc3,acc3,win1,0.5*win_len,win_len,Fs1);
[psd_p f_p]= cpsd(acc4,acc4,win1,0.5*win_len,win_len,Fs1);
%%
figure(5)
plot(f_x,10*log10(psd_x), 'linewidth',1.5);
hold on;
plot(f_y,10*log10(psd_y),'k', 'linewidth',1.5);
hold on,
plot(f_z,10*log10(psd_z),'g', 'linewidth',1.5);
hold on;
plot(f_p,10*log10(psd_p*t),'r', 'linewidth',1.5);
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold')
ylabel('Acceleration PSD dB rel [(m/s^2)^2*s]/Hz','fontsize',12,'FontWeight','bold')
title('Power spectral densities for Acceleration','fontsize',18,'FontWeight','bold');
grid on;
xlim([0 500])
legend('foot 3','floor near foot 3','foot 7','floor near foot 7');
%% Cross Power Spectral Density
[csd_xy f_xy]= cpsd(acc1,acc2,win1,0.5*win_len,win_len,Fs1);
[csd_yz f_yz]= cpsd(acc2,acc3,win1,0.5*win_len,win_len,Fs1);
[csd_xz f_xz]= cpsd(acc1,acc3,win1,0.5*win_len,win_len,Fs1);
[csd_xp f_xp]= cpsd(acc1,acc4,win1,0.5*win_len,win_len,Fs1);
[csd_yp f_yp]= cpsd(acc2,acc4,win1,0.5*win_len,win_len,Fs1);
[csd_zp f_zp]= cpsd(acc3,acc4,win1,0.5*win_len,win_len,Fs1);
% a = angle(csd_xy);
% figure
% plot(f_xy,a);
figure
plot(f_xy,10*log10(csd_xy),'k', 'linewidth',1.5);
hold on,
plot(f_yz,10*log10(csd_yz),'r', 'linewidth',1.5);
hold on;
plot(f_xz,10*log10(csd_xz),'g', 'linewidth',1.5);
hold on;
plot(f_xp,10*log10(csd_xp),'m', 'linewidth',1.5);
hold on;
plot(f_yp,10*log10(csd_yp),'b', 'linewidth',1.5);
hold on;
plot(f_zp,10*log10(csd_zp),'y', 'linewidth',1.5);
xlim([5 250])
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold');
ylabel('CSD dB rel [(Nm/s^2)^2]/Hz','fontsize',12,'FontWeight','bold');
title('Cross spectral densities of Force and Acceleration','fontsize',18,'FontWeight','bold');
grid on;

%% Theoritical FRF calculation
% m=3200;    
% fd=274;
% % br=5.1127;
% br = 0.82;       %-3 dB frequency.
% dam=br/(2*fd);
% wr=2*pi*fd;
% k=(wr^2)*m;
% c=2*dam*(sqrt(m*k));
% f=0:0.1:200;
% 
% for j=0: length(f)-1
%     %Hu_f(j+1)=[1/k]/[1-((f(j+1)/fd).^2)+(i*2*dam*(f(j+1)/fd))];
%     %Ha_f(j+1)= (-(2*pi*f(j+1)).^2)*Hu_f(j+1);
% 
% %Ha_f(j+1)=[(-(2*pi*f(j+1)).^2)/k]/[1-((f(j+1)/fd).^2)+(i*2*dam*(f(j+1)/fd))];
% Ha_f(j+1)=[1+j*2*dam*(f(j+1)/fd)]/[1-((f(j+1)/fd).^2)+(i*2*dam*(f(j+1)/fd))];
% end
%% Transfer function
[trans_fun trans_f]= tfestimate(acc1,acc2,win1,0.5*win_len,win_len,Fs1);
[trans_fun2 trans_f2]= tfestimate(acc3,acc4,win1,0.5*win_len,win_len,Fs1);
figure
plot(trans_f,20*log10(abs(trans_fun)),'k', 'linewidth',1.5);
hold on;
plot(trans_f2,20*log10(abs(trans_fun2)),'r', 'linewidth',1.5);
hold on;
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold');
ylabel('Magnitude of transfer function estimate','fontsize',12,'FontWeight','bold');
title('FRF plot of acceleration signals','fontsize',18,'FontWeight','bold');
legend('foot 3(A1-A2)','foot 6(A4-A3)');
xlim([0 500])
grid on;
% grid minor;

%% coherence
[Cxy1 f_c1]= mscohere(acc1,acc2,win1,0.5*win_len,win_len,Fs1);
[Cxy2 f_c2]= mscohere(acc3,acc4,win1,0.5*win_len,win_len,Fs1);
% [Cxy3 f_c3]= mscohere(acc2,acc4,win1,0.5*win_len,win_len,Fs1);
% [Cxy4 f_c4]= mscohere(acc3,acc1,win1,0.5*win_len,win_len,Fs1);
figure
plot(f_c2,Cxy2,'k','linewidth',1.5);
hold on;
plot(f_c1,Cxy1,'r', 'linewidth',1.5);
hold on;
axis([0 500 0 1])
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold')
ylabel('Coherence','fontsize',12,'FontWeight','bold')
title('Coherence function','fontsize',18,'FontWeight','bold');
legend('Foot 3(A1-A2)','Foot 6(A3-A4)');
grid on;
%axis([0 250 0 1])