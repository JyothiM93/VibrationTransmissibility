clc;
close all;
clear all;
d1=importdata('typeA_4.txt');  %h6
%%
% d1= d1(20:90);
fs=4096;
a=180;
b=290;
x=d1.data(:,2);       %acceleration
x=x(a*fs:b*fs);
y=d1.data(:,3);       %acceleration (big)
y=y(a*fs:b*fs);
z=d1.data(:,4);       %acceleration new 
z=z(a*fs:b*fs);
p=d1.data(:,5);       %point acceleration
p=p(a*fs:b*fs);
time = d1.data(:,1);   %time
time= time(a*fs:b*fs);
%%
sen_acc1 = 10.2*10^(-3);     %acc1-cube V
sen_acc2 = 51.3*10^(-3);     %acc2-cyl V 
sen_acc3 = 1;%0.996*10^(-12);   %acc3- ch
sen_acc4 = 1;%3.091*10^(-12);   %acc4- ch
%acc4
acc1=x/sen_acc1; 
acc2= y/sen_acc2;           %acceleration=output data/sensitivity of accelerator
acc3= z/sen_acc3;
acc4= p/sen_acc4;
figure(1)
plot(time,acc1,'k', 'linewidth',1.5);
xlabel('Time (s)','fontsize',12,'FontWeight','bold');
ylabel('Amplitude(m/s^2)','fontsize',12,'FontWeight','bold');
title('Acceleration signal at foot-1 of machine 1','fontsize',18,'FontWeight','bold');
grid minor;
figure(2)
plot(time,acc2,'k', 'linewidth',1.5);
xlabel('Time (s)','fontsize',12,'FontWeight','bold');
ylabel('Amplitude(m/s^2)','fontsize',12,'FontWeight','bold');
title('Acceleration signal near floor of foot-1 of machine 1','fontsize',18,'FontWeight','bold');
grid minor
figure(3)
plot(time,acc3,'k', 'linewidth',1.5);
xlabel('Time (s)','fontsize',12,'FontWeight','bold');
ylabel('Amplitude(m/s^2)','fontsize',12,'FontWeight','bold');
title('Acceleration signal at foot of machine 2','fontsize',18,'FontWeight','bold');
grid minor;
figure(4)
plot(time,acc4,'k', 'linewidth',1.5);
xlabel('Time (s)','fontsize',12,'FontWeight','bold');
ylabel('Amplitude(m/s^2)','fontsize',12,'FontWeight','bold');
title('Acceleration signal near floor of foot of machine 2','fontsize',18,'FontWeight','bold');
grid minor;
%% Power spectral density
win_len = 2^12;
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
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold')
ylabel('Acceleration PSD dB rel [(m/s^2)^2*s]/Hz','fontsize',12,'FontWeight','bold')
title('Power spectral densities for Acceleration signals of machine 1','fontsize',18,'FontWeight','bold');
legend('machine 1-foot','machine1-floor');
grid on;
xlim([0 500])
figure(6)
plot(f_z,10*log10(psd_z),'g', 'linewidth',1.5);
hold on;
plot(f_p,10*log10(psd_p*t),'r', 'linewidth',1.5);
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold')
ylabel('Acceleration PSD dB rel [(m/s^2)^2*s]/Hz','fontsize',12,'FontWeight','bold')
title('Power spectral densities for Acceleration signals of machine 2','fontsize',18,'FontWeight','bold');
grid on;
xlim([0 500])
legend('machine 2-foot','machine2-floor');
%% Cross Power Spectral Density
[csd_1 f_1]= cpsd(acc1,acc2,win1,0.5*win_len,win_len,Fs1);
[csd_2 f_2]= cpsd(acc1,acc3,win1,0.5*win_len,win_len,Fs1);
[csd_3 f_3]= cpsd(acc3,acc4,win1,0.5*win_len,win_len,Fs1);
[csd_4 f_4]= cpsd(acc2,acc4,win1,0.5*win_len,win_len,Fs1);
% a = angle(csd_xy);
% figure
% plot(f_xy,a);
figure
plot(f_1,10*log10(csd_1),'k', 'linewidth',1.5);
hold on,
plot(f_2,10*log10(csd_2),'r', 'linewidth',1.5);
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold');
ylabel('CSD dB rel [(Nm/s^2)^2]/Hz','fontsize',12,'FontWeight','bold');
title('Cross spectral densities of Force and Acceleration','fontsize',18,'FontWeight','bold');
legend('machine1 foot1-machine1 floor','machine1 foot1-machine2 foot');
grid on;
figure
plot(f_3,10*log10(csd_3),'g', 'linewidth',1.5);
hold on;
plot(f_4,10*log10(csd_4),'m', 'linewidth',1.5);
xlim([5 250])
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold');
ylabel('CSD dB rel [(Nm/s^2)^2]/Hz','fontsize',12,'FontWeight','bold');
title('Cross spectral densities of Force and Acceleration','fontsize',18,'FontWeight','bold');
legend('machine2 foot-machine2 floor','machine1 floor-machine2 floor');
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
[trans_fun3 trans_f3]= tfestimate(acc2,acc4,win1,0.5*win_len,win_len,Fs1);
[trans_fun4 trans_f4]= tfestimate(acc3,acc1,win1,0.5*win_len,win_len,Fs1);
figure
plot(trans_f,20*log10(abs(trans_fun)),'k', 'linewidth',1.5);
hold on;
plot(trans_f2,20*log10(abs(trans_fun2)),'r', 'linewidth',1.5);
hold on;
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold');
ylabel('Magnitude of transfer function estimate','fontsize',12,'FontWeight','bold');
title('FRF plot of acceleration signals','fontsize',18,'FontWeight','bold');
legend('machine1 foot1-machine 1 floor','machine2 foot-machine2 floor');
xlim([0 500])
grid on;
figure
plot(trans_f3,20*log10(abs(trans_fun3)),'g', 'linewidth',1.5);
hold on;
plot(trans_f4,20*log10(abs(trans_fun4)),'m', 'linewidth',1.5);
xlim([0 500])
grid on;

xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold');
ylabel('Magnitude of transfer function estimate','fontsize',12,'FontWeight','bold');
title('FRF plot of acceleration signals','fontsize',18,'FontWeight','bold');
legend('machine1 floor- machine2 floor','machine1 foot1-machine2 foot');
% grid minor;

%% coherence
[Cxy1 f_c1]= mscohere(acc1,acc2,win1,0.5*win_len,win_len,Fs1);
[Cxy2 f_c2]= mscohere(acc3,acc4,win1,0.5*win_len,win_len,Fs1);
[Cxy3 f_c3]= mscohere(acc2,acc4,win1,0.5*win_len,win_len,Fs1);
[Cxy4 f_c4]= mscohere(acc3,acc1,win1,0.5*win_len,win_len,Fs1);
figure
plot(f_c1,Cxy1,'r', 'linewidth',1.5);
hold on;
plot(f_c2,Cxy2,'g', 'linewidth',1.5);
hold on;
axis([0 500 0 1])
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold')
ylabel('Coherence','fontsize',12,'FontWeight','bold')
title('Coherence function','fontsize',18,'FontWeight','bold');
legend('Machine1 foot1-machine1 floor','Machine2 foot-machine2 floor');
grid on;
figure
plot(f_c3,Cxy3,'m', 'linewidth',1.5);
hold on;
plot(f_c4,Cxy4,'b', 'linewidth',1.5);
xlabel('Frequency [Hz]','fontsize',12,'FontWeight','bold')
ylabel('Coherence','fontsize',12,'FontWeight','bold')
title('Coherence function','fontsize',18,'FontWeight','bold');
grid on;
legend('Machine1 floor-machine2 floor','Machine1 foot1-machine2 foot2');
axis([0 500 0 1])
%%
ZZ(:,1)=20*log10(abs(trans_fun));
ZZ(:,2)=Cxy1;
YY(:,1)=20*log10(abs(trans_fun2));
YY(:,2)=Cxy2;
