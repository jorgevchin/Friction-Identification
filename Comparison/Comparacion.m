clear all;
close all;
clc;
set(0,'defaultTextInterpreter','latex');

file = strcat('pos1_20_exp2_gan.t');
fileID = fopen(file,'r');
A = textscan(fileID,'%f %f');
fclose('all');

q1_dahl = load('q1_dahl.mat');
q1_lugre = load('q1_lugre.mat');
q1_sin = load('q1_sin.mat');
q1_vcs2 = load('q1_vcs2.mat');

q1 = A{2};
t = A{1};

figure(1)
plot(t,q1,'k','LineWidth',2);
hold on
plot(q1_dahl.ans.Time,q1_dahl.ans.Data*180/pi,'b','LineWidth',2);
plot(q1_lugre.ans.Time,q1_lugre.ans.Data*180/pi,'r','LineWidth',2);
plot(q1_vcs2.ans.Time,q1_vcs2.ans.Data*180/pi,'g','LineWidth',2);
plot(q1_sin.ans.Time,q1_sin.ans.Data*180/pi,'m','LineWidth',2);
xlim([0 20])
grid on
legend({'Experiment','Dahl','LuGre','S+V','without friction'},'Location','northoutside','Orientation','horizontal','FontSize',10)
set(gca,'FontSize',16)
xlabel('time [s]','fontsize',16);
ylabel('$q_1$ [$^\circ$]','fontsize',16);
saveas(gcf,'q1_comp','epsc')

for i = 1:length(q1)
    q1e_dahl(i) = q1(i) - q1_dahl.ans.Data(ceil(t(i)*10000))*180/pi;
    q1e_lugre(i) = q1(i) - q1_lugre.ans.Data(ceil(t(i)*10000))*180/pi;
    q1e_vcs2(i) = q1(i) - q1_vcs2.ans.Data(ceil(t(i)*10000))*180/pi;
end

disp('Valores rms para q1')
disp('dahl =')
rms(q1e_dahl)
disp('lugre =')
rms(q1e_lugre)
disp('vc2')
rms(q1e_vcs2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = strcat('pos2_2_exp2_gan.t');
fileID = fopen(file,'r');
A = textscan(fileID,'%f %f');
fclose('all');

q2_dahl = load('q2_dahl.mat');
q2_lugre = load('q2_lugre.mat');
q2_sin = load('q2_sin.mat');
q2_vcs2 = load('q2_vcs2.mat');

q2 = A{2};
t = A{1};

figure(2)
plot(t,q2,'k','LineWidth',2);
hold on
plot(q2_dahl.ans.Time,q2_dahl.ans.Data*180/pi,'b','LineWidth',2);
plot(q2_lugre.ans.Time,q2_lugre.ans.Data*180/pi,'r','LineWidth',2);
plot(q2_vcs2.ans.Time,q2_vcs2.ans.Data*180/pi,'g','LineWidth',2);
plot(q2_sin.ans.Time,q2_sin.ans.Data*180/pi,'m','LineWidth',2);
xlim([0 20])
ylim([-40 40])
grid on
legend({'Experiment','Dahl','LuGre','S+V','without friction'},'Location','northoutside','Orientation','horizontal','FontSize',10)
set(gca,'FontSize',16)
xlabel('time [s]','fontsize',16);
ylabel('$q_2$ [$^\circ$]','fontsize',16);
saveas(gcf,'q2_comp','epsc')

for i = 1:length(q2)
    q2e_dahl(i) = q2(i) - q2_dahl.ans.Data(ceil(t(i)*10000))*180/pi;
    q2e_lugre(i) = q2(i) - q2_lugre.ans.Data(ceil(t(i)*10000))*180/pi;
    q2e_vcs2(i) = q2(i) - q2_vcs2.ans.Data(ceil(t(i)*10000))*180/pi;
end

disp('Valores rms para q2')
disp('dahl =')
rms(q2e_dahl)
disp('lugre =')
rms(q2e_lugre)
disp('vcs2')
rms(q2e_vcs2)