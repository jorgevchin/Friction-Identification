clear all;
close all;
clc;
set(0,'defaultTextInterpreter','latex');

% Comparision of q1
% Load the experimental data
file = strcat('pos1_20_exp2_gan.t');
fileID = fopen(file,'r');
A = textscan(fileID,'%f %f');
fclose('all');

q1 = A{2};
t = A{1};

% Load the simulation data
q1_dahl = load('q1_dahl.mat');
q1_lugre = load('q1_lugre.mat');
q1_wof = load('q1_wof.mat');
q1_vcs = load('q1_vcs.mat');

% Creates and saves a graphic with the experiment ans simulation
figure(1)
plot(t,q1,'k','LineWidth',2);
hold on
plot(q1_dahl.ans.Time,q1_dahl.ans.Data*180/pi,'b','LineWidth',2);
plot(q1_lugre.ans.Time,q1_lugre.ans.Data*180/pi,'r','LineWidth',2);
plot(q1_vcs.ans.Time,q1_vcs.ans.Data*180/pi,'g','LineWidth',2);
plot(q1_wof.ans.Time,q1_wof.ans.Data*180/pi,'m','LineWidth',2);
xlim([0 20])
grid on
legend({'Experiment','Dahl','LuGre','S+V','without friction'},'Location','northoutside','Orientation','horizontal','FontSize',10)
set(gca,'FontSize',16)
xlabel('time [s]','fontsize',16);
ylabel('$q_1$ [$^\circ$]','fontsize',16);
saveas(gcf,'q1_comp','epsc')

% Calculates the root-mean square value of the error between the
% simulations with frioction and the experiment
for i = 1:length(q1)
    q1e_dahl(i) = q1(i) - q1_dahl.ans.Data(ceil(t(i)*10000))*180/pi;
    q1e_lugre(i) = q1(i) - q1_lugre.ans.Data(ceil(t(i)*10000))*180/pi;
    q1e_vcs(i) = q1(i) - q1_vcs.ans.Data(ceil(t(i)*10000))*180/pi;
end

disp('Valores rms para q1')
disp('dahl =')
rms(q1e_dahl)
disp('lugre =')
rms(q1e_lugre)
disp('vcs')
rms(q1e_vcs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparision of q2
% Load the experimental data
file = strcat('pos2_2_exp2_gan.t');
fileID = fopen(file,'r');
A = textscan(fileID,'%f %f');
fclose('all');

q2 = A{2};
t = A{1};

% Load the simulation data
q2_dahl = load('q2_dahl.mat');
q2_lugre = load('q2_lugre.mat');
q2_wof = load('q2_wof.mat');
q2_vcs = load('q2_vcs.mat');

% Creates and saves a graphic with the experiment ans simulation
figure(2)
plot(t,q2,'k','LineWidth',2);
hold on
plot(q2_dahl.ans.Time,q2_dahl.ans.Data*180/pi,'b','LineWidth',2);
plot(q2_lugre.ans.Time,q2_lugre.ans.Data*180/pi,'r','LineWidth',2);
plot(q2_vcs.ans.Time,q2_vcs.ans.Data*180/pi,'g','LineWidth',2);
plot(q2_wof.ans.Time,q2_wof.ans.Data*180/pi,'m','LineWidth',2);
xlim([0 20])
ylim([-40 40])
grid on
legend({'Experiment','Dahl','LuGre','S+V','without friction'},'Location','northoutside','Orientation','horizontal','FontSize',10)
set(gca,'FontSize',16)
xlabel('time [s]','fontsize',16);
ylabel('$q_2$ [$^\circ$]','fontsize',16);
saveas(gcf,'q2_comp','epsc')

% Calculates the root-mean square value of the error between the
% simulations with friction and the experiment
for i = 1:length(q2)
    q2e_dahl(i) = q2(i) - q2_dahl.ans.Data(ceil(t(i)*10000))*180/pi;
    q2e_lugre(i) = q2(i) - q2_lugre.ans.Data(ceil(t(i)*10000))*180/pi;
    q2e_vcs(i) = q2(i) - q2_vcs.ans.Data(ceil(t(i)*10000))*180/pi;
end

disp('Valores rms para q2')
disp('dahl =')
rms(q2e_dahl)
disp('lugre =')
rms(q2e_lugre)
disp('vcs')
rms(q2e_vcs)