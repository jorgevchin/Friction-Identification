%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Static friction estimation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
set(0,'defaultTextInterpreter','latex');

% These arrays alawys to automate the load of the experimental results
var = char('pos1_','vel1_','tau1_','pos2_','vel2_','tau2_');
color = char('b','r','k','g');
expe = char('1','2','3','4');
pen = char('05','01');

for m  = 1:length(expe)  
    % Load the experimental results
    for j = 1:6
        if j > 3
            file = strcat('est_',var(j,:),pen(ceil(j/3),:),'_exp',expe(m)+4,'_gan.t');
        else
            file = strcat('est_',var(j,:),pen(ceil(j/3),:),'_exp',expe(m),'_gan.t');
        end
        fileID = fopen(file,'r');
        A = textscan(fileID,'%f %f');
        for k = 1:2
            for l = 1:length(A{1})
                switch (j)
                    case 1
                        q1(k,l) = A{k}(l);
                    case 2
                        q1Dot(k,l) = A{k}(l);
                    case 3
                        tau1(k,l) = A{k}(l);
                    case 4
                        q2(k,l) = A{k}(l);
                    case 5
                        q2Dot(k,l) = A{k}(l);
                    case 6
                        tau2(k,l) = A{k}(l);
                end
            end
        end
    end
    
    % Calculate the velocities using a seven-points center difference
    clear q1Dot
    clear q2Dot
    q1Dot(1,:) = q1(1,:);
    q2Dot(1,:) = q2(1,:);
    q1Dot(2,:) = numDiff(q1);    
    q2Dot(2,:) = numDiff(q2);
    
    % Creates a graphic with the results
    figure(1)
    plot(q1Dot(1,:),q1Dot(2,:),color(m),'LineWidth',2);
    hold on
    figure(2)
    plot(q2Dot(1,:),q2Dot(2,:),color(m),'LineWidth',2);
    hold on
    
    % Searchs for the first spike, using some threshold and calculates the
    % torque applied at the moment the spike appeared
    flag = 0;
    i = 1;
    while (flag==0)
        if (q1Dot(2,i) > 0.2)
            tau1_break(m) = q1Dot(1,i)*0.5;
            flag = 1;
        else
            i = i + 1;
        end
    end
    
    flag = 0;
    i = 1;
    while (flag==0)
        if (q2Dot(2,i) > 0.2)
            tau2_break(m) = q2Dot(1,i)*0.1;
            flag = 1;
        else
            i = i + 1;
        end
    end
    
    clear q1 q2 q1Dot q2Dot tau1 tau2    
end

% Calculates the static frcition coeffcients as the mean value of the
% results
fs1p = mean(tau1_break)
fs2p = mean(tau2_break)

% Adds a visual indication of the static frcition coeffcients to the
% graphics
figure(1)
xlim([0 25])
ylim([-1 3])
plot([fs1p/0.5 fs1p/0.5], [0 3],'--b','LineWidth',2)
txt = strcat('$\tau_1$ = ',num2str(fs1p),' $\rightarrow$');
text(fs1p/0.5,1.5,txt,'HorizontalAlignment','right','fontsize',16)
grid on
set(gca,'FontSize',16)
title('$\dot{q}_1$ - positive torque','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\dot{q}_1$ [$^\circ$]','fontsize',16);
saveas(gcf,'q1_est_pos','epsc')

figure(2)
ylim([-1 5])
plot([fs2p/0.1 fs2p/0.1], [0 5],'--b','LineWidth',2)
txt = strcat('$\tau_2$ = ',num2str(fs2p),' $\rightarrow$');
text(fs2p/0.1,2.5,txt,'HorizontalAlignment','right','fontsize',16)
grid on
set(gca,'FontSize',16)
title('$\dot{q}_2$ - positive torque','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\dot{q}_2$ [$^\circ$]','fontsize',16);
saveas(gcf,'q2_est_pos','epsc')

clear tau1_break tau2_break

for m  = 1:length(expe)    
    % Load the experimental results
    for j = 1:6
        file = strcat('est_',var(j,:),'n',pen(ceil(j/3),:),'_exp',expe(m),'_gan.t');
        fileID = fopen(file,'r');
        A = textscan(fileID,'%f %f');
        for k = 1:2
            for l = 1:length(A{1})
                switch (j)
                    case 1
                        q1(k,l) = A{k}(l);
                    case 2
                        q1Dot(k,l) = A{k}(l);
                    case 3
                        tau1(k,l) = A{k}(l);
                    case 4
                        q2(k,l) = A{k}(l);
                    case 5
                        q2Dot(k,l) = A{k}(l);
                    case 6
                        tau2(k,l) = A{k}(l);
                end
            end
        end
    end
    
    % Calculate the velocities using a seven-points center difference
    clear q1Dot
    clear q2Dot
    q1Dot(1,:) = q1(1,:);
    q2Dot(1,:) = q2(1,:);
    q1Dot(2,:) = numDiff(q1);    
    q2Dot(2,:) = numDiff(q2);
    
    % Creates a graphic with the results
    figure(3)
    plot(q1Dot(1,:),q1Dot(2,:),color(m),'LineWidth',2);
    hold on
    figure(4)
    plot(q2Dot(1,:),q2Dot(2,:),color(m),'LineWidth',2);
    hold on
    
    % Searchs for the first spike, using some threshold and calculates the
    % torque applied at the moment the spike appeared
    flag = 0;
    i = 1;
    while (flag==0)
        if (q1Dot(2,i) < -0.2)
            tau1_break(m) = q1Dot(1,i)*-0.5;
            flag = 1;
        else
            i = i + 1;
        end
    end
    
    flag = 0;
    i = 1;
    while (flag==0)
        if (q2Dot(2,i) < -0.2)
            tau2_break(m) = q2Dot(1,i)*-0.1;
            flag = 1;
        else
            i = i + 1;
        end
    end
    
    clear q1 q2 q1Dot q2Dot tau1 tau2
end

% Calculates the static frcition coeffcients as the mean value of the
% results
fs1n = mean(tau1_break)
fs2n = mean(tau2_break)

% Adds a visual indication of the static frcition coeffcients to the
% graphics
figure(3)
plot([fs1n/-0.5 fs1n/-0.5], [-6 0],'--b','LineWidth',2)
txt = strcat('$\tau$ = ',num2str(fs1n),' $\rightarrow$');
text(fs1n/-0.5,-5.5,txt,'HorizontalAlignment','right','fontsize',16)
ylim([-6 0]);
xlim([0 25])
grid on
set(gca,'FontSize',16)
title('$\dot{q}_1$ - negative torque','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\dot{q}_1$ [$^\circ$]','fontsize',16);
saveas(gcf,'q1_est_neg','epsc')

figure(4)
plot([fs2n/-0.1 fs2n/-0.1], [-20 0],'--b','LineWidth',2)
txt = strcat('$\tau$ = ',num2str(fs2n),' $\rightarrow$');
text(fs2n/-0.1,-10,txt,'HorizontalAlignment','right','fontsize',16)
ylim([-20 0])
xlim([0 20])
grid on
set(gca,'FontSize',16)
title('$\dot{q}_2$ - negative torque','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\dot{q}_2$ [$^\circ$]','fontsize',16);
saveas(gcf,'q2_est_neg','epsc')

% Seven-points center difference
function der = numDiff(x)
    h = 0.015;
    n = length(x(1,:));
    
    der = zeros(1,n);
    
    for i = 4:n-3
        der(i) = (-x(2,i-3) + 9*x(2,i-2) - 45*x(2,i-1) + 45*x(2,i+1) - 9*x(2,i+2) + x(2,i+3))/(60*h);
    end
    
    der(n-2:n) = x(2,n-2:n);
end