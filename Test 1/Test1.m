%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Friccion estatica
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
set(0,'defaultTextInterpreter','latex');

var = char('pos1_','vel1_','tau1_','pos2_','vel2_','tau2_');
color = char('b','r','k','g');
expe = char('1','2','3','4');
pen = char('05','01');

for m  = 1:length(expe)    
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
                        q1p(k,l) = A{k}(l);
                    case 3
                        tau1(k,l) = A{k}(l);
                    case 4
                        q2(k,l) = A{k}(l);
                    case 5
                        q2p(k,l) = A{k}(l);
                    case 6
                        tau2(k,l) = A{k}(l);
                end
            end
        end
    end
    
    clear q1p
    clear q2p
    q1p(1,:) = q1(1,:);
    q2p(1,:) = q2(1,:);
    q1p(2,:) = numDiff(q1);    
    q2p(2,:) = numDiff(q2);
    
    figure(1)
    plot(q1p(1,:),q1p(2,:),color(m),'LineWidth',2);
    hold on
    figure(2)
    plot(q2p(1,:),q2p(2,:),color(m),'LineWidth',2);
    hold on
    
    clear q1 q2 q1p q2p tau1 tau2    
end

figure(1)
xlim([0 25])
ylim([-1 3])
plot([10.75/0.5 10.75/0.5], [0 3],'--b','LineWidth',2)
txt = '$\tau_1$ = 10.75 $\rightarrow$';
text(10.75/0.5,1.5,txt,'HorizontalAlignment','right','fontsize',16)
grid on
set(gca,'FontSize',16)
title('$\dot{q}_1$ - par positivo','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\dot{q}_1$ [$^\circ$]','fontsize',16);
saveas(gcf,'q1_est_pos','epsc')

figure(2)
ylim([-1 5])
plot([1.49/0.1 1.49/0.1], [0 5],'--b','LineWidth',2)
txt = '$\tau$ = 1.49 $\rightarrow$';
text(1.49/0.1,2.5,txt,'HorizontalAlignment','right','fontsize',16)
grid on
set(gca,'FontSize',16)
title('$\dot{q}_2$ - par positivo','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\dot{q}_2$ [$^\circ$]','fontsize',16);
saveas(gcf,'q2_est_pos','epsc')

for m  = 1:length(expe)    
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
                        q1p(k,l) = A{k}(l);
                    case 3
                        tau1(k,l) = A{k}(l);
                    case 4
                        q2(k,l) = A{k}(l);
                    case 5
                        q2p(k,l) = A{k}(l);
                    case 6
                        tau2(k,l) = A{k}(l);
                end
            end
        end
    end
    
    clear q1p
    clear q2p
    q1p(1,:) = q1(1,:);
    q2p(1,:) = q2(1,:);
    q1p(2,:) = numDiff(q1);    
    q2p(2,:) = numDiff(q2);
    
    figure(3)
    plot(q1p(1,:),q1p(2,:),color(m),'LineWidth',2);
    hold on
    figure(4)
    plot(q2p(1,:),q2p(2,:),color(m),'LineWidth',2);
    hold on
    
    clear q1 q2 q1p q2p tau1 tau2
end
% 
figure(3)
plot([12.0/0.5 12.0/0.5], [-6 0],'--b','LineWidth',2)
txt = '$\tau$ = -12.0 $\rightarrow$';
text(12.0/0.5,-5.5,txt,'HorizontalAlignment','right','fontsize',16)
ylim([-6 0]);
xlim([0 25])
grid on
set(gca,'FontSize',16)
title('$\dot{q}_1$ - par negativo','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\dot{q}_1$ [$^\circ$]','fontsize',16);
saveas(gcf,'q1_est_neg','epsc')

figure(4)
plot([1.45/0.1 1.45/0.1], [-20 0],'--b','LineWidth',2)
txt = '$\tau$ = -1.45 $\rightarrow$';
text(1.45/0.1,-10,txt,'HorizontalAlignment','right','fontsize',16)
ylim([-20 0])
xlim([0 20])
grid on
set(gca,'FontSize',16)
title('$\dot{q}_2$ - par negativo','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\dot{q}_2$ [$^\circ$]','fontsize',16);
saveas(gcf,'q2_est_neg','epsc')

function der = numDiff(x)
    h = 0.015;
    n = length(x(1,:));
    
    der = zeros(1,n);
    
    for i = 4:n-3
        der(i) = (-x(2,i-3) + 9*x(2,i-2) - 45*x(2,i-1) + 45*x(2,i+1) - 9*x(2,i+2) + x(2,i+3))/(60*h);
    end
    
    der(n-2:n) = x(2,n-2:n);
end