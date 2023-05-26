%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Test 2 analysis. Viscous, Coulomb and Stribeck effect parameters estimation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear q1 q2 q1Dot q2Dot q1DotDot q2DotDot tau1 tau2 tauf1 tauf2 tauf1r tauf2r;
close all;
clc;
set(0,'defaultTextInterpreter','latex');

% fs1p = 10.75;
% fs1n = -12.00;
% fs2p = 1.49;
% fs2n = -1.45;

% These string arrays are used to automate the load of the experimental
% data
veldc = char('01','02','03','04','05','1','2','3','4','5','10','15','20','30','40','50','60','70','80','90');
veld = [0.1 0.2 0.3 0.4 0.5 1 2 3 4 5 10 15 20 30 40 50 60 70 80 90];
var = char('pos1_','vel1_','tau1_','pos2_','vel2_','tau2_');

% These array indicates the start and end time of the steady states in each 
% run of test 2 with positive speeds. 
t_ini = [19.0 6.0 3.3 2.6 1.5 1.0 1.0 0.8 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.3];
t_fin = [38.0 38.0 38.0 38.0 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.4 3.5 2.9 2.5 2.2 1.9];

g=9.81;

% This value is used to index the figures
fig = 0;

% Load the experimental data for positive speed of both q1 and q2,
% calculates f_i and the average of qDot and f_i in the steady state for
% each run of test 2.
for i = 1:length(veldc)
    % Load the experimental data
    for j = 1:6
        file = strcat('cin_',var(j,:),veldc(i,:),'_exp_gan.t');
        fileID = fopen(file,'r');
        A = textscan(fileID,'%f %f');
        
        for k = 1:2
            for l = 1:length(A{1})
                exp_res(i,j,k,l) = A{k}(l);
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
    
    %calculates qDot using a seven-points center difference    
    q1Dot(2,:) = numDiff(q1);
    q2Dot(2,:) = numDiff(q2);
    
    % calculates f_i
    n = 0;
    Sq1p = 0;
    Sq2p = 0;
    for m = 1:length(q1Dot)
        if (q1Dot(1,m) >= t_ini(i) && q1Dot(1,m) <= t_fin(i))
            n = n + 1;
            Sq1p = q1Dot(2,m) + Sq1p;
            Sq2p = q2Dot(2,m) + Sq2p;
            m11 = 2.351 + 0.168*cos(q2(2,m)*pi/180);
            m12 = 0.102 + 0.084*cos(q2(2,m)*pi/180);
            m21 = m12;
            m22 = 0.102;
            q1e = q1(1,m)*veld(i) - q1(2,m);
            q2e = q2(1,m)*veld(i) - q2(2,m);
            q1pe = veld(i) - q1Dot(2,m);
            q2pe = veld(i) - q2Dot(2,m);
            c11=-0.168*sin(q2(2,m)*pi/180)*q2Dot(2,m)*pi/180;
            c12=-0.084*sin(q2(2,m)*pi/180)*q2Dot(2,m)*pi/180;
            c21=0.084*sin(q2(2,m)*pi/180)*q1Dot(2,m)*pi/180;
            c22=0.0;
            g1=g*(3.921*sin(q2(2,m)*pi/180)+0.186*sin((q1(2,m) + q2(2,m))*pi/180));
            g2=g*0.186*sin((q1(2,m) + q2(2,m))*pi/180);
            tauf1r(n) = tau1(2,m) - c11*q1Dot(2,m)*pi/180 - c12*q2Dot(2,m)*pi/180 - g1;
            tauf2r(n) = tau2(2,m) - c21*q1Dot(2,m)*pi/180 - c22*q2Dot(2,m)*pi/180 - g2;
            
            % stores the f_i for a desired speed of 40 degrees per second
            % in q1 to show it as an example of the results in a run of
            % test 2
            if i == 15
                ff(n) = tauf1r(n);
                tf(n) = q1Dot(1,m);
            end
        end
    end
    
    % Average f_i and qDot
    q1Dot_ave(i) = Sq1p/n;
    q2Dot_ave(i) = Sq2p/n;
    tauf1_ave(i) = median(tauf1r);
    tauf2_ave(i) = median(tauf2r);
    hold off
    clear q1 q2 q1Dot q2Dot q1DotDot q2DotDot tau1 tau2 tauf1 tauf2 tauf1r tauf2r;
end

fclose('all');

syms x 'real'

% Estimates the viscous and Coulomb frcition parameters using polyfit with
% a polynomial of first degree
[p1p,S1] = polyfit(q1Dot_ave(11:length(veldc)),tauf1_ave(11:length(veldc)),1);
[p2p,S2] = polyfit(q2Dot_ave(11:length(veldc)),tauf2_ave(11:length(veldc)),1);
disp('Positives')
disp('fv1 = ')
disp(p1p(1)*180/pi)
disp('fv2 = ')
disp(p2p(1)*180/pi)
disp('fc1 = ')
disp(p1p(2))
disp('fc2 = ')
disp(p2p(2))

% evaluates the first order polynomial computed above and calculates the
% RMS value of the error between the polynomial and the experients data
y_fit = polyval(p1p,q1Dot_ave(11:length(veldc)),S1)
RMSE = sqrt(1/length(y_fit)*dot((tauf1_ave(11:length(veldc)) - y_fit),(tauf1_ave(11:length(veldc)) - y_fit)))
y_fit = polyval(p2p,q2Dot_ave(11:length(veldc)),S2)
RMSE = sqrt(1/length(y_fit)*dot((tauf2_ave(11:length(veldc)) - y_fit),(tauf2_ave(11:length(veldc)) - y_fit)))

% Creates a custom model for the fitting of the experimental f_k_i curve
% into the stribeck effect curve
ft = fittype(@(vs,ds,x) p1p(2) + (fs1p - p1p(2))*exp(-(x/vs).^ds) + p1p(1)*x);
options = fitoptions(ft);
options.StartPoint = [6, 0.5];
options.Lower =      [0.1, 0];
options.Upper =      [10, 2];
[p1,gof] = fit(q1Dot_ave(5:10)',tauf1_ave(5:10)',ft,options)

ft = fittype(@(vs,ds,x) p2p(2) + (fs2p - p2p(2))*exp(-(x/vs).^ds) + p2p(1)*x);
options = fitoptions(ft);
options.StartPoint = [6, 0.5];
options.Lower =      [2, 0];
options.Upper =      [10, 10];
[p2,gof] = fit(q2Dot_ave(5:10)',tauf2_ave(5:10)',ft,options)

% Generates curves using the parameters estimations
x = linspace(0,90,300);
y1 = p1p(2) + (fs1p - p1p(2))*exp(-(x/p1.vs).^p1.ds) + p1p(1)*x;
y2 = p2p(2) + (fs2p - p2p(2))*exp(-(x/p2.vs).^p2.ds) + p2p(1)*x;

% Creates the graphics that depicts the esprimentally obtained f_k_i curve
% and the one calculated with the estimated parameters
figure(1)
plot(q1Dot_ave,tauf1_ave,'*','LineWidth',2);
hold on
plot(x,y1,'k','LineWidth',2);
set(gca,'FontSize',16)

figure(2)
plot(q2Dot_ave,tauf2_ave,'*','LineWidth',2);
hold on
plot(x,y2,'k','LineWidth',2);
set(gca,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the experimental data for positive speed of both q1 and q2,
% calculates f_i and the average of qDot and f_i in the steady state for
% each run of test 2.

% These array indicates the start and end time of the steady states in each 
% run of test 2 with negative speeds. 
t_ini = [8.0 4.5 3.0 2.0 1.5 1.0 0.5 0.5 0.4 0.4 0.3 0.3 0.3 0.5 0.5 0.5 0.5 0.5 0.5 0.6];
t_fin = [38.0 38.0 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.3 3.4 2.8 2.5 2.1 1.9];

g=9.81;

% Load the experimental data for negative speed of both q1 and q2,
% calculates f_i and the average of qDot and f_i in the steady state for
% each run of test 2.
for i = 1:length(veldc)
    % Load the experimental data
    for j = 1:6
        file = strcat('cin_',var(j,:),'n',veldc(i,:),'_exp_gan.t');
        fileID = fopen(file,'r');
        A = textscan(fileID,'%f %f');
        
        for k = 1:2
            for l = 1:length(A{1})
                exp_res(i,j,k,l) = A{k}(l);
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
    
    % Calculates qDot using a seven-points center difference   
    q1Dot(2,:) = numDiff(q1);
    q2Dot(2,:) = numDiff(q2);
        
    % calculates f_i
    n = 0;
    Sq1p = 0;
    Sq2p = 0;
    
    for m = 1:length(q1Dot)
        if (q1Dot(1,m) >= t_ini(i) && q1Dot(1,m) <= t_fin(i))
            n = n + 1;
            Sq1p = q1Dot(2,m) + Sq1p;
            Sq2p = q2Dot(2,m) + Sq2p;
            m11 = 2.351 + 0.168*cos(q2(2,m)*pi/180);
            m12 = 0.102 + 0.084*cos(q2(2,m)*pi/180);
            m21 = m12;
            m22 = 0.102;
            q1e = -q1(1,m)*veld(i) - q1(2,m);
            q2e = -q2(1,m)*veld(i) - q2(2,m);
            q1pe = -veld(i) - q1Dot(2,m);
            q2pe = -veld(i) - q2Dot(2,m);
            c11=-0.168*sin(q2(2,m)*pi/180)*q2Dot(2,m)*pi/180;
            c12=-0.084*sin(q2(2,m)*pi/180)*q2Dot(2,m)*pi/180;
            c21=0.084*sin(q2(2,m)*pi/180)*q1Dot(2,m)*pi/180;
            c22=0.0;
            g1=g*(3.921*sin(q2(2,m)*pi/180)+0.186*sin((q1(2,m) + q2(2,m))*pi/180));
            g2=g*0.186*sin((q1(2,m) + q2(2,m))*pi/180);
            tauf1r(n) = tau1(2,m) - c11*q1Dot(2,m)*pi/180 - c12*q2Dot(2,m)*pi/180 - g1;
            tauf2r(n) = tau2(2,m) - c21*q1Dot(2,m)*pi/180 - c22*q2Dot(2,m)*pi/180 - g2;
        end
    end
    
    % Average f_i and qDot
    q1Dot_ave(i) = Sq1p/n;
    q2Dot_ave(i) = Sq2p/n;
    tauf1_ave(i) = median(tauf1r);
    tauf2_ave(i) = median(tauf2r);
    
    clear clear q1 q2 q1Dot q2Dot q1DotDot q2DotDot tau1 tau2 tauf1 tauf2 tauf1r tauf2r;
end

fclose('all');

syms x 'real'

% Estimates the viscous and Coulomb frcition parameters using polyfit with
% a polynomial of first degree
[p1n,S1] = polyfit(q1Dot_ave(11:length(veldc)),tauf1_ave(11:length(veldc)),1);
[p2n,S2] = polyfit(q2Dot_ave(11:length(veldc)),tauf2_ave(11:length(veldc)),1);
disp('Negatives')
disp('fv1 = ')
disp(p1n(1)*180/pi)
disp('fv2 = ')
disp(p2n(1)*180/pi)
disp('fc1 = ')
disp(p1n(2))
disp('fc2 = ')
disp(p2n(2))

% evaluates the first order polynomial computed above and calculates the
% RMS value of the error between the polynomial and the experients data
y_fit = polyval(p1n,q1Dot_ave(11:length(veldc)),S1)
RMSE = sqrt(1/length(y_fit)*dot((tauf1_ave(11:length(veldc)) - y_fit),(tauf1_ave(11:length(veldc)) - y_fit)))
y_fit = polyval(p2n,q2Dot_ave(11:length(veldc)),S2)
RMSE = sqrt(1/length(y_fit)*dot((tauf2_ave(11:length(veldc)) - y_fit),(tauf2_ave(11:length(veldc)) - y_fit)))
ft = fittype(@(vs,ds,x) p1n(2) + (fs1n - p1n(2))*exp(-(-x/vs).^ds) + p1n(1)*x);
options = fitoptions(ft);
options.StartPoint = [6, 0.5];
options.Lower =      [0.1, 0.1];
options.Upper =      [10, 10];
[p1,gof] = fit(q1Dot_ave(1:10)',tauf1_ave(1:10)',ft,options)

% Creates a custom model for the fitting of the experimental f_k_i curve
% into the stribeck effect curve
ft = fittype(@(vs,ds,x) p2n(2) + (fs2n - p2n(2))*exp(-(-x/vs).^ds) + p2n(1)*x);
options = fitoptions(ft);
options.StartPoint = [6, 0.5];
options.Lower =      [2, 0.1];
options.Upper =      [10, 10];
[p2,gof] = fit(q2Dot_ave(5:10)',tauf2_ave(5:10)',ft,options)

% Generates curves using the parameters estimations
x = linspace(-90,0,300);
y1 = p1n(2) + (fs1n - p1n(2))*exp(-(-x/p1.vs).^p1.ds) + p1n(1)*x;
y2 = p2n(2) + (fs2n - p2n(2))*exp(-(-x/p2.vs).^p2.ds) + p2n(1)*x;

% Creates the graphics that depicts the esprimentally obtained f_k_i curve
% and the one calculated with the estimated parameters
figure(1)
hold on
plot(q1Dot_ave,tauf1_ave,'*','LineWidth',2);
plot(x,y1,'k','LineWidth',2);
xlim([-100 100]);
ylim([-20 20]);
grid on
set(gca,'FontSize',16)
title('$\dot{q}_1$ vs $f_1$','fontsize',16,'interpreter','latex');
xlabel('$\dot{q}_1$ [$^\circ$/s]','fontsize',16,'interpreter','latex');
ylabel('$f_{1}$ [Nm]','fontsize',14,'interpreter','latex');
saveas(gcf,'q1_fk','epsc')

figure(2)
plot(q2Dot_ave,tauf2_ave,'*','LineWidth',2);
hold on
plot(x,y2,'k','LineWidth',2);
xlim([-100 100]);
ylim([-3 3]);
grid on
set(gca,'FontSize',16)
title('$\dot{q}_2$ vs $f_2$','fontsize',16,'interpreter','latex');
xlabel('$\dot{q}_2$ [$^\circ$/s]','fontsize',16);
ylabel('$f_{2}$ [Nm]','fontsize',16);
saveas(gcf,'q2_fk','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a graphic depicting one example of the data obtained from one run 
% of test 2, specifically the velocity q1 with desired speed of 40 degrees 
% per second, and another graphic with the torque tau_1 applied to the
% joint and the estimated f_1.

clear q1 q2 q1p q2p tau1 tau2

% Load the experimental data
file = strcat('cin_vel1_40_exp_gan.t');
fileID = fopen(file,'r');
A = textscan(fileID,'%f %f');
for k = 1:2
    for l = 1:length(A{1})
        q1Dot(k,l) = A{k}(l);
    end
end
file = strcat('cin_tau1_40_exp_gan.t');
fileID = fopen(file,'r');
A = textscan(fileID,'%f %f');
for k = 1:2
    for l = 1:length(A{1})
        tau1(k,l) = A{k}(l);
    end
end

% Creates a graphic of qDot for q1
figure(3)
plot(q1Dot(1,:),q1Dot(2,:),'LineWidth',2);
set(gca,'FontSize',16)
grid on
title('$\dot{q}_1$','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\dot{q}_1$ [$^\circ$/s]','fontsize',16);
saveas(gcf,'q1p_test2_example','epsc')

% Creates a graphic of tau_1 and f_1 for q1
figure(4)
plot(tau1(1,:),tau1(2,:),'LineWidth',2);
hold on
plot(tf,ff,'LineWidth',2)
set(gca,'FontSize',16)
grid on
title('$\tau_1$ and $f_1$','fontsize',16);
xlabel('time [s]','fontsize',16);
ylabel('$\tau_1$, $f_1$ [Nm]','fontsize',16);
legend({'$\tau_1$','$f_1$'},'Location','southwest','fontsize',16,'interpreter','latex')
saveas(gcf,'tau1_test2_example','epsc')

% seven-points center difference
function der = numDiff(x)
    h = 0.015;
    n = length(x(1,:));
    
    der = zeros(1,n);
    
    for i = 4:n-3
        der(i) = (-x(2,i-3) + 9*x(2,i-2) - 45*x(2,i-1) + 45*x(2,i+1) - 9*x(2,i+2) + x(2,i+3))/(60*h);
    end
    
    der(n-2:n) = x(2,n-2:n);
end