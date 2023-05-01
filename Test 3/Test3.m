%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Friccion Dinamica
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
set(0,'defaultTextInterpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sigma_0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:3
    file = strcat(strcat('sig0_pos1_6_exp',int2str(i),'_gan.t'));
    fileID = fopen(file,'r');
    A = textscan(fileID,'%f %f');
    fclose('all');

    q1 = A{2};
    t = A{1};

    figure(1)
    plot(t,q1,'LineWidth',2);
    hold on
end

grid on
set(gca,'FontSize',16)
title('$q_1$');
xlabel('time [s]','fontsize',16);
ylabel('$q_1$ [$^\circ$]','fontsize',16);
saveas(gcf,'q1_sigma0','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:3    
    file = strcat(strcat('sig0_pos2_08_exp',int2str(i),'_gan.t'));
    fileID = fopen(file,'r');
    A = textscan(fileID,'%f %f');
    fclose('all');

    q2 = A{2};
    t = A{1};

    figure(2)
    plot(t,q2,'LineWidth',2);
    hold on
end

grid on
set(gca,'FontSize',16)
title('$q_2$');
xlabel('time [s]','fontsize',16);
ylabel('$q_2$ [$^\circ$]','fontsize',16);
saveas(gcf,'q2_sigma0','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta1u = [0.0046  0.0039 -0.0004 0.0046 0.0039 -0.0007 0.0046 0.0039 -0.001 0.0046 0.0039 -0.0011];
theta1d = [-0.0109 -0.012 -0.0165 -0.0113 -0.0123 -0.0169 -0.0113 -0.0123 -0.0172 -0.0116 -0.0123 -0.0172];

Dtheta1 = theta1u - theta1d;

sigma01 = 6/mean(Dtheta1)*180/pi

theta2u = [0.0346 0.0286 -0.0159 0.0401 0.0319 -0.0137 0.0412 0.033 -0.0126 0.0428 0.0313 -0.0121 0.0417 0.033 -0.011];
theta2d = [-0.0077 -0.0187 -0.0571 -0.0077 -0.017 -0.056 -0.0071 -0.0132 -0.0549 -0.0044 -0.0137 -0.0549 -0.0038 -0.0143 -0.0549];

Dtheta2 = theta2u - theta2d;

sigma02 = 0.8/mean(Dtheta2)*180/pi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sigma_1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J1 = 2.519;
omegan1 = sqrt(sigma01/J1);
fv1 = 3.1465;
A1 = 6;

for i = 1:3  
    file = strcat('sig1_pos1_6_exp',int2str(i),'_gan.t');
    fileID = fopen(file,'r');
    A = textscan(fileID,'%f %f');
    fclose('all');

    q1 = A{2};
    t = A{1};
    
    k = 1;
    for j = 1:length(t)
        if (t(j) > 1 && t(j) < 1.5)
            t_resp(k) = t(j) - 1;
            q_resp(k) = q1(j);
            k = k + 1;
        end
    end

    figure(3)
    plot(t,q1,'LineWidth',2);
    hold on
    
    q_resp0 = q_resp(1);
    q_resp = q_resp - q_resp(1);
    qf = q_resp(end);
    K = qf/A1;

    ft = fittype(@(zeta, x) K*A1*(1 - (zeta + sqrt(zeta.^2 - 1))/sqrt(zeta.^2 - 1)*exp(-(zeta - sqrt(zeta.^2 - 1)).*x*omegan1) + (zeta - sqrt(zeta.^2 - 1))/sqrt(zeta.^2 - 1)*exp(-(zeta + sqrt(zeta.^2 - 1)).*x*omegan1)));
    options = fitoptions(ft);
    options.StartPoint = [10];
    options.Lower = [3];
    options.Upper = [20];
    [FIT,gof] = fit(t_resp',q_resp',ft,options)
    q0*(1 - exp(-FIT.zeta.*(t_resp + 1)*omegan1)).*(cosh(omegan1*sqrt(FIT.zeta^2 - 1).*(t_resp + 1)) + FIT.zeta/sqrt(FIT.zeta^2 - 1)*sinh(omegan1*sqrt(FIT.zeta^2 - 1).*(t_resp + 1))),'k*')
    a = (FIT.zeta - sqrt(FIT.zeta^2 - 1))*omegan1;
    b = (FIT.zeta + sqrt(FIT.zeta^2 - 1))*omegan1;
    q_est = K*6*(1 - b/(b - a)*exp(-a*t_resp) + a/(b - a)*exp(-b*t_resp))+q_resp0;
    plot(t_resp + 1,q_est,'k*')
    
    sigma11(i) = 2*FIT.zeta*sqrt(J1*sigma01) - fv1;
    clear q_resp t_resp
end

xlim([0 3])
grid on
set(gca,'FontSize',16)
title('$q_1$');
xlabel('time [s]','fontsize',16);
ylabel('$q_1$ [$^\circ$]','fontsize',16);
legend({'Run 1','Fitted function','Run2','Fitted function','Run3','Fitted function'},'Location','east')
saveas(gcf,'q1_sigma1','epsc')

disp('sigma11 = ')
vpa(mean(sigma11))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J2 = 0.102;
omegan2 = sqrt(sigma02/J2);
fv2 = 0.1237;
A2 = 0.8;

for i = 1:3  
    file = strcat('sig1_pos2_08_exp',int2str(i),'_gan.t');
    fileID = fopen(file,'r');
    A = textscan(fileID,'%f %f');
    fclose('all');

    q2 = A{2};
    t = A{1};
    
    k = 1;
    for j = 1:length(t)
        if (t(j) > 1 && t(j) < 1.5)
            t_resp(k) = t(j) - 1;
            q_resp(k) = q2(j);
            k = k + 1;
        end
    end

    figure(20)
    plot(t,q2,'LineWidth',2);
    hold on
    
    q_resp0 = q_resp(1);
    q_resp = q_resp - q_resp(1);
    qf = q_resp(end);
    K = qf/A2;

    ft = fittype(@(zeta, x) K*A2*(1 - (zeta + sqrt(zeta.^2 - 1))/sqrt(zeta.^2 - 1)*exp(-(zeta - sqrt(zeta.^2 - 1)).*x*omegan2) + (zeta - sqrt(zeta.^2 - 1))/sqrt(zeta.^2 - 1)*exp(-(zeta + sqrt(zeta.^2 - 1)).*x*omegan2)));
    options = fitoptions(ft);
    options.StartPoint = [10];
    options.Lower = [1];
    options.Upper = [20];
    [FIT,gof] = fit(t_resp',q_resp',ft,options)

    a = (FIT.zeta - sqrt(FIT.zeta^2 - 1))*omegan2;
    b = (FIT.zeta + sqrt(FIT.zeta^2 - 1))*omegan2;
    q_est = K*A2*(1 - b/(b - a)*exp(-a*t_resp) + a/(b - a)*exp(-b*t_resp))+q_resp0;
    plot(t_resp + 1,q_est,'k*')
    
    sigma12(i) = 2*FIT.zeta*sqrt(J2*sigma02) - fv2;
    clear q_resp t_resp
end
xlim([0 5])
grid on
set(gca,'FontSize',16)
title('$q_2$');
xlabel('time [s]','fontsize',16);
ylabel('$q_2$ [$^\circ$]','fontsize',16);
legend({'Run 1','Fitted function','Run2','Fitted function','Run3','Fitted function'},'Location','southeast')
saveas(gcf,'q2_sigma1','epsc')

disp('sigma12 = ')
vpa(mean(sigma12))