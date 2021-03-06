% Far Field

%% Variables
alpha = 2000; %m/s P wave velocity
beta = 1000; %m/s S wave velocity
rho = 2500; %kg/m3 density
fs = 25; %Hz source frequency
w = 2*pi*fs; %angular frequency
T = 1/fs; %time period

% Time and Angle Vector
t = 0:0.001:1.5; %time vector go to maximum R/beta of 1.5 s
theta = 0:10:90; %angle in 10 deg increments

% Receiver Distances (change R value if nf or ff)
%R_nf = 80; %m for near field
%R_ff = 1000; %m for far field
R = 1000; 

% Functions are defined for source wavelets
% x_o is for Far Field
% x_o_nf is for Near Field

%% Initialize Vectors
% Far Field
u1_ffs = zeros(length(theta),length(t));
u1_ffp = zeros(length(theta),length(t));
u3_ffs = zeros(length(theta),length(t));
u3_ffp = zeros(length(theta),length(t));
% Near Field
u1_nf = zeros(length(theta),length(t));
u3_nf = zeros(length(theta),length(t));
% Total Displacement
u1 = zeros(length(theta),length(t));
u3 = zeros(length(theta),length(t));
% P/S Ratio
u1_ps = zeros(1,length(theta));
u3_ps = zeros(1,length(theta));

%% Far Field Displacement Vectors
% u1 Far Field S
for j=1:length(theta)
    u1_ffs(j,:)=(1-(cosd(theta(j))^2))*(1/R)*x_o(R,beta,w,T,t)*(1/(4*pi*rho*(beta^2)));
end

% u1 Far Field P
for i=1:10
    u1_ffp(i,:)=(cosd(theta(i))^2)*(1/R)*x_o(R,alpha,w,T,t)*(1/(4*pi*rho*alpha^2));
end

figure;
subplot(2,1,1);
for i=1:length(theta)
    hold on;
    plot(t,u1_ffs(i,:));
end
title('u1 Direction Far Field S for S-R distance of 1000 m')
xlabel('Time(s)');
ylabel('Amplitude(m)');
legend('0 \circ', '10 \circ','20 \circ','30 \circ','40 \circ','50 \circ','60 \circ','70 \circ','80 \circ', '90 \circ');
grid on;
xlim([0,1.2])
subplot(2,1,2);
for i=1:length(theta)
    hold on;
    plot(t,u1_ffp(i,:));
end
title('u1 Direction Far Field P for S-R distance of 1000 m')
xlabel('Time(s)');
ylabel('Amplitude(m)');
legend('0 \circ', '10 \circ','20 \circ','30 \circ','40 \circ','50 \circ','60 \circ','70 \circ','80 \circ', '90 \circ');
grid on;
xlim([0,1.2])

% u3 Far Field S
for i=1:length(theta)
    u3_ffs(i,:)=(-cosd(theta(i))*sind(theta(i)))*(1/R)*x_o(R,beta,w,T,t)*(1/(4*pi*rho*beta^2));
end

% u3 Far Field P
for i= 1:10
   u3_ffp(i,:)=(cosd(theta(i))*sind(theta(i)))*(1/R)*x_o(R,alpha,w,T,t)*(1/(4*pi*rho*alpha^2));
end

figure;
subplot(2,1,1);
for i=1:length(theta)
    hold on;
    plot(t,u3_ffs(i,:));
end
title('u3 Direction Far Field S for S-R distance of 1000 m')
xlabel('Time(s)');
ylabel('Amplitude(m)');
legend('0 \circ', '10 \circ','20 \circ','30 \circ','40 \circ','50 \circ','60 \circ','70 \circ','80 \circ', '90 \circ');
grid on;
xlim([0,1.2])
subplot(2,1,2);
for i=1:length(theta)
    hold on;
    plot(t,u3_ffp(i,:));
end
title('u3 Direction Far Field P for S-R distance of 1000 m')
xlabel('Time(s)');
ylabel('Amplitude(m)');
legend('0 \circ', '10 \circ','20 \circ','30 \circ','40 \circ','50 \circ','60 \circ','70 \circ','80 \circ', '90 \circ');
grid on;
xlim([0,1.2])


%% Near Field Displacement Vectors
% u1 Near Field (increment over theta)
for i=1:10 
    u1_nf(i,:)=((3*(cosd(theta(i)))^2)-1)*((1/R^3)* x_o_nf(R,alpha,beta,w,T,t))*(1/(4*pi*rho));
end    

figure;
plot(t,u1_nf(1,:),t,u1_nf(10,:))
title('u1 Direction Near Field for S-R distance of 1000 m');
legend('Parallel to Force (0 \circ)', 'Perpendicular to Force (90 \circ)');
xlabel('Time(s)');
ylabel('Amplitude(m)');
grid on;


% u3 Near Field (increment over theta)
for i=1:10 
   u3_nf(i,:)=3*cosd(theta(i))*sind(theta(i))*1/(R^3)*x_o_nf(R,alpha,beta,w,T,t)*1/(4*pi*rho); 
end

figure;
plot(t,u1_nf(1,:),t,u1_nf(10,:))
title('u3 Direction Near Field for S-R distance of 80 m');
legend('0 \circ', '90 \circ');
xlabel('Time(s)');
ylabel('Amplitude(m)');

%% Total Displacements
% u1 Total Displacement
u1 = u1_ffs + u1_ffp + u1_nf;

figure;
title('u1 Displacement')
for i=1:length(theta) % increment over theta
    subplot(2,5,i);
    plot(t,u1(i,:));
    title(['angle \theta =' num2str(theta(i)) ' \circ']);
    xlabel('Time (s)');
    ylabel('Displacement (m)');
    xlim([0,1.5]);
    ylim([-3.5e-14,3.5e-14]);
end

% u3 Total Displacement
u3 = u3_ffs + u3_ffp + u3_nf;

figure;
title('u3 Displacement');
for i=1:length(theta) % increment over theta
    subplot(2,5,i);
    plot(t,u3(i,:));
    title(['angle \theta =' num2str(theta(i)) ' \circ']);
    xlabel('Time (s)');
    ylabel('Displacement (m)');
    xlim([0,1.5]);
    ylim([-3.5e-14,3.5e-14]);
end

%% Particle Motions Hodograms
figure;
title('Particle Motions Hodograms');
for i=1:length(theta) % increment over theta
    subplot(2,5,i)
    plot(u1(i,:),u3(i,:));
    title(['angle \theta =' num2str(theta(i)) ' \circ']);
    xlabel('u1 Displacement (m)');
    ylabel('u3 Displacement (m)');
    xlim([-3e-14,3e-14]);
    ylim([-2e-14,2e-14]);
end

%% P/S Ratio

for i=1:length(theta) % increment over theta
    u1_ps(i)=max(u1_ffp(i,:))/max(u1_ffs(i,:));
    u3_ps(i)=max(u3_ffp(i,:))/max(u3_ffs(i,:));
end

figure;
subplot(2,1,1);
plot(theta, u1_ps);
title('P/S Amplitude Ratios for u1 Direction')
xlabel('Angle (degrees)');
ylabel('Amplitude Ratio');
grid on;

subplot(2,1,2);
plot(theta, u3_ps);
hold on;
title('P/S Amplitude Ratios for u3 Direction')
xlabel('Angle (degrees)');
ylabel('Amplitude Ratio');
grid on;