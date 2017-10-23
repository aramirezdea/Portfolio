%% MECH454 Vehicle Dynamics and Control
% Final Project Script for Professor Craig Beal
% Written by Alejandro Ramirez de Arellano in December 2016

%% Cleanup 
clc
% clear all
close all

%% Parameters 
velocity = 18;      % m/s
wheelbase = 0.7;    % meters
Ktf = 2*1.5e5 ;     % spring constant of front tire
Ktr = Ktf;          % spring constant of rear tire
mf = 1.81;          % mass of front wheel in kg
mr = 2 ;            % mass of rear wheel in kg
mcf = 0.6*(145-3.81); % mass of front half of bike
mb = 18;            % battery weight in kg
mh = 82;            % mass of human rider in kg
Ksf = 0.37*1000;    % kg/m
Csf = 3500;       % damping coeff front
Ksr = Ksf;      % spring constant rear
Csr = 0.85*Csf;      % damping coeff rear
M = 145-mf-mr;      % mass of entire bicycle in kg
Iyy = (((2/5)*mf*0.21^2)+mf*0.47^2) + (((2/5)*mr*0.21^2)+mr*0.62^2)+...
    ((mb*(0.3^2+0.2^2)/12)+mb*0.2^2)+(((mh*(1.83^2+0.2^2)/12)+mh*0.57^2));          % Moment of Inertia in SI
ha = 0.25 ;         % meters
hb = 0.5;           % meters
a = 0.3 ;           % meters
b = 0.7;            % meters
asquared = a^2;     % meters squared
bsquared = b^2;     % meters squared
hasquared = ha^2;   % meters squared
hbsquared = hb^2;   % meters squared
alpha =80*(pi/180); % radians
beta = 40*(pi/180); % radians
zroad = 0.03;       % displacement of a speed bump

%% Equations of Motion in Matrix Form
% Setting zero matrices as foundation
A = zeros(8);
B = zeros(8,2);
C = zeros(2,8);
D = zeros(2);
% Adding in the equivalence statements
A(1,2) = 1;
A(3,4) = 1;
A(5,6) = 1;
A(7,8) = 1;

% Modifying the A zeros vector row by row based on EOMs
A(2,:) = [ (-Ksf-Ktf)/mf -Csf/mf 0 0 Ksf/mf Csf/mf (-Ksf*a)/(mf*sin(alpha)) (-Csf*a)/(mf*sin(alpha))];  
A(4,:) = [ 0 0 (-Ksr-Ktr)/mr -Csr/mr Ksr/mr Csr/mr (Ksr*b)/(mr*sin(beta)) (Csr*b)/(mr*sin(beta))];
A(6,:) = [ Ksf/M Csf/M Ksr/M Csr/M (-Ksf-Ksr)/M (-Csf-Csr)/M (Ksf*a)/(M*sin(alpha))+(Ksr*b)/(M*sin(beta)) (Csf*a)/(M*sin(alpha))-(Csr*b)/(M*sin(beta))];
A(8,:) = (1/Iyy).*[ ((-Ksf*a*sin(alpha))-(-Ksf*ha*cos(alpha))) ((-Csf*a*sin(alpha))-(Csf*ha*cos(alpha))) ((Ksr*b*sin(beta))-(Ksr*hb*cos(beta)))...
    ((Csr*b*sin(beta))-Csr*hb*cos(beta)) ((-Ksr*b*sin(beta)+Ksf*a*sin(alpha)+Ksr*hb*cos(beta)+Ksf*ha*cos(alpha))) ...
    ((-Csr*b*sin(beta)+Csf*a*sin(alpha)+Csr*hb*cos(beta)+Csf*ha*cos(alpha))) (-Ksr*bsquared-Ksf*asquared-Ksr*hbsquared-Ksf*hasquared) ...
    (-Csr*bsquared-Csf*asquared-Csr*hbsquared-Csf*hasquared)];
% Modifying the B matrix based on EOMs
B(2,1) = Ktf/mf;
B(4,2) = Ktr/mr;
% Modifying the C matrix based on desired outputs
C(1,5) = 1;
C(2,7) = 1;
%Modifying D based on desired outputs
% Assigning these matrices as a state-space
SYS = ss(A,B,C,D);

%Creating input function of road surface
Tstep = 0.01;
T = 0:Tstep:100;
heightSpeedBump = 0.03;
widthSpeedBump = 0.08;
duration = (velocity/widthSpeedBump)^-1;
sizeB = duration/Tstep;
sizeB = round(sizeB,0);
partA = zeros(1,5000);
partB = heightSpeedBump*ones(1,4);
partC = zeros(1,100001-5000);
% Creating the U vectors for each road surface case
    % speed bump
    heightSpeedBump = 0.03;
    widthSpeedBump = 0.08;
    U_speedBump = zeros(10001,2);
    for i = 1:4
        U_speedBump(5000+i,1) = heightSpeedBump;
        U_speedBump(5004+i,2) = heightSpeedBump;
    end

    % Up a street curb
    heightCurb = 0.08;
    U_UpCurb = zeros(10001,2);
    for i = 1:5001
        U_UpCurb(5000+i,1) = heightCurb;
        if i < (5001-4)
            U_UpCurb(5004+i,2) = heightCurb;
        end
    end
    
        % Down a street curb
    heightCurb = -0.08;
    U_DownCurb = zeros(10001,2);
    for i = 1:5001
        U_DownCurb(5000+i,1) = heightCurb;
        if i < (5001-4)
            U_DownCurb(5004+i,2) = heightCurb;
        end
    end
        % Hitting a pothole
    heightPothole = -0.3;
    U_pothole = zeros(10001,2);
    for i = 1:2                                 %  8 time steps in the pothole 
        U_pothole(5000+i,1) = heightPothole; 
        U_pothole(5004+i,2) = heightPothole;
    end 

%% Calculations

% PART 1: Determining ideal rear shock stiffness parameters

% FRONT SUSPENSISON CONFIGURATION
% The front suspension will remain at fixed parameters
Nccoefs = [Ktf*Csf Ksf*Ktf]/1e5;
Nwcoefs = Ktf*[mcf Csf Ksf]/1e5;
Dcoefs = [(mcf*mf), (Csf*mcf + Csf*mf), (mcf*(Ksf + Ktf) + Ksf*mf),...
    (Csf*(Ksf + Ktf) - Csf*Ksf), (Ksf*(Ksf + Ktf) - Ksf^2)]/1e4;

Gf = tf(Nccoefs,Dcoefs);
Tf = 0:0.001:2+wheelbase/velocity; 
Yf = impulse(Gf*zroad,Tf); 

% REAR SUSPENSION CONFIGURATIONS
% Vary the rear (manually) to determine the ideal combination given front suspension
% specs

% Constant Parameters
mcr = 0.4*(145-3.81);   %mass in kg with 180 lb rider
mwr = 2;                %mass of rear wheel assembly in kg ~5 lbs
kwr = 2*1.5e5;          %wheel stiffness rear

%Varying Rear Damping Rate

ksr = Ksf; %Equal front and rear
bsr = Csf;
    Nccoefs = [kwr*bsr ksr*kwr]/1e5;
    Dcoefs  = [(mcr*mwr), (bsr*mcr + bsr*mwr), (mcr*(ksr + kwr) + ksr*mwr),...
        (bsr*(ksr + kwr) - bsr*ksr), (ksr*(ksr + kwr) - ksr^2)]/1e4;
    GrEqual = tf(Nccoefs,Dcoefs);
    
ksr = Ksf; % 30% softer rear
bsr = .7*Csf;
    Nccoefs = [kwr*bsr ksr*kwr]/1e5;
    Dcoefs  = [(mcr*mwr), (bsr*mcr + bsr*mwr), (mcr*(ksr + kwr) + ksr*mwr),...
        (bsr*(ksr + kwr) - bsr*ksr), (ksr*(ksr + kwr) - ksr^2)]/1e4;
    GrSofter = tf(Nccoefs,Dcoefs);
    
ksr = Ksf; % 15% softer rear
bsr = 0.85*Csf;
    Nccoefs = [kwr*bsr ksr*kwr]/1e5;
    Dcoefs  = [(mcr*mwr), (bsr*mcr + bsr*mwr), (mcr*(ksr + kwr) + ksr*mwr),...
        (bsr*(ksr + kwr) - bsr*ksr), (ksr*(ksr + kwr) - ksr^2)]/1e4;
    Gr10 = tf(Nccoefs,Dcoefs);
    
ksr = Ksf; % 15% stiffer rear
bsr = 1.15*Csf;
    Nccoefs = [kwr*bsr ksr*kwr]/1e5;
    Dcoefs  = [(mcr*mwr), (bsr*mcr + bsr*mwr), (mcr*(ksr + kwr) + ksr*mwr),...
        (bsr*(ksr + kwr) - bsr*ksr), (ksr*(ksr + kwr) - ksr^2)]/1e4;
    Gr30 = tf(Nccoefs,Dcoefs);
    
ksr = Ksf; % 30% stiffer rear
bsr = 1.3*Csf;
    Nccoefs = [kwr*bsr ksr*kwr]/1e5;
    Dcoefs  = [(mcr*mwr), (bsr*mcr + bsr*mwr), (mcr*(ksr + kwr) + ksr*mwr),...
        (bsr*(ksr + kwr) - bsr*ksr), (ksr*(ksr + kwr) - ksr^2)]/1e4;
    Gr50 = tf(Nccoefs,Dcoefs);
    
% Setting up Tr and Yr
Tr = 0:0.001:2-0.001; %What is this??
array = wheelbase/velocity/0.001;
integer = round(array,0);
YrEqual = [zeros(integer,1); impulse(GrEqual*zroad,Tr)];
YrSofter = [zeros(integer,1); impulse(GrSofter*zroad,Tr)];
Yr10 = [zeros(integer,1); impulse(Gr10*zroad,Tr)];
Yr30 = [zeros(integer,1); impulse(Gr30*zroad,Tr)]; %
Yr50 = [zeros(integer,1); impulse(Gr50*zroad,Tr)];

% PART 2. Calculating the system response to particular input
SOLN_speedBump = lsim(SYS,U_speedBump,T);
SOLN_UpCurb = lsim(SYS,U_UpCurb,T);
SOLN_DownCurb = lsim(SYS,U_DownCurb,T);
SOLN_pothole = lsim(SYS,U_pothole,T);

% Adding Bounce and Pitch to model based on SOLN=lsim(SYS,U,T)
bounce_speedBump = SOLN_speedBump(:,1);
pitch_speedBump = SOLN_speedBump(:,2);

% PART 3. Finding the Eigenvalues to ensure stability, Natural Frequency, and Damping Ratio
[Eigenvectors,Eigenvalues] = eig(A);
naturalFrequency = sqrt(real(Eigenvalues)^2+imag(Eigenvalues)^2);
RealParts = real(Eigenvalues);
dampingRatio = sin((-1./RealParts)/naturalFrequency);

%% Plots

% Plotting Vertical Suspension Motion to tune rear stiffness
figure(1)
plot(Tf,Yf,'b',Tf,YrEqual,'c-.',Tf,YrSofter,'k-.',Tf,Yr10,'m-',Tf,Yr30, 'g-',Tf,Yr50, 'r-')
legend('Front','Equal Rear','30% softer','15% softer','15% stiffer','30% stiffer');
xlim([0 0.25])
ylabel('Displacement (m)');
xlabel('Time (s)');
title('Response to Speed Bump at Range of Rear Damp Rates')
TfOLD = Tf;
YfOLD = Yf;
YrSofterOLD = YrSofter;

% % For Report Purposes Only
% figure(1)
% hold on
% plot(TfOLD,Yf,'--',Tf,YrSofterOLD,'--')
% plot(Tf,Yf,'b',Tf,YrSofter,'k-.')
% legend('Front 60/40','Rear 60/40','Front 50/50','Rear 50/50');
% hold off
% xlim([0 0.25])
% ylabel('Displacement (m)');
% xlabel('Time (s)');
% title('Response to Speed Bump at Range of Rear Damp Rates')

%Plotting response to range of urban environment inputs 
figure(2);
subplot(411)
hold on
plot(SOLN_speedBump)
plot(zeros(1,length(SOLN_speedBump)),'--k')
hold off
legend('Body Bounce','Pitch Angle');
title('Speed bump with 40° rear shock angle');
xlim([5000 5020])
ylabel('Displacement (m) or Angle (rads)');
xlabel('Time (s)');

subplot(412)
hold on
plot(SOLN_UpCurb)
plot(zeros(1,length(SOLN_speedBump)),'--k')
hold off
legend('Body Bounce','Pitch Angle');
title('Riding up a curb with 40° rear shock angle');
xlim([5000 5020])
ylabel('Displacement (m) or Angle (rads)');
xlabel('Time (s)');

subplot(413)
hold on
plot(SOLN_DownCurb)
plot(zeros(1,length(SOLN_speedBump)),'--k')
hold off
legend('Body Bounce','Pitch Angle');
title('Riding down a Curb with 40° rear shock angle');
xlim([5000 5020])
ylabel('Displacement (m) or Angle (rads)');
xlabel('Time (s)');

subplot(414)
hold on
plot(SOLN_pothole)
plot(zeros(1,length(SOLN_speedBump)),'--k')
hold off
legend('Body Bounce','Pitch Angle');
title('Hitting pothole with 40° rear shock angle');
xlim([5000 5020])
ylabel('Displacement (m) or Angle (rads)');
xlabel('Time (s)');

%S-Plane plot for demonstration of Stability, Speed of Response, and
Natural Frequency
figure(3)
plot(real(Eigenvalues),imag(Eigenvalues),'o')
title('S-Plane')
xlabel('Real Axis')
ylabel('Imaginary Axis')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

%Manually rearrange these to store variables for range of Betas
Beta80_1 = SOLN_speedBump; 
Beta80_2 = SOLN_UpCurb;
Beta80_3 = SOLN_DownCurb;
Beta80_4 = SOLN_pothole;

% % For Report Purposes Only, uses stored variables from range of beta to display findings
%on one graph
% 
% figure(4);
% hold on
% plot(Beta20_1)
% plot(Beta40_1)
% plot(Beta60_1)
% plot(Beta80_1)
% plot(zeros(1,length(SOLN_speedBump)),'--k')
% hold off
% legend('Bounce at 20° ','Pitch at 20°','Bounce at 40° ','Pitch at 40°','Bounce at 60° ','Pitch at 60°','Bounce at 80° ','Pitch at 80°');
% title('Speed Bump Response at a Range of Rear Shock Angles');
% xlim([5000 5020])
% ylabel('Displacement (m) or Angle (rads)');
% xlabel('Time (s)');
% 
% figure(5)
% hold on
% plot(Beta20_2)
% plot(Beta40_2)
% plot(Beta60_2)
% plot(Beta80_2)
% plot(zeros(1,length(SOLN_speedBump)),'--k')
% hold off
% legend('Bounce at 20° ','Pitch at 20°','Bounce at 40° ','Pitch at 40°','Bounce at 60° ','Pitch at 60°','Bounce at 80° ','Pitch at 80°','location','east');
% title('Riding up a curb at a Range of Rear Shock Angles');
% xlim([5000 5020])
% ylabel('Displacement (m) or Angle (rads)');
% xlabel('Time (s)');
% 
% figure(6)
% hold on
% plot(Beta20_3)
% plot(Beta40_3)
% plot(Beta60_3)
% plot(Beta80_3)
% plot(zeros(1,length(SOLN_speedBump)),'--k')
% hold off
% legend('Bounce at 20° ','Pitch at 20°','Bounce at 40° ','Pitch at 40°','Bounce at 60° ','Pitch at 60°','Bounce at 80° ','Pitch at 80°','location','east');
% title('Riding down a Curb at a Range of Rear Shock Angles');
% xlim([5000 5020])
% ylabel('Displacement (m) or Angle (rads)');
% xlabel('Time (s)');
% 
% figure(7)
% hold on
% plot(Beta20_4)
% plot(Beta40_4)
% plot(Beta60_4)
% plot(Beta80_4)
% plot(zeros(1,length(SOLN_speedBump)),'--k')
% hold off
% legend('Bounce at 20° ','Pitch at 20°','Bounce at 40° ','Pitch at 40°','Bounce at 60° ','Pitch at 60°','Bounce at 80° ','Pitch at 80°');
% title('Hitting pothole at a Range of Rear Shock Angles');
% xlim([5000 5020])
% ylabel('Displacement (m) or Angle (rads)');
% xlabel('Time (s)');
% END OF SCRIPT