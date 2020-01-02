%% FIRST HYPERBOLIC FLYBY OF VENUS
clear all

% Define the Orbit: Values obtained from Backup_Cassini_Phase_1.m
mu = 324900;
e = 1.70909;
delta = 71.6212; % Degrees
h = 74678.2;
Vp = 11.7863;

% Plot the orbit
figure(2)
theta = linspace(0,2*pi,100);
r = (h^2/mu)*(1./(1+e.*cos(theta)));
[theta,r] = pol2cart(theta,r); 
plot(0,0,'ko','MarkerSize',10,'MarkerFaceColor','k'), hold on, grid on
plot(theta,r,'k')
set(gca, 'fontsize', 14);xlim([-2e+04 0.2e+05])
% title('Venus Flyby 1');legend('Venus','Flyby Trajectory','Location','Best')
xlabel('Distance (km)');ylabel('Distance (km)')
% print -djpeg 'Venus_flyby_1_trajectory.jpg'

%% VENUS FLYBY 2
clear all

% Define the Orbit: Values obtained from Backup_Cassini_Phase_1.m
mu = 324900;
e = 2.72681;
delta = 43.028; % Degrees
h = 87589.3;
Vp = 13.8241;

% Plot the orbit
theta = linspace(0,2*pi,100);
r = (h^2/mu)*(1./(1+e.*cos(theta)));
[theta,r] = pol2cart(theta,r); 
plot(-319,0,'mo','MarkerSize',10,'MarkerFaceColor','m'), hold on, grid on
plot(theta,r,'m')
set(gca, 'fontsize', 14);xlim([-2e+04 0.2e+05])
title('Venus Flyby 1 & 2');
legend('Venus at 1st Flyby','Flyby 1 Trajectory','Venus at 2nd Flyby','Flyby 2 Trajectory','Location','Best')
xlabel('Distance (km)');ylabel('Distance (km)')
print -djpeg 'Both_Venus_flyby_trajectorys.jpg'