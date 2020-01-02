% Cassini launched from earth on Oct 15, 1997 on a hyperbolic exit
% trajectory headed for Venus. At Venus, the spacecraft performed it's
% first GA with a flyby on April 26, 1998 at an altitude of 284km and a 
% speed of 11.8km/sec. 
%
%
%
% Non-Hohmann Interplanetary Trajectories: 
% Interplanetary missions are broken into three sections - (1) Departure 
% (2) Cruise (3) Arrival. We start with the cruise phase. The frame of
% reference that we use is the heliocentric ecliptic frame. The first step
% is to obtain the state vector of planet 1 at departure time (t) and the
% state vector of planet 2 at arrival time (t+t12). This is accomplished by 
% means of algorithm 8.1.
% The next step is to determine the spacecraft's transfer trajectory from
% planet 1 to planet 2. 

close all;clear all;clc
%% Launch
% Hyperbolic Exit Trajectory
muEarth = 398600;
C3 = 16.6;
VinfLaunch = sqrt(C3);
aLaunch = muEarth/C3;

% Planetary Rendevouz Venus1 (Curtis 8.8)
muVenus = 324900;
V1Alt = 284;
V1Speed = 11.8;

% Sun
muSun = 1.32712e+11;

%% Use Julian_Day_Function to confirm radii for Earth/Venus on launch day
%PLANETS:
%   Earth = 1
%   Mars = 2
%   Venus = 3
y = 1997;
UT = 12;
planet = 1;

% Call the state vector function to obtain r1 (Earth's distance from sun on launch day)
[Jan,Feb,March,April,May,June,July,Aug,Sep,Oct,Nov,Dec,Theta,listh,lista,liste,Omega,omega,I] = Cassini_State_Vector_Function(y,UT,planet);
r1 = Oct(15)
hEarthLaunch = listh(15,10);
eEarthLaunch = liste(15,10);

% Obtain the theta's of each day in october
ThetaOct = Theta(1:end,10);
ThetaLaunchEarth = ThetaOct(15);

% Generate the state vector for R1
rbar1 = (hEarthLaunch^2/muSun)*(1/(1+eEarthLaunch*cosd(ThetaLaunchEarth)))
SinCos = [cosd(ThetaLaunchEarth) sind(ThetaLaunchEarth) 0];
R1 = rbar1*SinCos;

% Generate the state vector for V1


% Send perifocal vector through DCM to get heliocentric state vector
% QxbarX = [-sin(Omega)*cos(I)*sin(omega)+cos(Omega)*cos(omega) -sin(Omega)*cos(I)*cos(omega)-cos(Omega)*sin(omega) sin(Omega)*sin(I);
%     cos(Omega)*cos(I)*sin(omega)+sin(Omega)*cos(omega) cos(Omega)*cos(I)*cos(omega)-sin(Omega)*sin(omega) -cos(Omega)*sin(I);
%     sin(I)*sin(omega) sin(I)*cos(omega) cos(I)];
% QXxbar = [-sin(Omega)*cos(I)*sin(omega)+cos(Omega)*cos(omega) cos(Omega)*cos(I)*sin(omega)+sin(Omega)*cos(omega) sin(I)*sin(omega);
%     -sin(Omega)*cos(I)*cos(omega)-cos(Omega)*sin(omega) cos(Omega)*cos(I)*cos(omega)-sin(Omega)*sin(omega) sin(I)*cos(omega);
%     sin(Omega)*sin(I) -cos(Omega)*sin(I) cos(I)];
% R1 = R1*QxbarX







% Re-define year and planet for Venus flyby
y = 1998;
planet = 3;

% Call the state vector function to obtain r2 (Venus' distance from sun at
% arrival time)
[Jan,Feb,March,April,May,June,July,Aug,Sep,Oct,Nov,Dec,Theta] = Cassini_State_Vector_Function(y,UT,planet);
r2 = April(26)

% Obtain the theta's for each day in April
ThetaApril = Theta(1:end,4);
ThetaArriveVenus = ThetaApril(26);

% Generate the state vector for R2
%R2 = r2*SinCos

%%
% Now we have the two state vectors for Earth's position, relative to the
% sun on, launch day and Venus' position, also relative to the sun, on 
% arrival day.
% Now, we need to calculate the spacecraft's planetary trajectory. This is
% done by means of Algorithm 5.2 (pg.253):
% 1. Using R1, R2, and delta T, calculate r1 and r2 using Eq.5.24
% 2. Choose either a prograde or retrograde trajectory and calculate delta
% Theta using Eq.5.26
% 3. Calculate A in Eq.5.35
% 4. By iteration, using Eqns 5.40, 5.43, and 5.45, solve Eq.5.39 for z.
% The sign of z tells us whether the orbit is a hyperbola(z<0),
% parabola(z=0), or ellipse(z>0).
% 5. Calculate y using Eq.5.38
% 6. Calculate the Lagrange f,g, and gdot functions using Eq.5.46
% 7. Calculate V1 and V2 from Eqns 5.28 and 5.29
% 8. Use R1 and V1 (or R2 and V2) in Algorithm 4.2 to obtain the orbital
% elements. 

% Total time in seconds from Earth launch to Venus flyby 1
deltaT = 16675200;

% Choosing PROGRADE trajectory, obtain delta Theta
% quadrantCheck = cross(R1,R2)
% if quadrantCheck(3) >= 0
%     dTheta = acos(dot(R1,R2)/(r1*r2));
% end
% if quadrantCheck(3) < 0
%     dTheta = 2*pi - acos(dot(R1,R2)/(r1*r2));
% end
% 
% % Calculate A from Eq.5.35 (pg.250)
% A = sin(dTheta)*sqrt((r1*r2)/(1-cos(dTheta)))
    
    
    
    





