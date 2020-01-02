%% Orbital parameters for Earth on day of launch
close all;clear all;clc
% Launch information from Cassini Mission Report
%
%
% Hyperbolic Exit Trajectory
muEarth = 398600;
C3 = 16.6;
VinfLaunch = sqrt(C3); %4.0743 km/s
aLaunch = muEarth/C3; %24012.048 km

% Planetary Rendevouz Venus1 (Curtis 8.8)
muVenus = 324900;
V1Alt = 284;
V1Speed = 11.8;

% Sun
muSun = 1.32712e+11;
%% EARTH ON LAUNCH DAY
d = 15;
m = 10;
y = 1997;
UT = 8.45;
planet = 1;

% Call Julian Day Function to obtain six orbital elements (VERIFIED)
[J0,T0,JD,h,a,e,I,Omega,omegaBar,L,omega,M] = Julian_Day_Function(d,m,y,UT,planet);
hEarthLaunch = h;
aEarthLaunch = a;
eEarthLaunch = e;
MEarthLaunch = M;
iEarthLaunch = I;
omegaEarthLaunch = omega;
OmegaEarthLaunch = Omega;

% Find Eccentric Anomaly for Earth on day of launch
[Eanomaly]=Newtons_Alg_Function(MEarthLaunch*(pi/180),eEarthLaunch*(pi/180));
Eanomaly = Eanomaly*(180/pi)

% Find the true anomaly
ThetaEarthLaunch = 2*atand(sqrt((1+eEarthLaunch)/(1-eEarthLaunch))*tand(Eanomaly/2));
if ThetaEarthLaunch < 0
    ThetaEarthLaunch = ThetaEarthLaunch + 360;
end
% Adjust to make relative to Vernal Equinox
ThetaEarthLaunch = ThetaEarthLaunch;

% Calculate the distance of earth from sun on launch day
mu = 1.32712e+11;
rEarthLaunch = (hEarthLaunch^2/mu)*(1/(1+eEarthLaunch*cosd(ThetaEarthLaunch)))

% Calculate the perifocal state vector's Rbar1 and Vbar1
PerifocalPositionVector = [cosd(ThetaEarthLaunch) sind(ThetaEarthLaunch) 0];
RbarEarth = rEarthLaunch*PerifocalPositionVector
PerifocalVelocityVector = [-sind(ThetaEarthLaunch) eEarthLaunch+cosd(ThetaEarthLaunch) 0];
VbarEarth = (muSun/hEarthLaunch)*PerifocalVelocityVector


%% VENUS ON ARRIVAL DAY
muSun = 1.32712e+11;
d = 26;
m = 4;
y = 1998;
UT = 13.7333;
planet = 3;

% Call Julian Day Function to obtain six orbital elements
[J0,T0,JD,h,a,e,I,Omega,omegaBar,L,omega,M] = Julian_Day_Function(d,m,y,UT,planet);
hVenus1 = h;
aVenus1 = a;
eVenus1 = e;
MVenus1 = M;
iVenus1 = I;
omegaVenus1 = omega;
OmegaVenus1 = Omega;


% Find Eccentric Anomaly for Earth on day of launch
Eanomaly = Newtons_Alg_Function(MVenus1*(pi/180),eVenus1*(pi/180));
EanomalyVenus = Eanomaly*(180/pi)

% Find the true anomaly
ThetaVenus1 = 2*atand(sqrt((1+eVenus1)/(1-eVenus1))*tand(EanomalyVenus/2));
if ThetaVenus1 < 0
    ThetaVenus1 = ThetaVenus1 + 360;
end
ThetaVenus1 = ThetaVenus1;

% Calculate the distance of earth from sun on launch day
rVenus1 = (hVenus1^2/muSun)*(1/(1+eVenus1*cosd(ThetaVenus1)))

% Calculate the perifocal state vector's Rbar1 and Vbar1
PerifocalPositionVector = [cosd(ThetaVenus1) sind(ThetaVenus1) 0];
RbarVenus = rVenus1*PerifocalPositionVector
PerifocalVelocityVector = [-sind(ThetaVenus1) eVenus1+cosd(ThetaVenus1) 0];
VbarVenus = (muSun/hVenus1)*PerifocalVelocityVector


%%
% Using a DCM, transform the perifocal state vectors into heliocentric
% state vecotrs using the heliocentric orbital elements found in the Julian
% Day function.
% 
% Perifocal to heliocentric ecliptic frame DCM for EarthLaunch
%incMatrix = [];
QxbarX = [-sin(OmegaEarthLaunch)*cos(iEarthLaunch)*sin(omegaEarthLaunch)+cos(OmegaEarthLaunch)*cos(omegaEarthLaunch) -sin(OmegaEarthLaunch)*cos(iEarthLaunch)*cos(omegaEarthLaunch)-cos(OmegaEarthLaunch)*sin(omegaEarthLaunch) sin(OmegaEarthLaunch)*sin(iEarthLaunch);
    cos(OmegaEarthLaunch)*cos(iEarthLaunch)*sin(omegaEarthLaunch)+sin(OmegaEarthLaunch)*cos(omegaEarthLaunch) cos(OmegaEarthLaunch)*cos(iEarthLaunch)*cos(omegaEarthLaunch)-sin(OmegaEarthLaunch)*sin(omegaEarthLaunch) -cos(OmegaEarthLaunch)*sin(iEarthLaunch);
    sin(iEarthLaunch)*sin(omegaEarthLaunch) sin(iEarthLaunch)*cos(omegaEarthLaunch) cos(iEarthLaunch)];
R1 = RbarEarth*QxbarX
VE = VbarEarth*QxbarX

% Perifocal to heliocentric ecliptic frame DCM for Venus flyby 1
QxbarX = [-sin(OmegaVenus1)*cos(iVenus1)*sin(omegaVenus1)+cos(OmegaVenus1)*cos(omegaVenus1) -sin(OmegaVenus1)*cos(iVenus1)*cos(omegaVenus1)-cos(OmegaVenus1)*sin(omegaVenus1) sin(OmegaVenus1)*sin(iVenus1);
    cos(OmegaVenus1)*cos(iVenus1)*sin(omegaVenus1)+sin(OmegaVenus1)*cos(omegaVenus1) cos(OmegaVenus1)*cos(iVenus1)*cos(omegaVenus1)-sin(OmegaVenus1)*sin(omegaVenus1) -cos(OmegaVenus1)*sin(iVenus1);
    sin(iVenus1)*sin(omegaVenus1) sin(iVenus1)*cos(omegaVenus1) cos(iVenus1)];
R2 = RbarVenus*QxbarX
V1V = VbarVenus*QxbarX


%%
% Use algorithm 5.2 to calculate spacecrafts velocity at departure of the
% planets sphere of influence VDv and the spacecrafts arrival at Venus'
% sphere of influence VAv (pg. 253)
%
%
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

% Time, in seconds, from launch dat to arrival day
deltaT = 16675200;

% 1.
% Calculate r1 and r2 using equation 5.24
rEarthLaunch = sqrt(dot(R1,R1));
rVenus1 = sqrt(dot(R2,R2));

% 2.
% Choose either prograde or retrograde trajectory
% Assuming Prograde:
Zcheck = cross(R1,R2);
Zcheck = Zcheck(3);
if Zcheck >= 0
    deltaTheta = acosd(dot(R1,R2)/(rEarthLaunch*rVenus1))
end
if Zcheck < 0
    deltaTheta = 360 - acosd(dot(R1,R2)/(rEarthLaunch*rVenus1))
end

% 3.
% Calculate A using equation 5.35
A1 = sind(deltaTheta)*sqrt((rEarthLaunch*rVenus1)/(1-cosd(deltaTheta)))

% 4.
% Use Lambert function to calculate V1 and V2 (Also, double check the A
% value calculated above
[A, V1, V2] = lambert_function(R1, R2, rEarthLaunch, rVenus1, deltaT, 'retro');
A = A
VLaunchEarth = V1
VarriveVenus = V2

% 5.
% Calculate the hyperbolic excess velocities at departure and arrival using
% equations 8.94 and 8.95
%
% Eq. 8.94
VinfDepartvec = VLaunchEarth - VE;  
VinfDepart = norm(VLaunchEarth)-norm(VE);
VinfDepart = VinfDepart/1000
fprintf('This is the Vinf relative to Earth at departure, meaning Earth \n')
fprintf('sees the spacecraft moving slower by %d km/s. \n',VinfDepart)
%
% norm(VE) = 29.879 km/s (compare to known ~30 km/s mean velocity)
% norm(VLaunchEarth) = 27.411 km/s (UNVERIFIED)

% Eq. 8.95
VinfArrivevec = VarriveVenus - V1V;
VinfArrive = norm(VarriveVenus)-norm(V1V);
VinfArrive = VinfArrive
fprintf('This is the Vinf relative to Venus on Cassini''s first flyby, \n')
fprintf('meaning Venus sees the spacecraft moving fast than it by %d km/s. \n',VinfArrive)
%
% norm(V1V) = 34.836 km/s (compare to known ~35 km/s mean velocity)
% norm(VarriveVenus) = 37.560 km/s (UNVERIFIED)



% BUILD THE TRANSFER ORBIT %
hTransfer = cross(R1,V1);
%Evec = 
%eTransfer = norm(Evec)




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%




%% VENUS FLYBY 1 ANALYSIS
% periapsis altitude of flyby 1 is 284 km

% Define known constants
muVenus = 324900;
radiusVenus = 6052;
aVenus = 108.2e+06;
eVenus = 0.0067774;

% Define approach velocity (V infinity)
VaV1 = VinfArrive;







