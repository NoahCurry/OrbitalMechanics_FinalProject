clear all; close all; clc

%% Compare to algorithm 8.2 - Example 8.8 .m file
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example_8_08
% ~~~~~~~~~~~~
%{
  This program uses Algorithm 8.2 to solve Example 8.8.
 
  mu           - gravitational parameter of the sun (km^3/s^2)
  deg          - conversion factor between degrees and radians
  pi           - 3.1415926...

  planet_id    - planet identifier:
                  1 = Mercury
                  2 = Venus
                  3 = Earth
                  4 = Mars
                  5 = Jupiter
                  6 = Saturn
                  7 = Uranus
                  8 = Neptune
                  9 = Pluto

  year         - range: 1901 - 2099
  month        - range: 1 - 12
  day          - range: 1 - 31
  hour         - range: 0 - 23
  minute       - range: 0 - 60
  second       - range: 0 - 60 

  depart       - [planet_id, year, month, day, hour, minute, second]
                 at departure
  arrive       - [planet_id, year, month, day, hour, minute, second]
                 at arrival

  planet1      - [Rp1, Vp1, jd1]
  planet2      - [Rp2, Vp2, jd2]
  trajectory   - [V1, V2]

  coe          - orbital elements [h e RA incl w TA]
                 where
                   h    = angular momentum (km^2/s)
                   e    = eccentricity
                   RA   = right ascension of the ascending
                          node (rad)
                   incl = inclination of the orbit (rad)
                   w    = argument of perigee (rad)
                   TA   = true anomaly (rad)
                   a    = semimajor axis (km)

  jd1, jd2     - Julian day numbers at departure and arrival
  tof          - time of flight from planet 1 to planet 2 (days)

  Rp1, Vp1     - state vector of planet 1 at departure (km, km/s)
  Rp2, Vp2     - state vector of planet 2 at arrival (km, km/s)
  R1, V1       - heliocentric state vector of spacecraft at
                 departure (km, km/s)
  R2, V2       - heliocentric state vector of spacecraft at
                 arrival (km, km/s)

  vinf1, vinf2 - hyperbolic excess velocities at departure
                 and arrival (km/s)

  User M-functions required: interplanetary, coe_from_sv,
                             month_planet_names
%}
% ---------------------------------------------

%clear all;
global mu
mu  = 1.327124e11;
deg = pi/180;

%...Data declaration for Example 8.8:

%...Departure
planet_id = 3;
year      = 1997;
month     = 10;
day       = 15;
hour      = 8;
minute    = 27;
second    = 0;
depart = [planet_id  year  month  day  hour  minute  second];

%...Arrival
planet_id = 2;
year      = 1998;
month     = 4;
day       = 26;
hour      = 13;
minute    = 45;
second    = 0;
arrive = [planet_id  year  month  day  hour  minute  second];

%...

%...Algorithm 8.2:
[planet1, planet2, trajectory] = interplanetary(depart, arrive);

R1  = planet1(1,1:3);
Vp1 = planet1(1,4:6);
jd1 = planet1(1,7);

R2  = planet2(1,1:3);
Vp2 = planet2(1,4:6);
jd2 = planet2(1,7);

V1  = trajectory(1,1:3);
V2  = trajectory(1,4:6);

tof = jd2 - jd1;

%...Use Algorithm 4.2 to find the orbital elements of the
%   spacecraft trajectory based on [Rp1, V1]...
coe  = coe_from_sv(R1, V1, mu);
%   ... and [R2, V2]
coe2 = coe_from_sv(R2, V2, mu);

%...Equations 8.94 and 8.95:
vinf1 = V1 - Vp1;
vinf2 = V2 - Vp2;

%...Echo the input data and output the solution to
%   the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Example 8.8')
fprintf('\n\n Departure:\n');
%fprintf('\n   Planet: %s', planet_name(depart(1)))
fprintf('\n   Year  : %g', depart(2))
%fprintf('\n   Month : %s', month_name(depart(3)))
fprintf('\n   Day   : %g', depart(4))
fprintf('\n   Hour  : %g', depart(5))
fprintf('\n   Minute: %g', depart(6))
fprintf('\n   Second: %g', depart(7))
fprintf('\n\n   Julian day: %11.3f\n', jd1)
fprintf('\n   Planet position vector (km)    = [%g  %g  %g]', ...
                                               R1(1),R1(2), R1(3))

fprintf('\n   Magnitude                      = %g\n', norm(R1))

fprintf('\n   Planet velocity (km/s)         = [%g  %g  %g]', ...
                                 Vp1(1), Vp1(2), Vp1(3))

fprintf('\n   Magnitude                      = %g\n', norm(Vp1))

fprintf('\n   Spacecraft velocity (km/s)     = [%g  %g  %g]', ...
                                               V1(1), V1(2), V1(3))

fprintf('\n   Magnitude                      = %g\n', norm(V1))

fprintf('\n   v-infinity at departure (km/s) = [%g  %g  %g]', ...
                                       vinf1(1), vinf1(2), vinf1(3))

fprintf('\n   Magnitude                      = %g\n', norm(vinf1))

fprintf('\n\n Time of flight = %g days\n', tof)

fprintf('\n\n Arrival:\n');
%fprintf('\n   Planet: %s', planet_name(arrive(1)))
fprintf('\n   Year  : %g', arrive(2))
%fprintf('\n   Month : %s', month_name(arrive(3)))
fprintf('\n   Day   : %g', arrive(4))
fprintf('\n   Hour  : %g', arrive(5))
fprintf('\n   Minute: %g', arrive(6))
fprintf('\n   Second: %g', arrive(7))
fprintf('\n\n   Julian day: %11.3f\n', jd2)
fprintf('\n   Planet position vector (km)   = [%g  %g  %g]', ...
                                              R2(1), R2(2), R2(3))

fprintf('\n   Magnitude                     = %g\n', norm(R1))

fprintf('\n   Planet velocity (km/s)        = [%g  %g  %g]', ...
                                  Vp2(1), Vp2(2), Vp2(3))

fprintf('\n   Magnitude                     = %g\n', norm(Vp2))

fprintf('\n   Spacecraft Velocity (km/s)    = [%g  %g  %g]', ...
                                              V2(1), V2(2), V2(3))

fprintf('\n   Magnitude                     = %g\n', norm(V2))

fprintf('\n   v-infinity at arrival (km/s)  = [%g  %g  %g]', ...
                                     vinf2(1), vinf2(2), vinf2(3))

fprintf('\n   Magnitude                     = %g', norm(vinf2))

fprintf('\n\n\n Orbital elements of flight trajectory:\n')

fprintf('\n  Angular momentum (km^2/s)                   = %g',...
                                                           coe(1))
fprintf('\n  Eccentricity                                = %g',...
                                                           coe(2))
fprintf('\n  Right ascension of the ascending node (deg) = %g',...
                                                       coe(3)/deg)
fprintf('\n  Inclination to the ecliptic (deg)           = %g',...
                                                       coe(4)/deg)
fprintf('\n  Argument of perihelion (deg)                = %g',...
                                                       coe(5)/deg)
fprintf('\n  True anomaly at departure (deg)             = %g',...
                                                       coe(6)/deg)
fprintf('\n  True anomaly at arrival (deg)               = %g\n', ...
                                                      coe2(6)/deg)
fprintf('\n  Semimajor axis (km)                         = %g',...
                                                           coe(7))
% If the orbit is an ellipse, output the period:
if coe(2) < 1
    fprintf('\n  Period (days)                               = %g', ...
                                      2*pi/sqrt(mu)*coe(7)^1.5/24/3600)
end
fprintf('\n-----------------------------------------------------\n')
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
fprintf('\n Earth relative excess velocity (Vinf) at departure: %g (km/s) \n\n',norm(V1)-norm(Vp1))
%
%
% Generate positions throughout the traveled transfer orbit
thetatheta = linspace(coe(6),coe2(6),100);
rr = (3832132648.38^2/mu)*(1./(1+0.210328.*cos(thetatheta)));
[thetatheta,rr] = pol2cart(thetatheta,rr);

% Plot the transfer
figure(1)
plot(thetatheta,rr,'k'),grid on; hold on;
title('Cassini Transfer Orbit 1: Non-Hohmann Earth to Venus Trajectory')

%% Plot Earth and Venus orbits
%
%
% Here, I use interplanetary.m function to calculate the departure Velocity of 
% Cassini. I then build the transfer orbit and attempt to plot it...
%
%
% Define departure and arrival components to pass to the .m function
depart = [3,1997,10,15,8,27,0];
arrive = [2,1998,4,26,13,45,0];

% Calculate the position and velocity vectors and magnitudes for Cassini on
% departure day and first arrival at Venus. 
[planet1, planet2, trajectory] = interplanetary(depart, arrive);
R1 = [planet1(1) planet1(2) planet1(3)];
r1 = norm(R1);
R2 = [planet2(1) planet2(2) planet2(3)];
r2 = norm(R2);
VdVvec = [trajectory(1) trajectory(2) trajectory(3)];
VdV = norm(VdVvec);
VaVvec = [trajectory(4) trajectory(5) trajectory(6)];
VaV = norm(VaVvec);

% Define the transfer orbit
[Theta,r,h,Vr,e] = Six_Orb_Elem_No_DCM(R1,VdVvec,mu);
VrE = Vr;
hTran1 = h;
eTran1 = e;
thetaEarth = Theta;
rrEarth = r;
[thetaEarth,rrEarth] = pol2cart(thetaEarth,rrEarth);
    % Calculate semi-major axis
    rp1 = (hTran1^2/mu)*(1/(1+eTran1));
    ra1 = (hTran1^2/mu)*(1/(1-eTran1));
    aTran1 = (rp1+ra1)/2;
    % Calculate period
    period1 = ((2*pi)/sqrt(mu))*aTran1^(3/2);
        % Find true anomaly for Venus at arrival
        [Theta,r,h,Vr,e] = Six_Orb_Elem_No_DCM(R2,VaVvec,mu);
        VrV = Vr;
        thetaVenus = Theta;
        rrVenus = r;
       [thetaVenus,rrVenus] = pol2cart(thetaVenus,rrVenus);
       
plot(thetaEarth,rrEarth,'bo','MarkerSize',8,'MarkerFaceColor','b')
plot(thetaVenus,rrVenus,'ro','MarkerSize',8,'MarkerFaceColor','r')

% Earth
hEarth = (30.4)*(1.46e+08);
eEarth = 0.0167;
ThetaEarth = linspace(0,2*pi,100);
rrEarth = (hEarth^2/mu)*(1./(1+eEarth*cos(ThetaEarth)));
[ThetaEarth,rrEarth] = pol2cart(ThetaEarth,rrEarth);
plot(ThetaEarth,rrEarth,'b--')


% Venus
hVenus = norm(cross([1.36218e+07  -1.07934e+08  -2.26085e+06],[34.51  4.26073  -1.93382]));
eVenus = 0.0067;
ThetaVenus = linspace(0,2*pi,100);
rrVenus = (hVenus^2/mu)*(1./(1+eVenus*cos(ThetaVenus)));
[ThetaVenus,rrVenus] = pol2cart(ThetaVenus,rrVenus);
plot(ThetaVenus,rrVenus,'r--')

% Sol
plot(0,0,'ro','MarkerSize',12,'MarkerFaceColor',[0.9100    0.4100    0.1700])
legend('Flown Transfer Orbit','Earth at Launch','Venus at Arrival','Earth Orbit','Venus Orbit','Sol','Location','SouthEast')
%% Venus Flyby 1
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% V infinity is found by norm(V1)-norm(Vp1). This is for leaving Earth and
% and going to Venus. I need to take this value and use it as my approach
% velocity for the flyby.
%
% The values that I need to find are listed in Rori's section of the
% report. Calculate those same ones and put the in the report. Also,
% compare those values with data that you can find in Venus trajectory
% papers.
%
% ORBITAL ELEMENTS EQUATIONS:
%           efly1 = 1 + ((rperiapsis*Vinf^2)/(muVenus))
%           delta (turn angle) = 2*asin(1/efly1)
%           h = rperiapsis*sqrt(vinf^2 + ((2*muVenus)/rperiapsis)
%           Vperiapsis = sqrt(norm(Vinf)^2+(2*muVenus/rperiapsis))
%
%
muVenus = 324900;
Vinf1 = 6.03; % Taken from Cassini Manuever Expreience: Launch and Early C.
rVenus = 6052;
Alt1 = 284;
rperiapsis = rVenus+Alt1;
efly1 = 1 + ((rperiapsis*Vinf1^2)/(muVenus));
delta1 = 2*asin(1/efly1)*(180/pi);
h1 = rperiapsis*sqrt(Vinf1^2 + ((2*muVenus)/rperiapsis));
Vperiapsis = sqrt(norm(Vinf1)^2+(2*muVenus/rperiapsis));

fprintf('Venus Flyby 1 Orbital Parameters: \n')
fprintf('Eccentricity =         %g \n',efly1)
fprintf('Turn Angle =           %g degrees \n',delta1)
fprintf('Angular Momentum =     %g (km^2/s^2) \n',h1)
fprintf('Periapsis Velocity =   %g (km/s) \n',Vperiapsis)

%
%
% NOTES FROM MEHRAN:
% Calculate the hyperbolic exit energy using energy eq with r going to
% zero and then from that you can get a (semi major). Use mu earth. 
%
%
%% VENUS-VENUS ORBIT
%
% The known period is 425 days, which = 36720000 seconds.
TVV = 36720000;

% Calculate semi-major axis
aVV = ((sqrt(mu)*TVV)/(2*pi))^(2/3);

% Use julian day function to get h and e and then calculate the distance of
% Venus from the sun on the day of flyby1
d = 26;
m = 4;
y = 1998;
UT = 12;
planet = 3;

% Call Julian Day Function to obtain six orbital elements (VERIFIED)
[J0,T0,JD,h,a,e,I,Omega,omegaBar,L,omega,M] = Julian_Day_Function(d,m,y,UT,planet);
hVVcheck = h;
aVVcheck = a;
eVVcheck = e;
MVV = M;

% Calcualte distance from sun
% Find Eccentric Anomaly for Earth on day of launch
[Eanomaly]=Newtons_Alg_Function(MVV*(pi/180),eVVcheck*(pi/180));
Eanomaly = Eanomaly*(180/pi)

% Find the true anomaly
ThetaVV = 2*atand(sqrt((1+eVVcheck)/(1-eVVcheck))*tand(Eanomaly/2));
if ThetaVV < 0
    ThetaVV = ThetaVV + 360;
end

% Calculate the distance of earth from sun on launch day
rVV = (hVVcheck^2/mu)*(1/(1+eVVcheck*cosd(ThetaVV)))
    
% Calculate radius of perigee and apogee
rpVV = (rVenus + 603) + rVV;
raVV = (2*aVV) - rpVV;

% Calculate eccentricity and angular momentum
eVV = (raVV - rpVV)/(raVV + rpVV);
hVV = sqrt(mu*r*(1+eVV));

% Generate postion values for both calculated values (check and regular)
 ThetaVV = linspace(0,2*pi,100);
% rVVcheck = (hVVcheck^2/mu)*(1./(1+eVVcheck.*cos(ThetaVV)));
% polarplot(ThetaVenus,rVVcheck,'m')
% ThetaVV = linspace(0,2*pi,100);
 rVV = (hVV^2/mu)*(1./(1+eVV.*cos(ThetaVV)));
[ThetaVV,rVV] = pol2cart(ThetaVV,rVV);
plot(ThetaVenus,rVV,'k')

