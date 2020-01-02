% THE SIX ORBITAL ELEMENTS
clear all; close all; clc
format long
R = [-5000 -8000 -2100]; %Position Vector (geocentric equatorial frame)
V = [-4 3.5 -3]; %Velocity vector (geocentric equatorial frame)
H = cross(R,V); %Angular moments vector
h = sqrt(dot(H,H)) %Specific angular momentum
I = acos(H(3)/h) %Inlcination
mu = 398600; %Earth constant
v = norm(V); %Velocity magnitude
r = norm(R); %Position magnitude
RK = [0 0 1]; %Unit vector in K direction
N = cross(RK,H); %Nodal vector
n = sqrt(dot(N,N)) %Magnitude of nodal vector
Omega = (acos(N(1)/n)) %Right ascention of the ascending node
Vr = dot(R,V)/r; %Radial velocity

%If E vector is given, then type it in here:
%E = [];
E = (1/mu)*(((v^2-(mu/r)).*R)-(r*Vr.*V));
e = norm(E)

w = (acos(dot(N,E)/(n*e))); %Argument of perigee
if E(3) < 0
    w = (2*pi) - w
end
if E(3) > 0 
    w
end
Theta = (acos(dot((E/e),(R/r)))); %True anomaly
if Vr < 0
    Theta = (2*pi) - Theta;
end
if Vr > 0
    Theta = Theta;
end
Theta


%%
%
%
%
% KEPLERS TIME PROBLEM (components)
%
% MUST RUN "SIX ORBITAL ELEMENTS" FIRST
T = ((2*pi)/mu^2)*(h/sqrt(1-e^2))^3 %Period
E0 = 2*atan(sqrt((1-e)/(1+e))*tan((Theta)/2)) %Initial eccentric anomaly
M0 = E0-(e*sin(E0)) %Initial mean anomaly
t0 = (M0/(2*pi))*T %Initial time
tstep = 50*60; %Change in time
t1 = t0 + tstep; %Final time
M = (2*pi)*(t1/T) %Final mean anomaly
%
%
% Using M and e, run Newton's Algorithm to obtain Ee
%
%
Ee = -0.5914;
ThetaFin = 2*(atan(tan(Ee/2)/(sqrt((1-e)/(1+e)))))
%
%
%
% MUST RUN KEPLERS TIME PROBLEM FIRST
%
% Q TRANSFORMATION MATRICES
PeriPosMatrix = [cos(ThetaFin); sin(ThetaFin); 0];
PeriVelMatrix = [-sin(ThetaFin); e+cos(ThetaFin); 0];
Rxbar = ((h^2/mu)*(1/(1+e*cos(ThetaFin))))*PeriPosMatrix;
Vxbar = (mu/h)*PeriVelMatrix;
QxbarX = [-sin(Omega)*cos(I)*sin(w)+cos(Omega)*cos(w) -sin(Omega)*cos(I)*cos(w)-cos(Omega)*sin(w) sin(Omega)*sin(I);
    cos(Omega)*cos(I)*sin(w)+sin(Omega)*cos(w) cos(Omega)*cos(I)*cos(w)-sin(Omega)*sin(w) -cos(Omega)*sin(I);
    sin(I)*sin(w) sin(I)*cos(w) cos(I)];
QXxbar = [-sin(Omega)*cos(I)*sin(w)+cos(Omega)*cos(w) cos(Omega)*cos(I)*sin(w)+sin(Omega)*cos(w) sin(I)*sin(w);
    -sin(Omega)*cos(I)*cos(w)-cos(Omega)*sin(w) cos(Omega)*cos(I)*cos(w)-sin(Omega)*sin(w) sin(I)*cos(w);
    sin(Omega)*sin(I) -cos(Omega)*sin(I) cos(I)];
Rx = QxbarX*Rxbar
Vx = QxbarX*Vxbar
