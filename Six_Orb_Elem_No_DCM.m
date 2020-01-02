function [Theta,r,h,Vr,e,E,I,H] = Six_Orb_Elem_No_DCM(R,V,Mu)

format long
% R = [-5000 -8000 -2100]; %Position Vector (geocentric equatorial frame)
% V = [-4 3.5 -3]; %Velocity vector (geocentric equatorial frame)
H = cross(R,V); %Angular moments vector
h = sqrt(dot(H,H)); %Specific angular momentum
I = acos(H(3)/h); %Inlcination
mu = Mu; %Earth constant
v = norm(V); %Velocity magnitude
r = norm(R); %Position magnitude
RK = [0 0 1]; %Unit vector in K direction
N = cross(RK,H); %Nodal vector
n = sqrt(dot(N,N)); %Magnitude of nodal vector
Omega = (acos(N(1)/n)); %Right ascention of the ascending node
Vr = dot(R,V)/r; %Radial velocity

%If E vector is given, then type it in here:
%E = [];
E = (1/mu)*(((v^2-(mu/r)).*R)-(r*Vr.*V));
e = norm(E);

w = (acos(dot(N,E)/(n*e))); %Argument of perigee
if E(3) < 0
    w = (2*pi) - w;
end
if E(3) > 0 
    w;
end
Theta = (acos(dot((E/e),(R/r)))); %True anomaly
if Vr < 0
    Theta = (2*pi) - Theta;
end
if Vr > 0
    Theta = Theta;
end
Theta;


end