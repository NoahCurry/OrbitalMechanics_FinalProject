function [J0,T0,JD,h,a,e,I,Omega,omegaBar,L,omega,M] = Julian_Day_Function(d,m,y,UT,planet)
% Calculate the Julian information and six orbital elements for the orbit of a
% planet if given the specific date of interest. 
% Input arguments: d (day), m (month), y (year), UT (universal time hour),
% planet: Earth = 1
%         Mars = 2
%         Venus = 3
% Output: J0, T0, JD, h, a, e, i, Omega, and omega.

format long
%% Calculate J0 (Julian Day Number for Specific Date)
A = (m+9)/12;
INT2 = floor(A);
B = (7*(y+INT2))/4;
INT1 = floor(B);
C = (275*m)/9;
INT3 = floor(C);
J0 = (367*y)-INT1+INT3+d+1721013.5;

%J0 = (367*y) - ((7*(y+((m+9)/12)))/4) + ((275*m)/9) + d + 1721013.5;

%% Calculate JD (Julian Day Number at x hr UT)
JD = J0 + (UT/24);

%% Calculate T0 (Time Elapsed Between J0 and J2000)
T0 = (JD - 2451545)/36525;


%% Calculate ThetaG (Greenwhich Sidereal Time at any UT)
ThetaG0 = 100.4606184 + (3600.77004*T0) + (0.000387933*T0^2) - ((2.583e-08)*T0^3);
while ThetaG0 > 360
    ThetaG0 = ThetaG0 - 360;
end
ThetaG = ThetaG0 + (360.98564724*(UT/24));

%% Calculate New Elements
mu = 1.32712e+11;
if planet == 1 % Earth
    aa = 1.00000261;
    ee = 0.01671123;
    ii = -0.00001531;
    OMEGA = 0;
    omegaBar1 = 102.93768193;
    adot = 0.00000562;
    edot = -0.00004392;
    idot = -0.01294668;
    OmegaDot = 0;
    OmegaBarDot = 0.32327364;
    LL = 100.46457166;
    Ldot = 35999.37244981;
end
if planet == 2 % Mars
    aa = 1.52371034;
    ee = 0.09339410;
    ii = 1.84969142;
    OMEGA = 49.55953891;
    adot = 0.0001847;
    edot = 0.00007882;
    idot = -0.008131313;
    omegaBar1 = -23.94362959;
    OmegaDot = -0.29257343;
    OmegaBarDot = 0.44441088;
    LL = -4.55343205;
    Ldot = 19140.30268499;
end
if planet == 3 % Venus
    aa = 0.72333566;
    ee = 0.00677672;
    ii = 3.39467605;
    OMEGA = 76.67984255;
    adot = 0.00000390;
    edot = -0.00004107;
    idot = -0.00078890;
    omegaBar1 = 131.60246718;
    OmegaDot = -0.27769418;
    OmegaBarDot = 0.00268329;
    LL = 181.97909950;
    Ldot = 58517.81538729;
end

a = aa + adot*T0;
a = a*(1.49597871e+08);
e = ee + edot*T0;
I = ii + idot*T0;
L = LL + Ldot*T0;
Omega = OMEGA + OmegaDot*T0;
omegaBar = omegaBar1 + OmegaBarDot*T0;
omega = omegaBar - Omega;
M = L - omegaBar;
h = sqrt(mu*a*(1-e^2));

%% Adjust new elements to fall within 0 < x < 360
if omegaBar > 360
    while omegaBar > 360
        omegaBar = omegaBar -360;
    end
end
if omegaBar < 0
    while omegaBar < 0
        omegaBar = omegaBar + 360;
    end
end
if I > 360
    while I > 0
        I = I -360;
    end
end
if I < 0
    while I < 0
        I = I + 360;
    end
end
if Omega > 360
    while Omega > 0
        Omega = Omega -360;
    end
end
if Omega < 0
    while Omega < 0
        Omega = Omega + 360;
    end
end
if L > 360
    while L > 0
        L = L -360;
    end
end
if L < 0
    while L < 0
        L = L + 360;
    end
end
if M > 360
    while M > 0
        M = M -360;
    end
end
if M < 0
    while M < 0
        M = M + 360;
    end
end




end
