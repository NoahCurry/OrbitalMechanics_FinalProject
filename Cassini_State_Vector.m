close all; clear all; clc    
% Define planet
%   Earth = 1
%   Mars = 2
%   Venus = 3
    planet = 1;
    
% Define initial day/time
    d = 1;
    m = 1;
    y = 2005;
    UT = 12;

% Create mean anomaly and eccentricity vectors
    listM = [];
    liste = [];
for month = 1:12
    for day = 1:30
        [J0,T0,JD,h,a,e,I,Omega,omegaBar,L,omega,M] = Julian_Day_Function(day,month,y,UT,planet);
        listM(day,month) = M;
        liste(day,month) = e;
    end
end


%% Find the eccentric anomaly (E)
Eanomaly = [];
for monthCount = 1:12
    for dayCount = 1:30
        Eanomaly(dayCount,monthCount)= Newtons_Alg_Function(listM(dayCount,monthCount)*(pi/180),liste(dayCount,monthCount));
    end
end
Eanomaly = Eanomaly.*(180/pi);


%% Find the true anomaly (Theta)
Theta = [];
for months = 1:12
    for days = 1:30
        Theta(days,months) = 2*atand(sqrt((1+liste(days,months))/(1-liste(days,months)))*tand(Eanomaly(days,months)/2));
        if Theta(days,months) > 0
            while Theta(days,months) > 0
                Theta(days,months) = Theta(days,months) - 360;
            end
        end
        if Theta(days,months) < 0
            while Theta(days,months) < 0
                Theta(days,months) = Theta(days,months) + 360;
            end
        end
    end
end


%% Calculate the position (r) and itemize the r matrix into months
mu = 1.32712e+11;
r = [];
for monthz = 1:12
    for dayz = 1:30
        r(dayz,monthz) = (h^2/mu)*(1/(1+(e*cosd(Theta(dayz,monthz)))));
    end
end
Jan = r(1:end,1);
Feb = r(1:end,2);
March = r(1:end,3);
April = r(1:end,4);
May = r(1:end,5);
June = r(1:end,6);
July = r(1:end,7);
Aug = r(1:end,8);
Sep = r(1:end,9);
Oct = r(1:end,10);
Nov = r(1:end,11);
Dec = r(1:end,12);


%% Sort through months to find maximum values and create maxMonths vector
maxJan = max(Jan);
maxFeb = max(Feb);
maxMarch = max(March);
maxApril = max(April);
maxMay = max(May);
maxJune = max(June);
maxJuly = max(July);
maxAug = max(Aug);
maxSep = max(Sep);
maxOct = max(Oct);
maxNov = max(Nov);
maxDec = max(Dec);

maxMonths = [maxJan maxFeb maxMarch maxApril maxMay maxJune maxJuly maxAug ...
    maxSep maxOct maxNov maxDec];


%% Determine which month has max r value
for maxCount = 1:length(maxMonths)
    Month = maxMonths(maxCount);
    if Month == max(maxMonths)
        maxMonth = maxCount;
    end
end


%% Determine which day in maxMonth has the max value
for maxDayCount = 1:30
    Day = July(maxDayCount);
    if Day == max(July)
        maxDay = maxDayCount;
    end
end