function [Jan,Feb,March,April,May,June,July,Aug,Sep,Oct,Nov,Dec,Theta,listh,a,liste,omega,Omega,I] = Cassini_State_Vector_Function(y,UT,planet)
%close all; clear all; clc    
% Define planet
%   Earth = 1
%   Mars = 2
%   Venus = 3
%   planet = 1;
    
% Define initial day/time
%   y = 2005;
%   UT = 12;

% Create mean anomaly and eccentricity vectors
    listM = [];
    liste = [];
    listh = [];
    lista = [];
for month = 1:12
    for day = 1:30
        [J0,T0,JD,h,a,e,I,Omega,omegaBar,L,omega,M] = Julian_Day_Function(day,month,y,UT,planet);
        listM(day,month) = M;
        liste(day,month) = e;
        lista(day,month) = a;
        listh(day,month) = h;
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
end