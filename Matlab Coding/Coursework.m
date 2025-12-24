% Henry_Hodgson_20752939.m. 
% egyhh13@nottingham.ac.uk


%% Q1 - VARIABLES [9 MARKS]
x = 2.3
% a)
a = sin(3*x)
% b)
b = cos(2*x) + sin(2*x)
% c)
c = x^6 - 4*x^4 + 2*x^2 - 7
% d)
d = exp(sin(x)*cos(x))
% e)
e = (x+2)/(x-2)+(x^2+1)/(x^3)
% f)
f = sin(3*x)^2 + cos(x)^2 + sin(x)^2
% g)
e = sqrt(x^2 + 1)*x^3 + (x+2)/(x-2)
%% Q2 - VECTORS 1 [9 MARKS]
% a)  
x = linspace(-3*pi, 3*pi, 150)
% b)
a = cos(x)
% c)
b = sin(3*x)
% d) 
c = exp(-cos(x)^2)
% e)
d = b * c
% f)
m = mean(b)
s = std(b)
% g)
% small rounding errors from the approximated value of pi used.
%% Q3 - VECTORS 2 [7 MARKS]
% a)
t = 0:2e-3:40e-3
% b)
V = 5 + 3*sin(100*pi*t)
% c)
i = 0.2 + 0.1*cos(100*pi*t)
% d)
p = i.*V
% e)
max(p)

%% Q4 - EQUATIONS [11 MARKS]
% a)
a = 1
b = -1.2e7
c = -2e13
ra = (-b + sqrt(b^2-4*a*c))/2*a
rb = (-b - sqrt(b^2-4*a*c))/2*a
% b)
% selection: ra, as it is the only positive root and a length can only be
% positive.
% c)
F = (6.674e-11*5.97e24*1000)/ra^2
% d)
v = sqrt((6.674e-11*5.97e24)/ra))
% e)
T = (2*pi*ra)/v


%% Q5 - TEXT FORMAT, PRINT AND LOOPS [12 MARKS]
% a) How many locations?

clear
clc

numLocations = input("Enter number of locations:")
% b) Info for each location:
locations = struct('name', cell(1,numLocations), ...
                   'temperature', cell(1,numLocations), ...
                   'humidity', cell(1,numLocations), ...
                   'windspeed', cell(1,numLocations));

for k = 1:numLocations
    prompt = sprintf("Enter name for location %d: ", k);
    name = input(prompt, 's');             % text

    temp = input(sprintf("Enter temperature for %s: ", name));   % numeric
    hum  = input(sprintf("Enter humidity for %s: ", name));      % numeric
    ws   = input(sprintf("Enter wind speed for %s: ", name));    % numeric

    locations(k).name = name;
    locations(k).temperature = temp;
    locations(k).humidity = hum;
    locations(k).windspeed = ws;
end


%% Q6 - PROGRAM FLOW [14 MARKS]
clear



%% Q7 - FOR LOOPS, SWITCH STATEMENTS AND DISPLAYING DATA [16 MARKS]
clear



%% PRESENTATION AND FORMATTING [15 MARKS]
% You do not need to put anything in this section
% Marks will be given for neat, tidy presentation, sensible and informative
% variable names, commenting and remembering to put your name and email address at the top of the
% script
