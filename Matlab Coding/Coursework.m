% Henry_Hodgson_20752939.m. 
% egyhh13@nottingham.ac.uk


%% Q1 - VARIABLES [9 MARKS]
clear
clc

x = 2.3;
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
clear
clc

% a)  
x = linspace(-3*pi, 3*pi, 150)
% b)
a = cos(x)
% c)
b = sin(3*x)
% d) 
c = exp(-cos(x).^2)
% e)
d = b .* c
% f)
m = mean(b)
s = std(b)
% g)
% small rounding errors from the approximated value of pi used.
%% Q3 - VECTORS 2 [7 MARKS]
clear
clc

% a)
t = 0:2e-3:40e-3
% b)
V = 5 + 3*sin(100*pi*t)
% c)
i = 0.2 + 0.1*cos(100*pi*t)
% d)
p = i.*V
% e)
MaxP = max(p);
fprintf('Max instantanious power: %d \n', MaxP)

%% Q4 - EQUATIONS [11 MARKS]
clear
clc

% a)
a = 1
b = -1.2e7
c = -2e13
ra = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
rb = (-b - sqrt(b^2 - 4*a*c)) / (2*a)
% b)
% selection: ra, as it is the only positive root and a length can only be
% positive.
% c)
F = (6.674e-11*5.97e24*1000)/ra^2
% d)
v = sqrt((6.674e-11*5.97e24)/ra)
% e)
T = (2*pi*ra)/v


%% Q5 - TEXT FORMAT, PRINT AND LOOPS [12 MARKS]
clear
clc

% a) How many locations?
numLocations = input("Enter number of locations: ");

% Preallocate variables
Location = strings(numLocations,1);
Temperature = zeros(numLocations,1);
Humidity = zeros(numLocations,1);
WindSpeed = zeros(numLocations,1);

% b) Info for each location:
for k = 1:numLocations
    Location(k) = input("Enter name for location: ", 's');
    
    Temperature(k) = input(sprintf("Enter temperature for %s: ", Location(k)));
    Humidity(k)    = input(sprintf("Enter humidity for %s: ", Location(k)));
    WindSpeed(k)   = input(sprintf("Enter wind speed for %s: ", Location(k)));
end

% Create table
WeatherTable = table(Location, Temperature, Humidity, WindSpeed);

% e) Calculate averages
avgTemp     = mean(Temperature);
avgHumidity = mean(Humidity);
avgWind     = mean(WindSpeed);

% Create a row for averages
AverageRow = table("Average", avgTemp, avgHumidity, avgWind, ...
    'VariableNames', WeatherTable.Properties.VariableNames);

% Add the average row to the table
WeatherTable = [WeatherTable; AverageRow];

% Display the full table
disp("Overall Weather Data:")
disp(WeatherTable)

% f) The script relies on being manually updated and then the data is only
% relevant for a certain amount of time before conditions change. It may be
% more useful to have the script automatically take current data from a source and
% create a table every hour, this would save having to key it in every
% time.
%% Q6 - PROGRAM FLOW [14 MARKS]
% b/c) Impliment code + comment

% NOTE: getAltitude, lowerLandingGear, activatebeacon are assumed system functions


clear;
clc;

TimeDelay = 3;
PitchAngle = 0;   % represents thrust
gearLowered = false;

while true
    Alt1 = getAltitude();   % altitude reading from computer

    % End loop when drone is below 50 m
    if Alt1 < 50
        break;
    end

    % Lower landing gear below 500 m (once)
    if Alt1 < 500 && ~gearLowered
        lowerLandingGear();
        gearLowered = true;
    end

    pause(TimeDelay);

    Alt2 = getAltitude();   % second altitude reading

    % Descent rate (positive = descending)
    Vd = -(Alt2 - Alt1) / TimeDelay;

    % Control logic
    if Vd > 7 && Vd < 8
        % correct descent rate → do nothing
    elseif Vd < 7
        PitchAngle = PitchAngle - 1;   % descending too slowly
    elseif Vd > 8
        PitchAngle = PitchAngle + 1;   % descending too quickly
    end
end

PitchAngle = PitchAngle + 4;

while true % new descent rate at 1m/s
    Alt1 = getAltitude();   % current altitude

    % Stop loop at 5 m (near touchdown)
    if Alt1 <= 5
        break;
    end

    pause(TimeDelay);

    Alt2 = getAltitude();   % altitude after delay

    % Compute descent rate (positive = descending)
    Vd = -(Alt2 - Alt1) / TimeDelay;

    % Control logic for slow descent
    if Vd > 0.5 && Vd < 1.5
        % Correct descent rate → do nothing
    elseif Vd < 1
        PitchAngle = PitchAngle - 1;   % descending too slowly
    elseif Vd > 1
        PitchAngle = PitchAngle + 1;   % descending too quickly
    end
end

PitchAngle = PitchAngle + 2;

while true % new decent rate at 0.5m/s
    Alt1 = getAltitude();   % current altitude

    % Stop loop at 5 m (near touchdown)
    if Alt1 <= 0.1
        break;
    end

    pause(TimeDelay);

    Alt2 = getAltitude();   % altitude after delay

    % Compute descent rate (positive = descending)
    Vd = -(Alt2 - Alt1) / TimeDelay;

    % Control logic for slow descent
    if Vd > 0.3 && Vd < 0.7
        % Correct descent rate → do nothing
    elseif Vd < 0.5
        PitchAngle = PitchAngle - 1;   % descending too slowly
    elseif Vd > 0.5
        PitchAngle = PitchAngle + 1;   % descending too quickly
    end
end

motor1 = 0;
motor2 = 0;
motor3 = 0;
motor4 = 0;

activatebeacon();

%% d) simulated estimate for landing gear deployment:
clear
clc
StartAlt = 3000;
LandingGearAlt = 500;
MaxDecentSpeed = 8;
MinDecentSpeed = 7;
MinTime = (StartAlt - LandingGearAlt)/MaxDecentSpeed;
MaxTime = (StartAlt - LandingGearAlt)/MinDecentSpeed;
MinMax = [MinTime, MaxTime];
AvgTime = mean(MinMax);

fprintf('Minimum time to deploy: %.2f seconds\n', MinTime);
fprintf('Maximum time to deploy: %.2f seconds\n', MaxTime);
fprintf('Average time to deploy: %.2f seconds\n', AvgTime);

%% e) simulated estimate for landing:
StartAlt = 3000;
Stage1Alt = 50;
MaxDecentSpeed = 8;
MinDecentSpeed = 7;
MinS1Time = (StartAlt - Stage1Alt)/MaxDecentSpeed;
MaxS1Time = (StartAlt - Stage1Alt)/MinDecentSpeed;

StartAlt = 50;
Stage2Alt = 5;
MaxDecentSpeed = 1.5;
MinDecentSpeed = 0.5;
MinS2Time = (StartAlt - Stage2Alt)/MaxDecentSpeed;
MaxS2Time = (StartAlt - Stage2Alt)/MinDecentSpeed;

StartAlt = 5;
Stage3Alt = 0;
MaxDecentSpeed = 0.7;
MinDecentSpeed = 0.3;
MinS3Time = (StartAlt - Stage3Alt)/MaxDecentSpeed;
MaxS3Time = (StartAlt - Stage3Alt)/MinDecentSpeed;

MaxTotalTime = MaxS1Time + MaxS2Time + MaxS3Time;
MinTotalTime = MinS1Time + MinS2Time + MinS3Time;
MaxMin = [MaxTotalTime, MinTotalTime];
AvgTotalTime = mean(MaxMin);

fprintf('Minimum time to landing: %.2f seconds\n', MinTotalTime);
fprintf('Maximum time to landing: %.2f seconds\n', MaxTotalTime);
fprintf('Average time to landing: %.2f seconds\n', AvgTotalTime);
%% Q7 - FOR LOOPS, SWITCH STATEMENTS AND DISPLAYING DATA [16 MARKS]
clear
clc

ExhaustVelocityChoice = menu('Choose an exhaust velocity (m/s)', '2000', '2500', '3000');

switch ExhaustVelocityChoice
    case 1
        ExhaustVelocity = 2000; % creates the largest possible array (10)
    case 2
        ExhaustVelocity = 2500;
    case 3
        ExhaustVelocity = 3000;
end

fprintf('Exhaust Velocity: %d m/s \n', ExhaustVelocity);

% b)
MassFR = 1;
i = 1;

% preallocate size of array for efficiency

ArraySize = ceil(19000 / ExhaustVelocity);

MassFRvals = zeros(1, ArraySize);
Thrustvals = zeros(1, ArraySize); 
Ratingvals = strings(1, ArraySize); 

while true
    Thrust = MassFR * ExhaustVelocity;

    % Classify thrust rating
    if Thrust < 10000
        Rating = "Low";
    elseif Thrust < 15000
        Rating = "Medium";
    else
        Rating = "High";
    end

    % Store results
    MassFRvals(1:i) = MassFR;
    Thrustvals(1:i) = Thrust;
    Ratingvals(1:i) = Rating;

    if Thrust > 19000
        break
    end
    fprintf('Mass Flow Rate: %d kg/s | Thrust: %d N | Level: %s \n', MassFR, Thrust, Rating);

    MassFR = MassFR + 1;
    i = i + 1;
end

PrintChoice = menu('save results as txt file?', 'Yes', 'No');

if PrintChoice == 1
     fileID = fopen('thrust_results.txt', 'w');

    % Write header
    fprintf(fileID, 'Exhaust Velocity: %d m/s\n\n', ExhaustVelocity);
    fprintf(fileID, 'Mass Flow Rate | Thrust (N) | Rating\n');

    % Write each row
    for k = 1:length(MassFRvals)
        fprintf(fileID, '%14d | %10d | %s\n', MassFRvals(k), Thrustvals(k), Ratingvals(k));
    end

    fclose(fileID);

    open('thrust_results.txt');   % opens the file

else
    fprintf('Results were not saved.\n');
end

% f) the switch statements uses fewer lines than an if stastement, and
% makes it easy to add more options for exaust velocity.


%% PRESENTATION AND FORMATTING [15 MARKS]]
% You do not need to put anything in this section
% Marks will be given for neat, tidy presentation, sensible and informative
% variable names, commenting and remembering to put your name and email address at the top of the
% script
