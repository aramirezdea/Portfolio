%% Ramirez Dynamometer, Incorporated's Dynamometer Modeling Program 

%  This program allows the user to select one of sixteen listed vehicles.
%  Once the make and model of the vehicle have been selected, it will
%  perform calculations for horsepower, torque, and quarter mile times. It
%  will then show these values graphically. In addition, the program will
%  then allow the user to either turbocharge or tune their engine, as
%  appropriate per engine type, and re-perform the aforementioned
%  calculations for comparison. 

%% Created by Alejandro Ramirez de Arellano on April 2nd, 2013; last modified on May 5th, 2013

%% Clean up

clear all          % Clear all variables
clc                % Clear the command window
close all          % Close all windows

%% Parameters for All Engines in All States

% Loads the data components for each of the sixteen vehicles, including
% make/model and corresponding engine displacement, compression rate, redline value, and curb
% weight. Also introduces a few constants of conversion for later
% calculations. 

vehicleDisplacements = load('engineDisplacementValues.txt')'; %varied units
compressionRates = load('compressionRateValues.txt.')'; %unitless
vehicleRedlines = load('redlineValues.txt.')'; %RPM
curbWeights = load('curbWeightValues.txt')'; % Pounds
litersToCubicFeet = 0.0353; %feet^3
ccToCubicFeet = .00003531466667; %feet^3
hpToTorqueConstant = 5252; %unitless
airDensity = 0.074887; %lbm/ft^3
airFuelRatio = 14.7 ; %unitless
lowerHeatingValue = 18999; %btu/lb
thermalEfficiencyConstant = 0.52; %unitless
eQuarterTimeConstant = 5.825; %unitless
btuToJoulesConstant = 1055.05585262; %unitless
forceIt = 0; %unitless
boostNumber = 0; %unitless

% The function "structCreate" takes the loaded vehicle data and creates a structure containing a structure for each vehicle
% that the code will later process in the calculations based off of user vehicle selection.

[ vehicleData ] = structCreate(vehicleDisplacements,compressionRates,vehicleRedlines,curbWeights);

% These facilitate the users vehicle selection. 

    vehicleData(1,1).Make = 'BMW';
    vehicleData(1,2).Make = 'BMW';
    vehicleData(1,3).Make = 'MINI';
    vehicleData(1,4).Make = 'Honda';
    vehicleData(1,5).Make = 'Scion';
    vehicleData(1,6).Make = 'Volkswagen';
    vehicleData(1,7).Make = 'Bugatti';
    vehicleData(1,8).Make = 'Aston Martin';
    vehicleData(1,9).Make = 'Aston Martin';
    vehicleData(1,10).Make = 'Ferrari';
    vehicleData(1,11).Make = 'Ferrari';
    vehicleData(1,12).Make = 'Lamborghini';
    vehicleData(1,13).Make = 'Lamborghini';
    vehicleData(1,14).Make = 'Subaru';
    vehicleData(1,15).Make = 'Mercedes Benz';
    vehicleData(1,16).Make = 'Mercedes Benz';
    
    vehicleData(1,1).Model = '330i';
    vehicleData(1,2).Model = 'm3';
    vehicleData(1,3).Model = 'Cooper';
    vehicleData(1,4).Model = 'Civic Si';
    vehicleData(1,5).Model = 'BR-Z';
    vehicleData(1,6).Model = 'GTI';
    vehicleData(1,7).Model = 'Veyron';
    vehicleData(1,8).Model = 'DB9 Coupe';
    vehicleData(1,9).Model = 'Vanquish';
    vehicleData(1,10).Model = '458 Italia';
    vehicleData(1,11).Model = 'Enzo';
    vehicleData(1,12).Model = 'Gallardo';
    vehicleData(1,13).Model = 'Murcielago';
    vehicleData(1,14).Model = 'WRX STI';
    vehicleData(1,15).Model = 'SLS AMG GT';
    vehicleData(1,16).Model = 'c63 AMG';

% Welcomes and allows the user to select a vehicle and assigns the selected vehicle to a
%variable, which will then be used to reference the corresponding values in
% it's structure within vehicleData.

disp('Welcome to Ramirez Dyno, Incorporated Dynamometer modeling program')
disp('We recommend you enlarge your command window for viewing ease.');
disp(' ');
disp('Please select your vehicle.')
disp(' ')
for x=1:length(vehicleData)
    disp([num2str(x),' = ' , vehicleData(1,x).Make,' ', vehicleData(1,x).Model]);
end
disp(' ');

wantedMake = input('Which vehicle make would you like to choose?')';


% This user defined function serves to convert the engine displacements in 
% the loaded files, which are in different units (cc or liters), to cubic
% feet. 

[ vehicleData ] = convertDisplacement( vehicleData,litersToCubicFeet,...
    ccToCubicFeet );
%% ENGINES IN NATURAL STATE

%% User Input for Engine in Natural State

% This loop provides the user with a chance to re-enter their selection if they 
%choose an invalid vehicle number 

count = 0;
    while count == 0
        if wantedMake > 16
            disp('ERROR: The answer you provided was invalid. Please enter a number between 1 and 16')
            wantedMake = input('Which vehicle make would you like to choose?')';
        else
                count = count +1;
        end
    end
    
% Once the vehicle is properly selected, the program proceeds. 
    
% Sets the RPM band for the vehicle based on the corresponding Redline
% value (taken from the structure's data)

engineSpeed = 1:1:vehicleData(1,wantedMake).redline;

%% Calculations for Engine in Natural State

% The function pureDyno now takes all of these parameters and runs 
% the power and torque calculations

[powerOutputHP,torqueOutput,maxPowerOutputHP,whereHP,maxTorqueOutput,...
    whereTQ,eQuarterTime,maxPowerRPM ,maxTorqueRPM]...
    = pureDyno(wantedMake,vehicleData,engineSpeed,litersToCubicFeet,...
    ccToCubicFeet,hpToTorqueConstant,airDensity, airFuelRatio,...
    lowerHeatingValue,thermalEfficiencyConstant,eQuarterTimeConstant,...
    btuToJoulesConstant);

%% Outputs for Engine in Natural State

% Output values for Engine in it's natural state, including a plot of
% horsepower and torque over the entire RPM band of the vehicle, annotated
% with the respective maximums. 

figure;
hold on
title([vehicleData(1,wantedMake).Make, ' ', vehicleData(1,wantedMake).Model]);
[AX,H1,H2] = plotyy(engineSpeed,powerOutputHP,engineSpeed,torqueOutput);
xlabel('RPM');
ylabel('Horsepower');
set(get(AX(2),'Ylabel'),'String','Torque (ft-lbs)'); 

% Annotates max horsepower and torque on the plot.  

legend(['Maximum HP = ',num2str(maxPowerOutputHP),' HP at ',num2str(maxPowerRPM),' RPM'],['Maximum Torque = ',num2str(maxTorqueOutput),' foot-pounds at ',num2str(maxTorqueRPM),' RPM'],'location','Northwest');
hold off

% Displays the natural state engine calculations for the user 

disp(' '); %for spacing
disp(['You have chosen the ',vehicleData(1,wantedMake).Make,' ', vehicleData(1,wantedMake).Model]);
disp(['Runs a quarter mile sprint in ' ,num2str(eQuarterTime),' seconds']);
disp(['Maximum power output of ' ,num2str(maxPowerOutputHP),' horsepower at ',num2str(maxPowerRPM),' RPM']);
disp(['Maximum torque output of ' ,num2str(maxTorqueOutput),' foot-pounds at ',num2str(maxTorqueRPM),' RPM']);
disp(' '); %for spacing


%% MODIFIED ENGINE STATES  

%% User input for Turbocharged Engine State

% This allows the users to turbocharge their engine, assuming their vehicle
% selection is not already turbocharged or is compatible with a
% turbocharger. The while loop will ensure that if the user types a number other 
% than 1 or 2, they are re-prompted after an error message. 

if wantedMake == 1 || wantedMake == 3 || wantedMake == 5  %all the others that dont already have turbos
    disp([num2str(1),'=' , 'Yes, Turbocharge it!']);
    disp(' ');
    disp([num2str(2),'=' , 'No thanks, I prefer reliability to power.']);
    disp(' ');
    count = 0;
    forceIt = input('Would you like to turbocharge your engine?');
        while count == 0
            if forceIt ~= 1 && forceIt ~= 2
            disp('ERROR: The answer you provided was invalid. Please enter a 1 or a 2')
            forceIt = input('Would you like to turbocharge your engine?');
            else
                count = count +1;
            end
        end
    
% The user now selects how much boost they would like from the turbo. Once 
% selected, the usual calculations will be performed for comparison. The
% while loop will ensure that if the user types a number other than 1,2, or 3, 
% they are re-prompted after an error message. 

    if forceIt == 1;
        disp([num2str(1),'=' , ' 13 psi (conservative novice)']);
        disp(' ');
        disp([num2str(2),'=' , ' 14 psi (boy racer)']);
        disp(' ');
        disp([num2str(3),'=' , ' 15 psi (Andretti status)']);
        disp(' ');
        
        count = 0;
        boostNumber = input('How much boost would you like the turbo to deliver?'); % unitless
        while count == 0
            if boostNumber > 3
            disp('The answer you provided was invalid. Please enter a 1,2,or 3')
            boostNumber = input('How much boost would you like the turbo to deliver?'); %unitless
            else
                count = count +1;
            end
        end

%% Calculations for Turbocharged Engine State

thermalEfficiencyConstant = 0.54; %Turbo-ing the engine will increase efficiency
        
        if boostNumber == 1;
            boostValue = 13; %psi
            effectiveCompressionRate = sqrt((boostValue+14.7)/14.7)*vehicleData(1,wantedMake).compressionRatio; % unitless
    % The function performs the usual power, torque, and quarter mile
    % calculations.
            [powerOutputHP,torqueOutput,maxPowerOutputHP,whereHP,maxTorqueOutput,whereTQ,eQuarterTime,maxPowerRPM ,maxTorqueRPM] = turboDyno(boostValue,effectiveCompressionRate,wantedMake,vehicleData,engineSpeed,litersToCubicFeet, ccToCubicFeet,hpToTorqueConstant,airDensity, airFuelRatio,lowerHeatingValue,thermalEfficiencyConstant,eQuarterTimeConstant,btuToJoulesConstant);
        
        elseif boostNumber ==2 
            boostValue = 14; %psi 
            effectiveCompressionRate = sqrt((boostValue+14.7)/14.7)*vehicleData(1,wantedMake).compressionRatio; % unitless
    % The function performs the usual power, torque, and quarter mile
    % calculations.
            [powerOutputHP,torqueOutput,maxPowerOutputHP,whereHP,maxTorqueOutput,whereTQ,eQuarterTime,maxPowerRPM ,maxTorqueRPM] = turboDyno(boostValue,effectiveCompressionRate,wantedMake,vehicleData,engineSpeed,litersToCubicFeet, ccToCubicFeet,hpToTorqueConstant,airDensity, airFuelRatio,lowerHeatingValue,thermalEfficiencyConstant,eQuarterTimeConstant,btuToJoulesConstant);
       
        else 
            boostValue = 15; %psi
            effectiveCompressionRate = sqrt((boostValue+14.7)/14.7)*vehicleData(1,wantedMake).compressionRatio; %unitless
    % The function performs the usual power, torque, and quarter mile
    % calculations.
            [powerOutputHP,torqueOutput,maxPowerOutputHP,whereHP,maxTorqueOutput,whereTQ,eQuarterTime,maxPowerRPM ,maxTorqueRPM] = turboDyno(boostValue,effectiveCompressionRate,wantedMake,vehicleData,engineSpeed,litersToCubicFeet, ccToCubicFeet,hpToTorqueConstant,airDensity, airFuelRatio,lowerHeatingValue,thermalEfficiencyConstant,eQuarterTimeConstant,btuToJoulesConstant);
        end
%% Outputs for Turbocharged engine state

% Output values for engine in it's turbocharged state, including a plot of
% horsepower and torque over the entire RPM band of the vehicle, annotated
% with the respective maximums. 

figure;
hold on
title(['Turbocharged ',vehicleData(1,wantedMake).Make, ' ', vehicleData(1,wantedMake).Model]);
[AX,H1,H2] = plotyy(engineSpeed,powerOutputHP,engineSpeed,torqueOutput);
xlabel('RPM');
ylabel('Horsepower');
set(get(AX(2),'Ylabel'),'String','Torque (ft-lbs)'); 

% Annotate max hp and torque on the graph 
legend(['Maximum HP = ',num2str(maxPowerOutputHP),' HP at ',num2str(maxPowerRPM),' RPM'],['Maximum Torque =',num2str(maxTorqueOutput),' foot-pounds at ',num2str(maxTorqueRPM),' RPM'],'location','Northwest');
hold off
% Displays the turbocharged engine calculations for the user 
disp(' ');
if boostNumber == 1 || boostNumber == 2 || boostNumber == 3
    disp(['You have chosen to Turbocharge your ',vehicleData(1,wantedMake).Make,' ', vehicleData(1,wantedMake).Model, ' with a boost of ',num2str(boostValue),' psi']);
    disp(['Now runs a quarter mile sprint in ' ,num2str(eQuarterTime),' seconds']);
    disp(['New maximum power output of ' ,num2str(maxPowerOutputHP),' horsepower at ',num2str(maxPowerRPM),' RPM']);
    disp(['New maximum torque output of ' ,num2str(maxTorqueOutput),' foot-pounds at ',num2str(maxTorqueRPM),' RPM']);
    disp(' ');
end
    else 
        disp(' ');
        disp('You were lame and decided to keep your car as is, so the previous power and torque numbers still apply. Consider adding more excitement into your lifestyle.');
    end
    
%% User Inputs for Tuned Engine State

% The user who selected a non-Turbo compatible car is instead able to select
% whether or not they want to tune their engine. If they choose to tune,
% the code will run the usual calculations for comparison. The
% while loop will ensure that if the user types a number other than 1 or a 2, 
% they are re-prompted after an error message. 

elseif wantedMake ~= 1 && wantedMake ~= 3 && wantedMake ~= 5
    disp('Your cars engine already incorporates forced induction');
    disp(' ');
    disp([num2str(1),'=' , 'Yes, I crave speed!']);
    disp(' ');
    disp([num2str(2),'=' , 'No thanks, Im too scared to void my warranty.']);
    disp(' ');
    wannaTune = input('Would you like to add a performance tune to your engine?');
    count = 0;
        while count == 0
            if wannaTune > 2
            disp('ERROR: The answer you provided was invalid. Please enter a 1 or a 2')
            wannaTune = input('Would you like to add a performance tune to your engine?');
            else
                count = count +1;
            end
        end
        
%% Calculations for Tuned Engine State

% Once a Tune is selected, the usual calculations will be performed for comparison.    
    
    if wannaTune == 1
        thermalEfficiencyConstant= 0.55;
        % The function calculates the power, torque, and quarter mile time
        % over the rpm band of the tune engine. 
        [powerOutputHP,torqueOutput,maxPowerOutputHP,whereHP,maxTorqueOutput,whereTQ,eQuarterTime,maxPowerRPM ,maxTorqueRPM] = pureDyno(wantedMake,vehicleData,engineSpeed,litersToCubicFeet, ccToCubicFeet,hpToTorqueConstant,airDensity, airFuelRatio,lowerHeatingValue,thermalEfficiencyConstant,eQuarterTimeConstant,btuToJoulesConstant);    
        
%% Outputs for Tuned Engine State

% Output values for Engine in it's tuned state, including a plot of
% horsepower and torque over the entire RPM band of the vehicle, annotated
% with the respective maximums. 

    figure;
    hold on
    title(['Tuned ',vehicleData(1,wantedMake).Make, ' ', vehicleData(1,wantedMake).Model]);
    [AX,H1,H2] = plotyy(engineSpeed,powerOutputHP,engineSpeed,torqueOutput);
    xlabel('RPM');
    ylabel('Horsepower');
    set(get(AX(2),'Ylabel'),'String','Torque (ft-lbs)'); 
% Annotates max hp and torque on the plot. 
    legend(['Maximum HP = ',num2str(maxPowerOutputHP),' HP at ',num2str(maxPowerRPM), ' RPM'],['Maximum Torque = ',num2str(maxTorqueOutput),' foot-pounds at ',num2str(maxTorqueRPM),' RPM'],'location','Northwest');
    hold off
    disp(' ');
    disp(['You have chosen to Tune your ',vehicleData(1,wantedMake).Make,' ', vehicleData(1,wantedMake).Model]);
    disp(['Now runs a quarter mile sprint in ' ,num2str(eQuarterTime),' seconds']);
    disp(['New maximum Power output of ' ,num2str(maxPowerOutputHP),' horsepower at ',num2str(maxPowerRPM),' RPM']);
    disp(['New maximum Torque output of ' ,num2str(maxTorqueOutput),' foot-pounds at ',num2str(maxTorqueRPM),' RPM']);
    disp(' ');
    else % If you decided not to Tune your car, you are shunned and dismissed
    disp(' ');
    disp('You were lame and decided to keep your car as is, so the previous power and torque numbers still apply. Consider adding more excitement into your lifestyle.');
    end
end
%% Farewell Message

disp(' ')
disp('Thank you for choosing Ramirez Dyno, Inc. We look forward to doing business in the future.')