
%% --- Baseline Properties ---------------------------------------------------------------------------------------------------------

t_analysis = 1000;
baseline_x_pos = 0.2;
baseline_wedge_angle = 10;
baseline_M1 = 5;
baseline_h = 20000;
baseline_k = 50;
baseline_cp = 500;
baseline_rho = 5000;
baseline_emissivity = 0.5;
baseline_thickness = 0.01;

%% --- Freestream Conditions ----------------------------------------------------------------------------------------------------------------------------------

[T1,~] = Altitude(baseline_h);
[M2,T2,P2] = TPG_Shocks_Solver(baseline_h,baseline_M1,baseline_wedge_angle);

%% --- Conductivity Variance ------------------------------------------------------------------------------------------------------------

conductivity_calculation_intervals = 100;
lower_power = -5;
upper_power = 5;
power_interval = (upper_power - lower_power)/conductivity_calculation_intervals;
walltemp_vs_conductivity = zeros(conductivity_calculation_intervals+1,2);

for row = 1:conductivity_calculation_intervals+1
    power = lower_power + power_interval * (row-1);  
    disp(power)
    base_10 = 10^power;
    walltemp_vs_conductivity(row,1) = power;
    wall_temp = Variable_HFT(t_analysis,baseline_x_pos,T1,M2,T2,P2,base_10,baseline_cp,baseline_rho,baseline_emissivity,baseline_thickness);
    walltemp_vs_conductivity(row,2) = wall_temp; 
end

plot(walltemp_vs_conductivity(:,1),walltemp_vs_conductivity(:,2))

%% --- Conductivity Variance 2 ------------------------------------------------------------------------------------------------------------

conductivity_calculation_intervals = 100;
lower_power = -5;
upper_power = 5;
power_interval = (upper_power - lower_power)/conductivity_calculation_intervals;
walltemp_vs_conductivity = zeros(conductivity_calculation_intervals+1,2);

for row = 1:conductivity_calculation_intervals+1
    power = lower_power + power_interval * (row-1);  
    disp(power)
    base_10 = 10^power;
    walltemp_vs_conductivity(row,1) = power;
    wall_temp = Variable_HFT(t_analysis,baseline_x_pos,T1,M2,T2,P2,base_10,baseline_cp,baseline_rho,baseline_emissivity,baseline_thickness);
    walltemp_vs_conductivity(row,2) = wall_temp; 
end

conductivity_comparison = length(walltemp_vs_conductivity);

for row = 1:length(walltemp_vs_conductivity)
        
    conductivity_1 = walltemp_vs_conductivity(row,1);
    wall_temp_1 = walltemp_vs_conductivity(row,2);
    wall_temp_2 = wall_temp_1 - 100;
    conductivity_2 = interp1(walltemp_vs_conductivity(:,2),walltemp_vs_conductivity(:,1),wall_temp_2);
    conductivity_comparison(row,1) = 10^conductivity_1;
    conductivity_comparison(row,2) = 10^conductivity_2;

end

plot(conductivity_comparison(row,2),conductivity_comparison(row,1))

%% --- Emissivity Variance ------------------------------------------------------------------------------------------------------------

emissivity_calculation_intervals = 30;
lower_emissivity = 0;
upper_emissivity = 1;
emissivity_interval = (upper_emissivity - lower_emissivity)/emissivity_calculation_intervals;
walltemp_vs_emissivity = zeros(emissivity_calculation_intervals+1,2);

for row = 1:emissivity_calculation_intervals+1
    emissivity = lower_emissivity + emissivity_interval * (row-1);  
    disp(emissivity)
    walltemp_vs_emissivity(row,1) = emissivity;
    wall_temp = Variable_HFT(t_analysis,baseline_x_pos,T1,M2,T2,P2,baseline_k,baseline_cp,baseline_rho,emissivity,baseline_thickness);
    walltemp_vs_emissivity(row,2) = wall_temp; 
end

plot(walltemp_vs_emissivity(:,1),walltemp_vs_emissivity(:,2))

%% --- Heat Capacity Variance ------------------------------------------------------------------------------------------------------------

capacity_calculation_intervals = 100;
lower_capacity = 100;
upper_capacity = 5000;
capacity_interval = (upper_capacity - lower_capacity)/capacity_calculation_intervals;
walltemp_vs_capacity = zeros(capacity_calculation_intervals+1,2);

for row = 1:capacity_calculation_intervals+1
    capacity = lower_capacity + capacity_interval * (row-1);  
    disp(capacity)
    walltemp_vs_capacity(row,1) = capacity;
    wall_temp = Variable_HFT(t_analysis,baseline_x_pos,T1,M2,T2,P2,baseline_k,capacity,baseline_rho,baseline_emissivity,baseline_thickness);
    walltemp_vs_capacity(row,2) = wall_temp; 
end

plot(walltemp_vs_capacity(:,1),walltemp_vs_capacity(:,2))

%% --- Thickness Variance ------------------------------------------------------------------------------------------------------------

thickness_calculation_intervals = 100;
lower_thickness = 0.0001;
upper_thickness = 1;
thickness_interval = (upper_thickness - lower_thickness)/thickness_calculation_intervals;
walltemp_vs_thickness = zeros(thickness_calculation_intervals+1,2);

for row = 1:thickness_calculation_intervals+1
    thickness = lower_thickness + thickness_interval * (row-1);  
    disp(thickness)
    walltemp_vs_thickness(row,1) = thickness;
    wall_temp = Variable_HFT(t_analysis,baseline_x_pos,T1,M2,T2,P2,baseline_k,baseline_cp,baseline_rho,baseline_emissivity,thickness);
    walltemp_vs_thickness(row,2) = wall_temp; 
end

plot(walltemp_vs_thickness(:,1),walltemp_vs_thickness(:,2))

%% --- Density Variance ------------------------------------------------------------------------------------------------------------

density_calculation_intervals = 100;
lower_density = 500;
upper_density = 10000;
density_interval = (upper_density - lower_density)/density_calculation_intervals;
walltemp_vs_density = zeros(density_calculation_intervals+1,2);

for row = 1:density_calculation_intervals+1
    density = lower_density + density_interval * (row-1);  
    disp(density)
    walltemp_vs_density(row,1) = density;
    wall_temp = Variable_HFT(t_analysis,baseline_x_pos,T1,M2,T2,P2,baseline_k,baseline_cp,density,baseline_emissivity,baseline_thickness);
    walltemp_vs_density(row,2) = wall_temp; 
end

plot(walltemp_vs_density(:,1),walltemp_vs_density(:,2))

%% --- Functions ---------------------------------------------------------------------------------------------------------------------------------

function [T_wall] = Variable_HFT(t_analysis,x_pos,T1,M2,T2,P2,k_wall,cp_wall,rho_wall,emis_wall,thick_wall)

% --- Freestream Conditions

% T1 = 215;
% M2 = 4;
% T2 = 300;
% P2 = 3E4;
% k_wall = 21;
% cp_wall = 500;
% rho_wall = 4500;
% emis_wall = 1;
% thickwall = 0.1;
% x_pos = 0.2;

global M_2 %#ok<*GVMIS> 
global T_2
global P_2
global T_1

T_1 = T1;
M_2 = M2;
T_2 = T2;
P_2 = P2;

% --- Transient Factors

t_step = 10;
t_end = t_analysis;

% --- Geometry  

reservoir_y = 0.2;

upper_wall_y = reservoir_y + thick_wall;

plot_y = 0.1+upper_wall_y; %#ok<*NASGU> 

thermalmodel = createpde('thermal','transient');

gdm1 = [3 4 0.1 1 1 0.1 0.0 0.0 reservoir_y reservoir_y]';

gdm2 = [3 4 0.1 1 1 0.1 0.2 0.2 upper_wall_y upper_wall_y]';

test = [gdm1,gdm2];
g = decsg(test);

geometryFromEdges(thermalmodel,g);

% figure; 
% pdegplot(thermalmodel,'EdgeLabels','on','FaceAlpha',0.25,'FaceLabel','on'); 
% axis([0 1.1 -.1 plot_y]);

% --- Material Properties 

global surface_emissivity
surface_emissivity = emis_wall;

k_Ti = @(~,state) Ti6Al4V_Conductivity(state.u);
cp_Ti = @(~,state) Ti6Al4V_HeatCapacity(state.u);

thermalProperties(thermalmodel,'Face',1,'ThermalConductivity',k_Ti, ...
                                        'MassDensity',4500, ...
                                        'SpecificHeat',cp_Ti);

thermalProperties(thermalmodel,'Face',2,'ThermalConductivity',k_wall, ...
                                        'MassDensity',rho_wall, ...
                                        'SpecificHeat',cp_wall);

% --- Initial Conditions and Boundary Conditions 

thermalBC(thermalmodel,'Edge',2, ... 
                       'HeatFlux',@heatflux);

thermalIC(thermalmodel,300);

% --- Mesh Generation

maximum_element_size = .05; 
generateMesh(thermalmodel,'Hmax',maximum_element_size);

%figure;
%pdemesh(thermalmodel)
%axis([0 1.1 -.1 0.3]);
%title('Mesh with Quadratic Triangular Elements')

% --- Solution

tlist = 0:t_step:t_end;
thermalresults = solve(thermalmodel,tlist);
[qx,qy] = evaluateHeatFlux(thermalresults); %#ok<*ASGLU> 

% --- Wall Temperature Interpolation

interpolation_resolution = 0.05;
%step_resolution = 50;

X = 0.1:interpolation_resolution:0.3;
Y = zeros(size(X));
%TABULATED_RESULTS = zeros(length(X),2,1);

for column=1:length(X)
    Y(1,column) = upper_wall_y;
end 

analysis_step = (t_analysis/t_step);
X_transpose = zeros(length(X),1);

for column = 1:length(X)
    X_transpose(column,1) = X(1,column);
end

% --- Interpolation and Storage 

tabulated_results = zeros(length(X),2);
tabulated_results(:,1) = X;
wall_temps = interpolateTemperature(thermalresults,X,Y,analysis_step);
tabulated_results(:,2) = wall_temps;

T_wall = interp1(tabulated_results(:,1),tabulated_results(:,2),x_pos);

end

function Q = heatflux(region,state)

    % --- Importing Flow Conditions. 

    global M_2
    global P_2 %#ok<*GPNES> 
    global T_2
    global T_1
    global surface_emissivity

    % --- Extraction of Temp and Distance. 
    
    T_wall = state.u;
    x = region.x;

    % --- Temps which are NaN are removed. 

    if isnan(T_wall) == true 
        Q = nan;
    else
        Q_aero = Aeroheating(T_wall,x,M_2,T_2,P_2);
        Q_rad = Radiation(T_wall,T_1,surface_emissivity);
        Q = Q_aero - Q_rad;
    end 
  
end 

function[Q] = Aeroheating(T_guess,x,M_2,T_2,P_2)

    % --- Independent Variables 

    T_wall = T_guess;
    x_position = x;

    % --- Parameters for Root Solving 

    Iteration_No = 5;
    Step_no = 5;
    T_lower = 10;
    T_upper = 5000;
    Bisection_Solving_Data = zeros((Step_no+2),2);

    % --- Root Solving for Reference Temperature 

    for iteration=1:Iteration_No

        for row=1:Step_no+2
            
            Bisection_Solving_Data(row,1) = T_lower + ((T_upper-T_lower)/(Step_no+1)) * (row-1);
            T_in = Bisection_Solving_Data(row,1);
            mu = (1.458E-6) * ((T_in^1.5)/(T_in+110.4));
           
            [gamma,cp,~,~,~,k] = ThermoProperties(T_in);
            Pr = (mu*cp)/k;
            recovery_factor = sqrt(Pr);

            T_ref_ratio = 0.45 + 0.55 * (T_wall/T_2) + 0.16 * recovery_factor * ((gamma - 1)/2)*(M_2*M_2);
            T_ref_out = T_ref_ratio * T_2;
            
            delta_T = T_ref_out - T_in;
        
            Bisection_Solving_Data(row,2) = delta_T;

        end 
        
        Bisection_Solving_Data_Reset = Bisection_Solving_Data;
        
        % --- Calculation of Positive Minimum 
        
        for row=1:length(Bisection_Solving_Data)
            if Bisection_Solving_Data(row,2) < 0
                Bisection_Solving_Data(row,2) = nan;
            end
        end
        
        Threshold_Values(1,2) = min(Bisection_Solving_Data(:,2));
        
        Bisection_Solving_Data = Bisection_Solving_Data_Reset;
        
        % --- Calculation of Negative Minimum 
        
        for row=1:length(Bisection_Solving_Data)
            if Bisection_Solving_Data(row,2) > 0
                Bisection_Solving_Data(row,2) = nan;
            end
        end
        
        Threshold_Values(2,2) = max(Bisection_Solving_Data(:,2));
        
        % --- Calculation of Threshold Temperatures 
        
        for threshold = 1:2
            for row = 1:length(Bisection_Solving_Data_Reset)
                if Threshold_Values(threshold,2) == Bisection_Solving_Data_Reset(row,2)
                    Threshold_Values(threshold,1) = Bisection_Solving_Data_Reset(row,1);
                end
            end
        end 
        
        T_ref = interp1(Threshold_Values(:,2),Threshold_Values(:,1),0);
        T_lower = Threshold_Values(1,1);
        T_upper = Threshold_Values(2,1);
        Threshold_Values = zeros(2,2);
    
    end

    % --- Thermodynamic Gas Properties
     
    [gamma,cp,~,R,~,~] = ThermoProperties(T_ref);

    % --- Air Velocity 
    
    [~,~,~,~,a,~] = ThermoProperties(T_2);
    V_e = a*M_2;
    
    % --- Stanton Number and Coefficient of Friction 
   
    P_ref = P_2;
    rho_ref = (P_ref)/(R*T_ref);
    mu_ref = (1.458E-6) * ((T_ref^1.5)/(T_ref+110.4));
    Re_ref_x = (rho_ref*V_e*x_position)/mu_ref;
    St_ref = 0.332 * ((Re_ref_x)^(-0.5)) * Pr^(-2/3);
    
    % --- Adiabatic Wall Temperature 

    T_stag = T_2 * (1 + ((gamma - 1)/2) * M_2^2);
    T_adiabatic = (recovery_factor * (T_stag - T_2)) + T_2;
    
    % --- Flat Plate Aerodynamic Heating 

    Q = rho_ref*V_e*cp*(T_adiabatic-T_wall)*St_ref;

end

function[gamma,cp,cv,R,a,k] = ThermoProperties(T)

    gamma_perf = 1.4;
    cp_perf = 1005;
    T_vibrations = 3055.5556;
    T_ratio = T_vibrations/T;
    x = gamma_perf - 1;
    y = (T_ratio^2)*(exp(T_ratio))/((exp(T_ratio)-1)^2);
    
    gamma = 1 + (x/(1+x*y));
    cp = cp_perf * (1 + (x/gamma_perf)*y);
    cv = cp/gamma;
    R = cp - cv;
    a_squared = R*T*(1+(x/(1+x*y)));
    a = sqrt(a_squared);

    k_0 = 0.0241;
    T_0 = 273;
    S_k = 194; 

    k_k0 = (T/T_0)^1.5 * (T_0+S_k)/(T+S_k);
    k = k_k0 * k_0;

end 

function[Q] = Radiation(T_large,T_small,emissivity)

    sigma = 5.67E-8; % Stefan Boltzmann constant 
    Q = emissivity * sigma * (T_large^4 - T_small^4);

end

function[delta,theta,M_2,T_2,P_2] = TPG_Shocks(h,M_1,T2_T1)
    
    % --- Initial Calculations

    [T_1,P_1] = Altitude(h);
    T_2 = T2_T1 * T_1;
    T1_T2 = T2_T1^-1; 
    [gamma_1,~,~,~,~,~] = ThermoProperties(T_1);
    [gamma_2,~,~,~,~,~] = ThermoProperties(T_2);
    gamma_perf = 1.4;
    thermal_constant = 3065;
    
    % --- Mach Number 

    x = (2/gamma_2)*(1/T2_T1); 
    y = 0.5*gamma_1*M_1^2;
    w = (gamma_perf/(gamma_perf-1))*(1-T2_T1);
    z = (thermal_constant/T_1)*((exp(thermal_constant/T_1)-1)^-1)*((exp(thermal_constant/T_2)-1)^-1);
    
    M2_squared = x*(y+w+z);
    M_2 = sqrt(M2_squared);
    M2_M1 = M_2/M_1;
    
    % --- Pressure Ratio 

    j = (1+gamma_2*M_2^2);
    k = (1+gamma_1*M_1^2);
    P1_P2 = 0.5 * (j - T1_T2 * k + sqrt(((j-T1_T2*k))^2 + 4 * T1_T2));
    P2_P1 = P1_P2^-1;
    P_2 = P2_P1 * P_1;
    
    % --- Density Ratio 

    rho1_rho2 = P1_P2 * T2_T1;
    
    % --- Angles 

    sin_theta_squared = ((gamma_2/gamma_1)*T2_T1*(M2_M1^2)-1)/((rho1_rho2)^2-1); 
    sin_theta = sqrt(sin_theta_squared);
    theta = asin(sin_theta);
    
    cot_delta = tan(theta) * (((gamma_1*M_1^2)/(P2_P1-1))-1);
    delta = acot(cot_delta);
    
    theta = rad2deg(theta);
    delta = rad2deg(delta);

end

function[M_out,T_out,P_out] = TPG_Shocks_Solver(h,M_in,wedge_angle)

% --- Solving Parameters 

T2_T1_Lower = 1;
T2_T1_Upper = 6.5 * M_in;
Intervals = 50000;
[T_in,~] = Altitude(h);

% --- Polar Generation 

TPG_Shock_Polar = zeros(Intervals+1,5);

for row = 1:length(TPG_Shock_Polar)

    TempRatio = T2_T1_Lower + ((T2_T1_Upper-T2_T1_Lower)/Intervals)*(row-1);

    [delta,theta,M_out,~,P_out] = TPG_Shocks(h,M_in,TempRatio);
    TPG_Shock_Polar(row,4) = TempRatio*T_in;
    TPG_Shock_Polar(row,1) = delta;
    TPG_Shock_Polar(row,2) = theta;
    TPG_Shock_Polar(row,3) = M_out;
    TPG_Shock_Polar(row,5) = P_out;

    Percentage_Progress = (row/Intervals)*100;
    if floor(Percentage_Progress/25) == ceil(Percentage_Progress/25)
        disp("Oblique Shock Solving (%): " + Percentage_Progress)
    end

end

% --- Removal of Impossible Flight Scenarios 

for row = 1:length(TPG_Shock_Polar)
    for column = 1:width(TPG_Shock_Polar)
        if isreal(TPG_Shock_Polar(row,column)) == false
            TPG_Shock_Polar(row,column) = nan;
        elseif TPG_Shock_Polar(row,column) <= 0
            TPG_Shock_Polar(row,column) = nan;
        end
    end 
end

for row = 1:length(TPG_Shock_Polar)
    if TPG_Shock_Polar(row,3) <= 1
        TPG_Shock_Polar(row,3) = nan;
    end
end

TPG_Shock_Polar(any(isnan(TPG_Shock_Polar), 2), :) = [];

% --- Plotting & Interpolation 

%plot(TPG_Shock_Polar(:,2),TPG_Shock_Polar(:,1))

M_out = interp1(TPG_Shock_Polar(:,1),TPG_Shock_Polar(:,3),wedge_angle);
T_out = interp1(TPG_Shock_Polar(:,1),TPG_Shock_Polar(:,4),wedge_angle);
P_out = interp1(TPG_Shock_Polar(:,1),TPG_Shock_Polar(:,5),wedge_angle);

end

function[T,P] = Altitude(h)

    if h <= 11000
        T = (15.04 - 0.00649*h)+273.1;
        P = (101.29 * ((T)/288.08)^5.256)*1000;
    elseif h > 11000 && h < 25000
        T = -56.46 + 273.1;
        P = (22.65 * exp(1.73 - 0.000157*h))*1000;
    else
        T = (-131.21 + 0.00299*h)+273.1;
        P = (2.488 * ((T)/216.6)^-11.388)*1000;
    end 
    
end 

function [k] = Ti6Al4V_Conductivity(T)

Ti_Conductivity = zeros(35,2);

Ti_Conductivity(:,1) = [194;214;233;248;273;287;296;322;372;419;477;...
                       519;569;612;656;714;765;812;858;905;924;954;...
                       998;1051;1105;1144;1169;1194;1237;1283;1312;...
                       1380;1430;1478;1532];

Ti_Conductivity(:,2) = [5;5.2;5.4;5.6;6;6.1;6.3;6.6;7.2;7.9;8.7;...
                       9.2;9.7;10.3;11.1;11.9;12.4;13.1;14;14.3;...
                       14.7;16.2;17.4;18.6;20;21;22.5;26.2;32;...
                       42.7;35.2;25.5;25.2;26.2;27];

k = interp1(Ti_Conductivity(:,1),Ti_Conductivity(:,2),T);

end 

function [Cp] = Ti6Al4V_HeatCapacity(T)

Ti_HeatCapacity = zeros(36,2);

Ti_HeatCapacity(:,1) = [194;250;273;300;350;400;450;500;550;600;...
                       650;700;750;800;850;900;950;1000;1050;1100;...
                       1150;1175;1200;1225;1250;1275;1300;1325;1350;...
                       1375;1400;1450;1500;1550;1600;1650];

Ti_HeatCapacity(:,2) = [489;514;523;534;551;565;580;598;611;626;632;648;...
                       658;670;681;693;710;736;771;810;862;893;968;1045;...
                       1200;1510;1257;993;906;788;747;750;746;757;754;757];

Cp = interp1(Ti_HeatCapacity(:,1),Ti_HeatCapacity(:,2),T);
         
end
