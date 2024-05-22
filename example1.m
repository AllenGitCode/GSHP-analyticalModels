%% Early Stage Design of BHE

% Li & Lai 2015 (DOI: 10.1016/j.apenergy.2015.04.070) formalize the key
% problem as either:

% a) what is the heat transfer rate (qL) of a GHE as a function of time,
% given a particular temp diff between the circulating fluid and the ground
% (T_f_av - Tgo)

% OR

% b) what is the temp diff between the circulating fluid and the ground
% (T_f_av - Tgo) as a function of time, given a required heat transfer rate
% (qL)

% The following code explores both of these formalizations.
% Formalization (b) can be validated against Li & Lai 2015 figure 4.


%% a) Calculate the qL
% based on a defined flow rate, resistance details, BHE dimensions, etc.,
% and a defined T_f_in:

% steps: 
% 1) set T_f_in, flow rate, and other stuff (ground initial temp, ground props, pipe dimensions, etc)
% 2) get Gb and Rb, and sum them together to get R(t)
% 3) calc qL via rearranged formula: qL = ( Tgo - T_f_in ) / ( H/(2mcp) - R(t) )
% 4) calc T_f_out via: T_f_out = 2( qL*R(t) + Tgo - T_f_in/2 )
% 5) calc T_f_av = (T_f_in + T_f_out) / 2 if needed (optional)

% 1)
T_f_in = 0;                     % Celcius

vol_flow_rate_f = 0.61 / 1000;  % m3/s
rho_f = 1000;                   % kg/m3
cp_f = 4000;                    % J/kg/K
Tgo = 8;                        % Celcius
lambda_g = 2.5;                 % thermal conductivity of ground, W/(m*K)
lambda_b = lambda_g/3.57;
poro = 0.26;
cp_water = 4186;                % J/(kg*K) or J/(kg*C)
rho_water = 1000;               % kg/m3
cp_solid = 880;
rho_solid = 2650;
rhoXcp_g = poro*(rho_water*cp_water) + (1-poro)*(rho_solid*cp_solid); % volumetric heat capacity, J/(m3*K)
alpha_g = lambda_g / rhoXcp_g;  % thermal diffusivity of ground, m3/s
H = 230;                        % meters
rb = 0.25;                      % meters
rp = 0.04 / 2;                  % meters
D = 2 * rp;                     % meters
C = 1;                          % either 1, 2, 3, or 4 (see the 4 cases in Gu & O'Neal 1998)
Ls = C * D;                     % meters
t = 365*24*60*60;               % seconds


% 2)
Gb = Gfunction_FLS(lambda_g, alpha_g, H, rb, H/2, t);
Rb = Rb_equivalent_diameter_single(lambda_b, rb, rp, Ls);
Rt = Rb + Gb;


% 3)
mass_flow_rate_f = vol_flow_rate_f * rho_f;
qL = ( T_f_in - Tgo ) / ( H/(2 * mass_flow_rate_f * cp_f) + Rt);


% 4)
T_f_out = 2 * ( qL * Rt + Tgo - T_f_in/2 );


% 5)
T_f_av = (T_f_in + T_f_out) / 2;


figure
plot([T_f_out T_f_av T_f_in],[0 -H 0], 'o--', 'DisplayName', [num2str(t/24/60/60),' days'])
ylim([-H-10 0])
xlabel('temp (C)'), ylabel('depth (m)')
title('temp-profiles in borehole using analytical method')
legend

clear

%% b) Calculate the (T_f_av - Tgo)
% based on a defined flow rate, resistance details, BHE dimensions, etc.,
% and a required qL:

% steps:
% 1) set qL, flow rate, and other stuff (ground initial temp, ground props, pipe dimensions, etc)
% 2) get Gb and Rb, and sum them together to get R(t) 
% 3) calc T_f_av via basic rearranging of qL*H = (T_f_av - Tgo) / R(t)
% 4) calc T_f_in and T_f_out via:
%       T_f_in  = T_f_av + qL*H / 2*rho*cp*vol_flow_rate
%       T_f_out = T_f_av - qL*H / 2*rho*cp*vol_flow_rate
% 5) calc "temperature response", as Li & Lai call it, via:
%       dT = (T_f_av - Tgo) (optional, for plotting)


% 1)
qL = 1; % W/m

vol_flow_rate_f = 0.61 / 1000;  % m3/s
rho_f = 1000;                   % kg/m3
cp_f = 4000;                    % J/kg/K
Tgo = 8;                        % Celcius
lambda_g = 2.5;                 % thermal conductivity of ground, W/(m*K)
lambda_b = lambda_g/3.57;
poro = 0.26;
cp_water = 4186;                % J/(kg*K) or J/(kg*C)
rho_water = 1000;               % kg/m3
cp_solid = 880;
rho_solid = 2650;
rhoXcp_g = poro*(rho_water*cp_water) + (1-poro)*(rho_solid*cp_solid); % volumetric heat capacity, J/(m3*K)
alpha_g = lambda_g / rhoXcp_g;  % thermal diffusivity of ground, m3/s
alpha_b = alpha_g * (0.45)^2;   % thermal diffusivity of grouting material (used for time scaling?)
H = 230;                        % meters
rb = H/770;                     % meters
rp = 0.04 / 2;                  % meters
D = 2 * rp;                     % meters
C = 1;                          % either 1, 2, 3, or 4 (see the 4 cases in Gu & O'Neal 1998)
Ls = C * D;                     % meters

% times to use
ts = [1e-3 1e-2 0.1 1] * 365*24*60*60; % seconds
for i = 1:length(ts)
    t = ts(i);

    % 2)
    Gb = Gfunction_FLS(lambda_g, alpha_g, H, rb, H/2, t);
    Rb = Rb_equivalent_diameter_single(lambda_b, rb, rp, Ls);
    Rt = Rb + Gb;
    
    
    % 3)
    T_f_av = Tgo + qL * Rt;
    
    
    % 4)
    A = qL * H / (2*rho_f*cp_f*vol_flow_rate_f);
    T_f_in  = T_f_av + A;
    T_f_out = T_f_av - A;
    
    
    % 5)
    dT = T_f_av - Tgo;

    dTs(i) = dT;
    T_f_avs(i) = T_f_av;

end
figure
semilogx(ts/24/60/60/365, T_f_avs-Tgo, 'o--', 'DisplayName','FLS G-function')
xlabel('time, years'), ylabel('dT = Tf - Tgo')
% compare this figure with Fig 4 in Li & Lai 2015's review!

