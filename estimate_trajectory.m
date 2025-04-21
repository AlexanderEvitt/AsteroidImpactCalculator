%% Estimates the trajectory of a MPC object
% Alexander Evitt
clear; clc; close all;

data = unpack_MPC("2024_YR4.txt");
data = data(1:490,:); % remove NaNs that show up at the end

% Times
n = 1000;
ts = linspace(0,365*86400,n).';
JDs = juliandate(2000,1,1,0,0,0) + (ts/86400);

% Pack ephemeris data to struct
ephemeris.dates = JDs;
ephemeris.sun = planetEphemeris(JDs,"SolarSystem","Sun").';
ephemeris.venus = planetEphemeris(JDs,"SolarSystem","Venus").';
ephemeris.earth = planetEphemeris(JDs,"SolarSystem","Earth").';
ephemeris.mars = planetEphemeris(JDs,"SolarSystem","Mars").';
ephemeris.jupiter = planetEphemeris(JDs,"SolarSystem","Jupiter").';

% Ephemeris 2024 YR4 at JD TT = 2460600.5
a = 2.5370741*1.495978707e8;
e = 0.6638562;
i = 3.45180;
argp = 134.64369;
RAAN = 271.41255;
M = 351.06547;
n = 0.24389582;
true_X0 = MPC_to_orbit(a,e,i,RAAN,argp,M,n);


%%
t = [0;data(:,1)];
y = data(:,2:3);
n = 2;
m = length(t);

% Initial conditions
X0 = true_X0 + 0*384*[1000;-2500;500;0.1;-0.2;0.01]; % slightly perturb the known ephemeris to check convergence

%% Carry out batch estimation
clc;
R = 0.01*eye(2);
P = 0.01*diag([10,10,10,0.00001,0.00001,0.00001]);

% Models
F = @(t,y) solar_system_force_model(t,y,ephemeris);
G = @(t,X) optical_obs_model(t,X,ephemeris);


[new_X0,nominal,new_P0] = batch_estimate(F,G,t,y,P,R,3,X0,1);
norm(X0(1:3) - true_X0(1:3))
pos_err = norm(new_X0(1:3) - true_X0(1:3))
vel_err = norm(new_X0(4:6) - true_X0(4:6))

%% Construct residuals
close all;
nominal_y = zeros(m-1,2);
for i = 2:m
    nominal_y(i-1,:) = optical_obs_model(t(i),nominal(:,i),ephemeris).';
end

% Plot observations
figure(1)
subplot(2,1,1)
plot(t(2:end), nominal_y(:,1))
hold on
plot(t(2:end),y(:,1))
xlabel("Time (s)")
ylabel("Right Ascension (deg)")
legend("Nominal trajectory", "Observations")

subplot(2,1,2)
plot(t(2:end), nominal_y(:,2))
hold on
plot(t(2:end),y(:,2))
xlabel("Time (s)")
ylabel("Declination (deg)")
legend("Nominal trajectory", "Observations")