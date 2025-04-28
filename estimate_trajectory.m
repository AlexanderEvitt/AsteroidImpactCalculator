%% Estimates the trajectory of a MPC object
% Alexander Evitt
clear; clc; close all;

% Times (epoch start to 2032 close approach)
n = 5000;
JDi = juliandate(2024,12,25,0,0,0); % start date
JDf = juliandate(2028,12,23,0,0,0); % end date
span = (JDf - JDi)*86400; % time span in seconds
%span = 60*86400;
ts = linspace(0,span,n).';
JDs = JDi + (ts/86400);

% Unpack data
data = unpack_MPC("2024_YR4.txt",JDi);
data = data(1:5:490,:); % remove NaNs that show up at the end

% Pack ephemeris data to struct
ephemeris = pack_ephemeris(JDs);

% Ephemeris 2024 YR4 at JD TT = 2460600.5
JPL_X0 = [1.490564028447319e8; -1.995902566836032e7; 1.210407097519749e6;
    -9.096151730398509; 3.355424389651756e1; 1.400752719658557e1];

true_X0 = JPL_X0;


%%
offset = 1;
t = [0;data(:,1)];
y = data(:,2:3);
m = length(t);

% Initial conditions
dX = 0.01*384*[1000;-2500;500;0.1;-0.2;0.01];
X0 = true_X0 + dX; % slightly perturb the known ephemeris to check convergence

%% Carry out batch estimation
clc;
sigma = 0.000277778;  % 1 arcsecond in radians
R = (sigma^2)*eye(2);
P = 1*diag(dX.^2);

% Models
F = @(t,y) solar_system_force_model(t,y,ephemeris,JDi);
G = @(t,X) optical_obs_model(t,X,ephemeris,JDi);

tic
[new_X0,nominal,new_P0] = batch_estimate(F,G,t,y,P,R,30,X0,offset);
toc
norm(dX(1:3))
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