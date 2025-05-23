%% Estimates the trajectory of a MPC object
% Alexander Evitt
clear; clc; close all;

% Times (epoch start to 2032 close approach)
n = 5000;
JDi = juliandate(2024,12,25,0,0,0); % start date
JDf = juliandate(2025,12,23,0,0,0); % end date
span = (JDf - JDi)*86400; % time span in seconds
ts = linspace(0,span,n).';
JDs = JDi + (ts/86400);

% Unpack data
data = unpack_MPC("2024_YR4.txt",JDi);
fin = 417; % 493 for all, 417 for max collision
data = data(1:1:fin,:); % remove NaNs that show up at the end

% Pack ephemeris data to struct
ephemeris = pack_ephemeris(JDs);

% Ephemeris 2024 YR4 at JD TT = 2460669.5
JPL_X0 = [-9.930297342883009e+06; 1.336508161385725e+08; 5.761474336593619e+07;
    -3.682792878261075e+01; 1.006291321097075e+01; 1.963022731732157e+00];
true_X0 = JPL_X0;


%%
offset = 1;
t = [0;data(:,1)];
y = data(:,2:3);
m = length(t);

% Initial conditions
% Perturb by 1 Earth-Moon dist and 0.1 km/s
theta = 33;
phi = -72;
dr = 384000*[cosd(theta)*sind(phi); sind(theta)*sind(phi); cosd(phi)];
theta = -5;
phi = 123;
dv = 0.1*[cosd(theta)*sind(phi); sind(theta)*sind(phi); cosd(phi)];
dX = [dr;dv];
X0 = true_X0 + dX;

%% Carry out batch estimation
clc;
timsec = 0.0041667; % 1 second in degrees
arcsec = 0.00027778;  % 1 arcsecond in degrees
% Nominally precision is 0.01s RA, 0.1" DEC
R = [(0.01*timsec)^2, 0;
    0, (0.1*arcsec)^2];
P = diag(dX.^2);
iters = 10;

% Models
F = @(t,y) solar_system_force_model(t,y,ephemeris,JDi);
G = @(t,X) optical_obs_model(t,X,ephemeris,JDi);

tic
[new_X0,nominal,new_P0,residuals] = batch_estimate(F,G,t,y,P,R,iters,X0,offset);
toc
norm(dX(1:3))
pos_err = norm(new_X0(1:3) - true_X0(1:3))
vel_err = norm(new_X0(4:6) - true_X0(4:6))

%% Construct residuals
close all;
nominal_y = zeros(m-1,2);
for i = 2:m
    nominal_y(i-1,:) = optical_obs_model(t(i),nominal(:,i),ephemeris,JDi).';
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

%% Plot residuals in each iteration

figure(2)
semilogy(1:iters, residuals(:,1))
hold on
semilogy(1:iters, residuals(:,2))
xlabel("Iteration")
ylabel("Sample variance (deg^2)")
legend("Right Ascension", "Declination")