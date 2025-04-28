clear; clc; close all;

% Initial conditions from JPL Horizons
a = 3.798426562794065e+08;
e = 6.641400641819333e-01;
i = 3.452785235470560;
argp = 1.346422366023097e+2;
RAAN = 2.714123388372419e+2;
theta = 3.074395859965221e+2;
[r, v] = orbital_elements_to_rv(a,e,i,RAAN,argp,theta,1.327124400189e11);
X0 = [r;v];

JPL_X0 = [1.490564028447319e8; -1.995902566836032e7; 1.210407097519749e6;
    -9.096151730398509; 3.355424389651756e1; 1.400752719658557e1];

X0 = JPL_X0;

% Times (epoch start to 2032 close approach)
n = 500;
JDi = juliandate(2024,10,17,0,0,0); % start date
JDf = juliandate(2025,05,01,0,0,0); % end date
span = (JDf - JDi)*86400; % time span in seconds
ts = linspace(0,span,n).';
JDs = JDi + (ts/86400);

% Pack ephemeris data to struct
ephemeris = pack_ephemeris(JDs);

% Unpack data
data = unpack_MPC("2024_YR4.txt",JDi);
data = data(1:1:490,:); % remove NaNs that show up at the end
t = [data(:,1)];
y = data(:,2:3);
m = length(t);

% Numerically integrate with ode45
tic
options = odeset('AbsTol',1e-2,'RelTol',1e-10);
[~, traj] = ode45(@(t,y) solar_system_force_model(t,y,ephemeris,JDi), ts, X0, options);
toc

%% Plot observations

nominal_y = zeros(n,2);
traj = traj.';
for i = 1:n
    nominal_y(i,:) = optical_obs_model(ts(i),traj(:,i),ephemeris,JDi).';
end

dist_to_earth = [];
for i = 1:n
    r = traj(1:3,:);
    R = ephemeris.earth(:,i);
    r_ECI = r - R;
    r_hat = r_ECI/norm(r_ECI);
    dist_to_earth(i) = norm(r);
end

% Plot observations
figure(1)
subplot(3,1,1)
plot(ts, nominal_y(:,1))
hold on
plot(t,y(:,1))
xlabel("Time (s)")
ylabel("Right Ascension (deg)")
legend("Nominal trajectory", "Observations")

subplot(3,1,2)
plot(ts, nominal_y(:,2))
hold on
plot(t,y(:,2))
xlabel("Time (s)")
ylabel("Declination (deg)")
legend("Nominal trajectory", "Observations")

% Add plot of distance to Earth
subplot(3,1,3)
plot(ts, dist_to_earth)
xlabel("Time (s)")
ylabel("Distance to Earth")