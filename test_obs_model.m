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
data = data(1:8:490,:); % remove NaNs that show up at the end
t = [data(:,1)];
y = data(:,2:3);
m = length(t);

% Numerically integrate with ode45
tic
options = odeset('AbsTol',1e-2,'RelTol',1e-10);
[~, traj] = ode45(@(t,y) solar_system_force_model(t,y,ephemeris,JDi), ts, X0, options);
toc

%% Plot observations

close all;
nominal_y = zeros(n,2);
for i = 1:n
    nominal_y(i,:) = optical_obs_model(ts(i),traj(i,:),ephemeris,JDi).';
end

dist_to_earth = [];
for i = 1:n
    dist_to_earth(i) = norm(traj(i,1:3) - ephemeris.earth(:,i).');
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


%% Plot trajectory
close all;
figure(1)
pa = plot3(traj(end,1),traj(end,2),traj(end,3),"ro","LineWidth",1);
hold on
plot3(traj(:,1),traj(:,2),traj(:,3),"r--","LineWidth",1)

% Plot Earth and Sun
pe = plot3(ephemeris.earth(1,end),ephemeris.earth(2,end),ephemeris.earth(3,end),"bo");
plot3(ephemeris.earth(1,:),ephemeris.earth(2,:),ephemeris.earth(3,:),"b-")

plot3(ephemeris.sun(1,1),ephemeris.sun(2,1),ephemeris.sun(3,1),"ko")

xlabel("x (km)")
ylabel("y (km)")
zlabel("z (km)")

legend("2024 YR4","2024 YR4 Trajectory","Earth","Earth Trajectory","Sun")
title("2024 YR4 Trajectory [17 Oct 2024 - 23 Dec 2032]")
axis equal

% Animated plot of markers
for k = n%1:10:n
    pa.XData = traj(k,1);
    pa.YData = traj(k,2);
    pa.ZData = traj(k,3);

    pe.XData = ephemeris.earth(1,k);
    pe.YData = ephemeris.earth(2,k);
    pe.ZData = ephemeris.earth(3,k);

    drawnow
end
