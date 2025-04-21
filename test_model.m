clear; clc; close all;

data = unpack_MPC("2024_YR4.txt");

% Ephemeris 2024 YR4 at JD TT = 2460600.5
a = 2.5370741*1.495978707e8;
e = 0.6638562;
i = 3.45180;
argp = 134.64369;
RAAN = 271.41255;
M = 351.06547;
n = 0.24389582;
X0 = MPC_to_orbit(a,e,i,RAAN,argp,M,n)

% Times
n = 1000;
ts = linspace(0,3*365*86400,n).';
JDs = juliandate(2024,10,17,0,0,0) + (ts/86400);

% Pack ephemeris data to struct
ephemeris.dates = JDs;
ephemeris.sun = planetEphemeris(JDs,"SolarSystem","Sun").';
ephemeris.venus = planetEphemeris(JDs,"SolarSystem","Venus").';
ephemeris.earth = planetEphemeris(JDs,"SolarSystem","Earth").';
ephemeris.mars = planetEphemeris(JDs,"SolarSystem","Mars").';
ephemeris.jupiter = planetEphemeris(JDs,"SolarSystem","Jupiter").';


% Numerically integrate with ode45
options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[~, traj] = ode45(@(t,y) solar_system_force_model(t,y,ephemeris), ts, X0, options);

%% Plot trajectory
close all;
figure(1)
plot3(traj(:,1),traj(:,2),traj(:,3),"r--","LineWidth",1)
hold on

% Plot Earth and Sun
plot3(ephemeris.earth(1,:),ephemeris.earth(2,:),ephemeris.earth(3,:),"b-")
plot3(ephemeris.earth(1,1),ephemeris.earth(2,1),ephemeris.earth(3,1),"bo")

plot3(ephemeris.sun(1,1),ephemeris.sun(2,1),ephemeris.sun(3,1),"ko")

xlabel("x (km)")
ylabel("y (km)")
zlabel("z (km)")

axis equal