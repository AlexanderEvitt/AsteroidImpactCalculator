clear; clc; close all;

% Initial conditions from JPL Horizons

JPL_X0 = [1.490564028447319e8; -1.995902566836032e7; 1.210407097519749e6;
    -9.096151730398509; 3.355424389651756e1; 1.400752719658557e1];

X0 = JPL_X0;

% Times (epoch start to 2032 close approach)
n = 5000;
JDi = juliandate(2024,10,17,0,0,0); % start date
JDf = juliandate(2032,12,23,0,0,0); % end date
span = (JDf - JDi)*86400; % time span in seconds
ts = linspace(0,span,n).';
JDs = JDi + (ts/86400);

% Pack ephemeris data to struct
ephemeris = pack_ephemeris(JDs);

% Numerically integrate with ode45
tic
options = odeset('AbsTol',1e-2,'RelTol',1e-10);
[~, traj] = ode45(@(t,y) solar_system_force_model(t,y,ephemeris,JDi), ts, X0, options);
toc

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

%% Minimum distance

j = n;
tq = linspace(ts(n-2),ts(n),10000);
rs_ast = interpn(ts,traj(:,1:3),tq);
rs_ear = interpn(ts,ephemeris.earth(:,:).',tq);

dists = [];
for i = 1:length(tq)
    dists(i) = norm(rs_ast(:,i) - rs_ear(:,i));
end

plot(tq,dists)
miss_distance = min(dists)