%% Estimate odds of collision
clear; clc; close all;
file = "converged_hit.mat";

%% First, get trajectory data

% Input x0, either from batch estimation or external source
JPL_X0 = [-9.930297342883009e+06; 1.336508161385725e+08; 5.761474336593619e+07;
    -3.682792878261075e+01; 1.006291321097075e+01; 1.963022731732157e+00];
load(file)
X0 = JPL_X0;

% Times (epoch start to 2032 close approach)
n = 5000;
JDi = juliandate(2024,12,25,0,0,0); % start date
JDf = juliandate(2032,12,23,0,0,0); % end date
span = (JDf - JDi)*86400; % time span in seconds
ts = linspace(0,span,n).';
JDs = JDi + (ts/86400);

% Pack ephemeris data to struct
ephemeris = pack_ephemeris(JDs);

% Numerically integrate with ode45
options = odeset('AbsTol',1e-2,'RelTol',1e-10);
[~, traj] = ode45(@(t,y) solar_system_force_model(t,y,ephemeris,JDi), ts, X0, options);


%% Get closest approach time
tq = linspace(ts(n-50),ts(n),10000);
rs_ast = interpn(ts,traj(:,1:3),tq,"spline");
rs_ear = interpn(ts,ephemeris.earth(:,:).',tq,"spline");

dists = [];
for i = 1:length(tq)
    dists(i) = norm(rs_ast(:,i) - rs_ear(:,i));
end

plot(tq,dists)
[miss_distance,closest_index] = min(dists);
min_time = tq(closest_index);

%% Get state transition at closest approach

% Integrate nominal trajectory
F = @(t,y) solar_system_force_model(t,y,ephemeris,JDi);
show_output = @(tc,y,flag)fprintf('Percent complete = %s\n',mat2str(tc/min_time))*0;
options = odeset("RelTol",1e-3,"AbsTol",1e-2,"OutputFcn",show_output);
Phi0 = eye(6);
Y0 = [new_X0; Phi0(:)];
[~, Yout] = ode45(@(tc, Y) combined_dynamics(tc, Y, F, 6), [0;min_time], Y0, options);

%% Get covariance at closest approach
load(file)
Phi = reshape(Yout(end,6+1:end), 6, 6);
P = Phi*new_P0*Phi.';


%% Get collision probability
clc;
% Interpolate asteroid position at closest approach
r_ast = interpn(ts,traj(:,1:3),min_time,"spline");
v_ast = interpn(ts,traj(:,4:6),min_time,"spline");
[r_ear, v_ear] = planetEphemeris(JDi + (min_time/86400),"SolarSystem","Earth","432t");
r_ear = r_ear.';
v_ear = v_ear.';

% Get in ECI, meters
r = (r_ast - r_ear);
v = (v_ast - v_ear);

% Plug into Foster's method
P_ast = P(1:3,1:3);
P_ear = 1e-6*eye(6); % assume Earth is precisely known, ~1 mm
Pc2D_Foster(r.',v.',P_ast,[0,0,0],[0,0,0],P_ear,6371,1e-12,"circle")