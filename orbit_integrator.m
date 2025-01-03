%% Euler Integration versus Kepler's Equation Orbit Comparison
% Raymann Singh
% Init
clear, clc, close all
%% Constants
mu = 398600; % km^3/s^2

%% Init Variables
% Change this Variables based on your specific initial position
% and velocity vectors. dt is the integration step size (smaller number
% means more precise orbit, but longer computation time).

% Format in :
%       |1131.34  |            |-5.64305|
% R_i = |-2282.343| km,  V_i = |4.30333 | km/s
%       |6672.423 |            |2.42879 |
%
Ri = [1131.34 -2282.343 6672.423]'; % km
Vi = [-5.64305 4.30333 2.42879]'; % km/s
dt = [0.1 1 10 60 300]; % Integration step sizes
time = 24*60*60; % Time to integrate over (period), seconds

%% Run Functions
[Rf, Vf, orbitPaths, Rmag, Vmag] = integrate(Ri, Vi, mu, dt, time);
[COE] = RV2COE(Ri, Vi, mu);
[nuf, Rnuf, Vnuf] = findTrueNu(COE, mu, time);
[speIntegrator, sphIntegrator] = findSpecifics(Rf, Vf, mu);
[speKepler, sphKepler] = findSpecifics(Rnuf, Vnuf, mu);

%% Comparison
% This is an inital comparison of the specific mechanical energy
% and specific angular momentum of the default orbital elements.
% It is important to note that these are the two elements that are being 
% compared since Kepler's Equation only gives final position and velocity
% vectors which through Euler Integration will only be reached as 
% dt approaches 0.

text = sprintf(['The final position and velocity vectors from the different' ...
    '\nstep sizes are drastic. Below are the final vectors for the given' ...
    '\ntime steps. The graph also shows the 3D plot of the different orbits.' ...
    '\nNote how the larger orbits does not close, that is because how the large' ...
    '\ntime step of 300s. The smaller orbits loop over multiple times because' ...
    '\nof the large amount of iterations with the given period (Period < Iteration Time).' ...
    '\nIn Figures 2 and 3 you can see the Radius over Time and Velocity over Time' ...
    '\ngraphs, note the large differences between dt of 300 and dt of 0.1. In the dt of' ...
    '\n0.1 curve, you are hardly able to discern the amplitude and frequency, this is' ...
    '\nbecause of how large the scale is due to the dt of 300 orbit.']);
disp(text)
disp(table(dt, 'VariableNames',{'dt'}))
disp(table(round(Rf,2), 'VariableNames',{'Final R vector (km)'}))
disp(table(round(Vf,4), 'VariableNames',{'Final V Vector (km/s)'}))
disp(table(round(Rnuf,2),round(Vnuf,4),'VariableNames',{'Kepler Equation R vector', 'Kepler Equation V Vector'}))
disp(table(round(COE(1),2),round(COE(2),2),round(COE(3),2),round(COE(4),2),round(COE(5),2),round(COE(6),2), 'VariableNames',{'a (km)', 'e', 'i (deg)', 'RA (deg)', 'w (deg)', 'nu (deg)'}))
disp(table(round(nuf,2), 'VariableNames',{'Final True Anomaly (deg)'}))
text2 = sprintf(['The R and V vectors for the Kepler Equation is different than the' ...
    '\nintegrator with a small step size. This is because the the integrator is not' ...
    '\nthat efficient so the value is off. The values are also off because of how long the' ...
    '\nintegrator is ran for. If the time to run is longer than the period the orbit will loop' ...
    '\naround itself, but with each completition will grow larger. Below are the Specific' ...
    '\nMechanical Energy values for the integrator and the Kepler Equation. Note how the small' ...
    '\nstep size is close to the kepler equation. Likewise with the Specific Angular Momentum values.']);
disp(text2)
disp(table(round(speIntegrator,3),'VariableNames',{'Specific Mechanical Energy for the Integrator'}))
disp(table(round(speKepler,3), 'VariableNames',{'Specific Mechanical Energy for Kepler'}))
disp(table(round(sphIntegrator,3),'VariableNames',{'Specific Angular Momentum for the Integrator'}))
disp(table(round(sphKepler,3), 'VariableNames',{'Specific Angular Momentum for Kepler'}))
text3 = sprintf(['The reason why the Angular momentum for both kepler and the small step size' ...
    '\nis so close is because both the integrator and kepler equation keeps the orbit in the same plane,' ...
    '\nwhich makes the angular momentum vector always normal to the plane, but the actual value of the vector' ...
    '\ndiffers based off how big. With smaller and smaller step sizes the values will converge to the same kepler' ...
    '\nvalue. Likewise for energy since the only value that matters is the final position and velocity vector.' ...
    '\nIf you properly run the correct amount of iterations, based off the calculated period, they will converge' ...
    '\nto the same kepler value.' ...
    '\n' ...
    '\nScroll up to see the entire comparison.']);
disp(text3)


%% Plotting
hold on
grid on
legend show
for i=1:length(dt)
    % Replace 0 values with Nan
    orbitPaths(orbitPaths==0) = NaN;
    % Plot all orbits
    txt = ['dt = ', num2str(dt(i))];
    plot3(orbitPaths(3*i-2,:), orbitPaths(3*i-1,:), orbitPaths(3*i,:), 'DisplayName', txt)
    title("Orbit in perifocal coordinates")
    % Need this to view in 3D because hold on breaks stuff
    view(3);
end
% Plot Starting Nu
plot3(orbitPaths(1,1), orbitPaths(2,1), orbitPaths(3,1), 'r.', 'MarkerSize', 20, 'DisplayName', 'Starting True Anomaly')
% Plot True Nu
plot3(Rnuf(1), Rnuf(2), Rnuf(3), 'b.', 'MarkerSize', 20, 'DisplayName', 'Final True Anomaly')
hold off

% Plot Radius over time
figure
hold on
grid on
legend show
for i=1:length(dt)
    plotTime = (0:dt(i):time);
    txt = ['dt = ', num2str(dt(i))];
    plot(plotTime, Rmag(i,1:length(plotTime)), 'DisplayName', txt)
    title("Radius over time")
    xlabel("Time (sec)")
    ylabel("Radius (km)")
end
hold off

% Plot Velocity over time
figure
hold on
grid on
legend show
for i=1:length(dt)
    plotTime = (0:dt(i):time);
    txt = ['dt = ', num2str(dt(i))];
    plot(plotTime, Vmag(i,1:length(plotTime)), 'DisplayName', txt)
    title("Velocity over time")
    xlabel("Time (sec)")
    ylabel("Velocity (km/s)")
end
hold off


%% Functions

function [spe, sph]=findSpecifics(R, V, mu)
    [~,n] = size(R);
    r = zeros(1,n);
    v = zeros(1,n);
    spe = zeros(1,n);
    for i=1:n
        r(i) = norm(R(:,i));
        v(i) = norm(V(:,i));
        spe(i) = (v(i)^2)/2 - (mu/r(i));
    end
    sph = cross(R,V);
end

function [COE]=RV2COE(Ri, Vi, mu)
    % COE 1 - a - Semimajor axis
    r = norm(Ri);
    v = norm(Vi);
    epsilon = (v^2)/2 - (mu/r);
    a = -mu/(2*epsilon);
    % COE 2 - e - Eccentricity
    rdotv = dot(Ri,Vi);
    ebar = (1/mu) * ((v^2 - mu/r)*Ri - rdotv*Vi);
    e = norm(ebar);
    % COE 3 - i - Inclination
    hbar = cross(Ri,Vi);
    h = norm(hbar);
    i = acosd(hbar(3)/h);
    % COE 4 - RA - Right Ascension
    K = [0 0 1];
    nbar = cross(K,hbar);
    n = norm(nbar);
    RA = acosd(nbar(1)/n);
    if nbar(2) < 0
        RA = 360 - RA;
    end
    % COE 5 - w - Argument of Perigee
    w = acosd(dot(nbar,ebar)/(n*e));
    if ebar(3) < 0
        w = 360 - w;
    end
    % COE 6 - nu - True Anomaly
    nu = acosd(dot(Ri,ebar)/(r*e));
    if rdotv < 0
        nu = 360 - nu;
    end
    COE = [a e i RA w nu];
end

function [nuf, Rnuf, Vnuf]=findTrueNu(COE, mu, time)
    % Unpack
    a = COE(1);
    e = COE(2);
    i = COE(3);
    RA = COE(4);
    w = COE(5);
    nui = COE(6) * pi / 180; % Deg to Rad

    % Find n
    n = sqrt(mu/a^3); % rad/s

    % Find Ei
    Ei = atan2(sqrt(1-e^2)*sin(nui),(e+cos(nui)));
    % Quadrant Check
    while Ei < 0
        Ei = Ei + 2*pi;
    end

    % Find Mi
    Mi = Ei - e*sin(Ei);
    % Find Mf
    Mf = Mi + n*time; 
    % Quadrant Check
    while Mf > 2*pi
        Mf = Mf - 2*pi;
    end
    
    % Iterate to find Ef
    % Init diff and guess
    % Guess assumes e ~ 0
    diff = 1;
    Ef = Mf;
    while diff > 10e-6
        Efnew = Mf + e*sin(Ef);
        diff = abs(Efnew - Ef);
        Ef = Efnew;
    end

    % Find nuf
    nuf = acos(cos(Ef)-e / (1-e*cos(Ef)));
    nuf = nuf*180/pi; % Degrees
    % Find Rnuf and Vnuf
    p = a*(1-e*e);
    cnu = cosd(nuf);
    snu = sind(nuf);
    Rmag = p/(1+e*cnu);
    Rpqw = zeros(1,3);
    Rpqw(1) = Rmag * cnu;
    Rpqw(2) = Rmag * snu;
    mup = sqrt(mu/p);
    Vpqw = zeros(1,3);
    Vpqw(1) = -mup*snu;
    Vpqw(2) = mup*(e+cnu);

    % Do A matrix and transforms
    cw = cosd(w);
    sw = sind(w);
    ROT3w = [cw -sw 0; sw cw 0; 0 0 1];
    ci = cosd(i);
    si = sind(i);
    ROT1 = [1 0 0; 0 ci -si; 0 si ci];
    cra = cosd(RA);
    sra = sind(RA);
    ROT3ra = [cra -sra 0; sra cra 0; 0 0 1];
    A = ROT3ra*ROT1*ROT3w;
    Rnuf = A*Rpqw';
    Vnuf = A*Vpqw';
end

function [Rf, Vf, orbitPaths, Rmag, Vmag]=integrate(Ri, Vi, mu, dt, time)
    % Preallocate
    orbits = zeros(6, length(dt));
    Rmag = zeros(5,time+1);
    Vmag = zeros(5,time+1);
    for i=1:length(dt)
        % Init state,r, and fxt value for each dt value
        state = [Ri; Vi];
        r = norm(Ri);
        fxt = [Vi; -mu/r^3*Ri];
        orbitPaths(3*i-2:3*i,1) = Ri;
        Rmag(i,1) = norm(Ri);
        Vmag(i,1) = norm(Vi);
        for j=1:time/dt(i)
            % Run Euler Integration
            X = state + fxt*dt(i);
            state = X;
            r = norm(X(1:3));
            fxt = [X(4:6); -mu/r^3*X(1:3)];
            % Put each iteration in an array
            orbitPaths(3*i-2:3*i,j + 1) = X(1:3);
            Rmag(i,j+1) = norm(X(1:3));
            Vmag(i,j+1) = norm(X(4:6));
        end
        % Return final values into an array
        orbits(:,i) = X;
    end
    % Return final values
    Rf = orbits(1:3, :);
    Vf = orbits(4:6, :);
end
