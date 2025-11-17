% Moving point load on a simply supported Euler-Bernoulli beam (modal method)
% Author: Chela (for Miguel)
clear; clc;

%% --- beam + load params (SI) ---
E    = 210e9;        % Young's modulus [Pa]
rho  = 7850;         % density [kg/m^3]
b    = 0.05;         % width [m]
h    = 0.008;        % thickness [m]
A    = b*h;          % area [m^2]
I    = b*h^3/12;     % second moment [m^4]
L    = 3.0;          % span [m]

P    = -1e3;          % moving point load [N]
v    = 3.0;          % load speed [m/s]
zeta = 0.1;         % modal damping ratio (same for all modes)
N    = 10;           % number of modes (tune for accuracy)

%% --- time + space discretization for reconstruction ---
t_in  = 0;
t_out = L/v;                 % load leaves at t=L/v
tail  = 10.5*(L/v);           % capture free vibration after exit
Tend  = t_out + tail;
dt    = 5e-4;                % output time step [s]
tout  = (t_in:dt:Tend).';    % output times

nx    = 201;                 % points along the beam for shape plots
xgrid = linspace(0,L,nx).';

%% --- modal properties ---
n     = (1:N).';
beta  = n*pi/L;
omega = (beta.^2) * sqrt(E*I/(rho*A));    % natural circular frequencies [rad/s]

% convenience vectors for vectorized state equations
OMEGA = omega;                   % Nx1
DAMP  = 2*zeta*OMEGA;            % Nx1

%% --- ODE state: y = [q; qdot], size 2N ---
y0 = zeros(2*N,1);

% forcing amplitude factor for each mode (RHS/Mn with Mn = rho*A*L/2)
Ffac = (2*P)/(rho*A*L);          % scalar

% ODE function (vectorized across modes)
odefun = @(t,y) beam_modal_rhs(t,y,OMEGA,DAMP,Ffac,v,L);

% Integrate (use stiff solver if you crank N high)
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[tsim,ysim] = ode45(odefun, tout, y0, opts);

q   = ysim(:, 1:N);              % modal coordinates over time [nt x N]
% rebuild physical deflection w(x,t) = sum phi_n(x)*q_n(t)
Phi = sin(xgrid * (n*pi/L).');   % [nx x N]
W   = Phi * q.';                 % [nx x nt]


%% --- animation of the response (deflected shape + moving load marker) ---
makeVideo = true;                 % set true to save MP4
videoName = 'beam_moving_load.mp4';
fps       = 60;                    % animation frame rate
subsample = max(1, round((1/dt)/fps));   % throttle frames to target FPS

% precompute y-limits for a stable camera
wmax = max(abs(W), [], "all");
yl   = 1.1 * [-wmax, wmax];
if diff(yl) < 1e-12, yl = [-1, 1]*1e-6; end  % fallback if trivial motion

% figure + plot objects
fig = figure('Name','Beam animation'); clf; hold on; box on; grid on;
shapeLine = plot(xgrid, W(:,1), 'LineWidth', 2);
loadPt    = plot(0, 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
xlabel('x [m]'); ylabel('w(x,t) [m]');
title('Moving load on simply supported beam — animation');
xlim([0, L]); ylim(yl);

% optional: "exaggeration" factor for visuals (1 = true scale)
visScale = 1.0;

% video writer
if makeVideo
    vw = VideoWriter(videoName, 'MPEG-4'); vw.FrameRate = fps; open(vw);
end

% animation loop
for k = 1:subsample:numel(tsim)
    t  = tsim(k);
    yk = W(:,k) * visScale;

    % update shape
    set(shapeLine, 'YData', yk);

    % current load position and deflection under the load (if on the span)
    xL = v*t;
    if xL >= 0 && xL <= L
        wL = interp1(xgrid, yk, xL, 'linear', 'extrap');
        set(loadPt, 'XData', xL, 'YData', wL, 'Visible', 'on');
    else
        set(loadPt, 'Visible', 'off');
    end

    drawnow limitrate;

    if makeVideo
        frame = getframe(fig);
        writeVideo(vw, frame);
    end
end

if makeVideo
    close(vw);
    fprintf('Video saved: %s\n', videoName);
end


% %% --- examples: plot midspan response + snapshots ---
% mid = round((nx+1)/2);
% 
% figure; % midspan time history
% plot(tsim, W(mid,:));
% xlabel('Time [s]'); ylabel('w(L/2,t) [m]');
% title('Midspan deflection vs time');
% 
% figure; % deflected shapes at a few times
% tmarks = [0.25*t_out, 0.5*t_out, t_out, Tend];
% hold on;
% for k = 1:numel(tmarks)
%     [~,idx] = min(abs(tsim - tmarks(k)));
%     plot(xgrid, W(:,idx), 'DisplayName', sprintf('t = %.3f s', tsim(idx)));
% end
% xlabel('x [m]'); ylabel('w(x,t) [m]');
% title('Deflected shapes at selected times');
% legend('Location','best'); grid on;
% 
% %% --- helper prints ---
% fprintf('Fundamental freq f1 = %.3f Hz\n', omega(1)/(2*pi));
% % Warn about potential resonance speeds (excitation freq n*pi*v/L ~ omega_n)
% res_v = omega .* L ./ (n*pi); % "critical" speeds by mode
% fprintf('First few resonance-ish speeds [m/s]:\n'); disp(res_v(1:min(6,N)).');

%% ================= helper function =================
function dydt = beam_modal_rhs(t,y,OMEGA,DAMP,Ffac,v,L)
    N = numel(OMEGA);
    q    = y(1:N);
    qdot = y(N+1:end);

    % moving load inside the span only
    if v*t >= 0 && v*t <= L
        s = sin((1:N)'*pi*v*t/L);   % phi_n evaluated at x=vt
        F = Ffac * s;               % Nx1
    else
        F = zeros(N,1);
    end

    qddot = -DAMP.*qdot - (OMEGA.^2).*q + F;

    dydt = [qdot; qddot];
end
