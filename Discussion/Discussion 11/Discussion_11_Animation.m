%% Square plate animation: combined view + 3 side views (placeholders for modes)
clear; clc; close all;

%% Plate geometry (ft)
Lx = 1.0;   % length in x
Ly = 1.0;   % length in y (square plate, change if you want a rectangle)

%% Time settings
t_end   = 20;           % total time [s]
nFrames = 500;          % number of frames
t       = linspace(0, t_end, nFrames);

%% Motion definitions (your original)
w1 = 6.549;
w2 = 5.876;
w3 = 6.791;

ScaleFactors = [2.0, 2.0, 10.0];

q1 = ScaleFactors(1)*1/12 * cos(w1*t);
q2 = ScaleFactors(2)*0.0426 * cos(w2*t);
q3 = ScaleFactors(3)*0.0408 * cos(w3*t);

ux    = 1.0 * q1;                    % translation in x [ft]
uy    = 1.0 * q2 + 1.0 * q3;         % translation in y [ft]
theta = -0.094 * q2 + 0.098 * q3;    % rotation [rad]

%% For now: copy the SAME motion into 4 "channels"
% Later you can edit these independently for each mode.

% Main (combined) motion
ux_main    = ux;
uy_main    = uy;
theta_main = theta;

% Mode 1 placeholder
ux_mode1    = 1*q1;
uy_mode1    = 0*q1;
theta_mode1 = 0*q1;

% Mode 2 placeholder
ux_mode2    = 0*q2;
uy_mode2    = 1.0*q2;
theta_mode2 = -0.094*q2;

% Mode 3 placeholder
ux_mode3    = 0*q3;
uy_mode3    = 1.0*q3;
theta_mode3 = 0.098*q3;

%% Figure + layout
fig = figure('Color','w');

% One big column on the left, three stacked on the right
tiledlayout(3, 2, 'TileSpacing','compact', 'Padding','compact');

% --- Main combined motion (big left plot) ---
axMain = nexttile([3 1]);
hold(axMain, 'on'); grid(axMain, 'on');
axis(axMain, 'equal');
axis(axMain, [-1 1 -1 1]);
xlabel(axMain, 'x [ft]');
ylabel(axMain, 'y [ft]');
title(axMain, 'Combined motion');

platePatchMain = patch(axMain, 'XData', [], 'YData', [], ...
    'FaceColor', [0.2 0.6 0.8], ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.5);

% --- Mode 1 (top right) ---
axM1 = nexttile;
hold(axM1, 'on'); grid(axM1, 'on');
axis(axM1, 'equal');
axis(axM1, [-1 1 -1 1]);
title(axM1, 'Mode 1');
platePatch1 = patch(axM1, 'XData', [], 'YData', [], ...
    'FaceColor', [0.8 0.4 0.4], ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.0);

% --- Mode 2 (middle right) ---
axM2 = nexttile;
hold(axM2, 'on'); grid(axM2, 'on');
axis(axM2, 'equal');
axis(axM2, [-1 1 -1 1]);
title(axM2, 'Mode 2');
platePatch2 = patch(axM2, 'XData', [], 'YData', [], ...
    'FaceColor', [0.4 0.8 0.4], ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.0);

% --- Mode 3 (bottom right) ---
axM3 = nexttile;
hold(axM3, 'on'); grid(axM3, 'on');
axis(axM3, 'equal');
axis(axM3, [-1 1 -1 1]);
title(axM3, 'Mode 3');
platePatch3 = patch(axM3, 'XData', [], 'YData', [], ...
    'FaceColor', [0.4 0.4 0.8], ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.0);

%% Animation loop
for k = 1:nFrames
    % Main
    [xMain, yMain] = plate_coords(ux_main(k),    uy_main(k),    theta_main(k),    Lx, Ly);
    % Mode 1
    [x1, y1]       = plate_coords(ux_mode1(k),   uy_mode1(k),   theta_mode1(k),   Lx, Ly);
    % Mode 2
    [x2, y2]       = plate_coords(ux_mode2(k),   uy_mode2(k),   theta_mode2(k),   Lx, Ly);
    % Mode 3
    [x3, y3]       = plate_coords(ux_mode3(k),   uy_mode3(k),   theta_mode3(k),   Lx, Ly);

    % Update patches
    set(platePatchMain, 'XData', xMain, 'YData', yMain);
    set(platePatch1,    'XData', x1,    'YData', y1);
    set(platePatch2,    'XData', x2,    'YData', y2);
    set(platePatch3,    'XData', x3,    'YData', y3);

    % Update title with current values on the main plot
    title(axMain, sprintf('Combined motion  |  t = %.2f s   u_x = %.3f ft   u_y = %.3f ft   \\theta = %.3f rad', ...
                          t(k), ux_main(k), uy_main(k), theta_main(k)));

    drawnow;
    pause(0.05);  % same as before
end

%% ===== Local function(s) below =====
function [xGlobal, yGlobal] = plate_coords(ux, uy, theta, Lx, Ly)
%PLATE_COORDS Return coordinates of a rotated/translated rectangular plate.
%   ux, uy : translations in ft
%   theta  : rotation in rad
%   Lx, Ly : plate side lengths in ft

    % Local coordinates of a rectangle centered at (0,0)
    % Order: 4 corners + first corner repeated to close the polygon
    xLocal = 0.5 * Lx * [-1  1  1 -1 -1];
    yLocal = 0.5 * Ly * [-1 -1  1  1 -1];

    % Rotation matrix
    c = cos(theta);
    s = sin(theta);
    R = [c -s;
         s  c];

    % Rotate
    coordsLocal = [xLocal; yLocal];
    coordsRot   = R * coordsLocal;

    % Translate
    xGlobal = coordsRot(1, :) + ux;
    yGlobal = coordsRot(2, :) + uy;
end
