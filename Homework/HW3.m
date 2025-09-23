% Problem 2

usto = 0.198;
zeta = 0.3;

Rd_1 = Rd(zeta, 0.1)
uo_1 = Rd_1 * usto

Rd_2 = Rd(zeta, 0.9)
uo_2 = Rd_2 * usto

Rd_3 = Rd(zeta, 3.0)
uo_3 = Rd_3 * usto

% Problem 2
f = 10; g = 386;

dst = 6 * g / (2*pi*f) ^ 2

function [rd] = Rd(zeta, r)
    rd = 1 ./ sqrt((1 - (r).^2).^2 + (2 .* zeta .* r).^2);
end