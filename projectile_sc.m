function [value,isterminal,direction] = projectile_sc(t,y)

value(1) = y(3)-0;     % Detect height = 0
isterminal(1) = 1;   % Stop the integration
direction(1) = -1;   % Negative direction only