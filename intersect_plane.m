function [flag,t]=intersect_plane( x, v, r, h, planeNorm, planeDist )
di = planeNorm'*x - (planeDist+r);
if( di <= 0 )
    % In contact at beginning of step. TODO: Should this use a tolerance?
    flag = true;
    t = 0;
else
    nv = planeNorm'*v;
    if( nv >= 0 )
        % Separating so guaranteed no collision.
        t = 0;
        flag = false;
    else
        t = -di/nv;
        flag = t <= h;
    end
end