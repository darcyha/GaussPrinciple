function [flag,t]=intersect_particle( x, v, r, h, x2 ,v2, r2 )
xd = x - x2;

c = xd'*xd - (r+r2)^2;
if( c <= 0 )
    % We are already in contact at time 0
    % TODO: Should this use a tolerance?
    flag = true;
    t = 0;
else
    vd = v - v2;
    
    b = 2*xd'*vd;
    if( b >= -1e-10 )
        % Separating so guaranteed no collision.
        flag = false;
        t = 0;
    else
        a = vd'*vd;
        d = b*b - 4*a*c;
        if( d < 0 )
            % No real solutions: Parallel trajectories?
            flag = false;
            t = 0;
        else
            t = [ (-b + sqrt(d))/(2*a), (-b - sqrt(d))/(2*a) ];
            t = min(t(t>=0));
            flag = t <= h;
        end
    end
end