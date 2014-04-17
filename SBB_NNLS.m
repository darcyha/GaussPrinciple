function x = SBB_NNLS( A, x0, b )
% Computes argmin_x 1/2 \|Ax - b\|^2 

M = 20;
sigma = 0.9;
beta = 1.0;
nu = 0.8;

f = @(A, b, x) ( 0.5 * norm(A * x - b)^2 );
fg = @(A, b, x) ( A' * ( A * x - b ) );

x = x0;
g = fg( A, b, x );
gp = g;
gp( x == 0 & gp > 0 ) = 0;

xOuter = x0;
vOuter = f( A, b, x0 );
gOuter = g;

xVec = {};
vVec = [];
gVec = {};

i = 0;
j = M;
while( norm(gp, 'inf') > 1e-8 )
    Ag = A * gp;
    
    % If norm(Ag) is very small, we will suffer from numerical problems so
    % call it a day and exit early.
    if( norm(Ag) < 1e-9 )
        break
    end
    
    AgAg = Ag' * Ag;
    
    if( mod(i,2) == 0 )
        gg = gp' * gp;
        alpha = gg / AgAg;
    else
        AAg = A' * Ag;
        % I'm not sure about this next line, but the sample code for SBB 
        % included it. Maybe its fixing numerical errors?
        AAg( x == 0 & g > 0 ) = 0; 
        AAgAAg = AAg' * AAg;
        alpha = AgAg / AAgAAg;
    end
    
    %xOld = x;
    x = x - (beta * alpha) * g;
    x( x < 0 ) = 0; % Project to non-negative orthant
    
    gp = g;
    g = fg( A, b, x );
    gp( x == 0 & g > 0 ) = 0;
    
    i = i + 1;
    j = j - 1;
    
    if( j == 0 )
        vCur = f( A, b, x );
        if( vOuter - vCur < sigma * gOuter' * (x - xOuter) )
            beta = beta * nu;
        end
        
        j = M;
        xOuter = x;
        vOuter = vCur;
        gOuter = g;
    end

    %xVec{end+1} = norm(x - xOld);
    %vVec = [vVec f( A, b, x )];
    %gVec{end+1} = norm(gp, inf);
end

return
figure(1);
clf();
semilogy( cell2mat(gVec) );

figure(2);
clf();
semilogy( vVec );

figure(3);
clf();
semilogy( cell2mat(xVec) );
