function [a,v]=compute_accelerations( x, v, r, Q, Qi, bN, bD, verbose )
global timing;

%Calculate acceleration. Gravity is a known acceleration [0; -1] here.
a = zeros( size(x) );
a( 2:2:end ) = -1;

%Apply contact constraints
assert( length(x) == 2*length(r) );
assert( length(bN) == 2*length(bD) );

rows = [];
cols = [];
vals = [];
nc = 0;

for i=1:length(r)
    idex = [2*i-1,2*i];
    xCur = x(idex);
    rCur = r(i);
    for j=1:length(bD)
        jdex = [2*j-1,2*j];
        bNCur = bN(jdex);
        bDCur = bD(j);
        flag = contact_plane( xCur, rCur, bNCur, bDCur );
        if( flag )
            vCur = v(idex);
            nv = bNCur'*vCur;
            if( nv >= 0 )
                %debug_plot( x );
                %v(idex) = vCur - nv*bNCur;
            end
            
            %v(idex) = vCur * .99;
            
            nc = nc + 1;
            rows = [rows, nc, nc];
            cols = [cols, idex];
            vals = [vals, bNCur'];
        end
    end
    
    for j=i+1:length(r)
        jdex = [2*j-1,2*j];
        xOther = x(jdex);
        rOther = r(j);
        flag = contact_particle( xCur, rCur, xOther, rOther );
        if( flag )
            n = xCur - xOther;
            n = n / sqrt(n'*n);
            
            nc = nc + 1;
            rows = [rows, nc, nc, nc, nc];
            cols = [cols, idex, jdex];
            vals = [vals, n', -n'];
        end
    end
end

J = sparse( rows, cols, vals, nc, size(x,1) );

if( verbose )
    fprintf( '\tJ is %dx%d\n', size(J,1), size(J,2) );
end

if( nc > 0 )
    k = zeros( size(a) ); %Jk = c, c is zeros.
    s = Q*(k - a);
    
    t = cputime;
    %lambda = lsqnonneg( Qi'*J', s );
    lambda = SBB_NNLS(Qi'*J',zeros( size(J,1),1 ),s );
    e = cputime-t;
    
    if( verbose )
        timing = [timing [size(J,1); size(J,2); e]];
        fprintf( '\tElaspsed time: %f\n', e );
    end

    a = a + Q' \ (Q \ (J' * lambda));
end

function flag=contact_plane( x, r, planeNorm, planeDist )
%d = planeNorm'*x - (planeDist+r);
%flag = d <= 1e-4;
flag = planeNorm'*x - (planeDist+r) <= 1e-6;

function flag=contact_particle( x, r, x2, r2 )
d = x - x2;
s = r + r2;
flag = d'*d - s*s <= 1e-6;
