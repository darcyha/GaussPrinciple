function [a,v]=compute_accelerations( x, v, r, Q, Qi, bN, bD )
global opts;

%Calculate acceleration. Gravity is a known acceleration [0; -1] here.
a = zeros( size(x) );
a( 2:2:end ) = -1;

%Apply contact constraints
J = [];
c = [];

assert( length(x) == 2*length(r) );
assert( length(bN) == 2*length(bD) );

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
            
            Jrow = zeros( 1, length(x) );
            Jrow( idex ) = bNCur;
            J = [J; Jrow];
            c = [c; 0];
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
            Jrow = zeros( 1, length(x) );
            Jrow( idex ) = n;
            Jrow( jdex ) = -n;
            J = [J; Jrow];
            c = [c; 0];
        end
    end
end

if( size(J,1) > 0 )
    k = zeros( size(a) ); %Jk = c;
    s = Q*(k - a);
    out = solnls(Qi'*J',s,zeros( size(J,1),1 ),opts);
    lambda = out.x;

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
