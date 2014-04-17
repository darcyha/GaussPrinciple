function [x,v] = process_collisions( h, kr, x, v, r, Q, Qi, bN, bD )

assert( length(x) == 2*length(r) );
assert( length(bN) == 2*length(bD) );

stepLeft = h;

while( stepLeft > 0 )
    % TODO: Collect the set of objects colliding at the same time! We also need
    % objects in contact with colliding objects too.
    tMin = Inf;
    iMins = {};

    for i=1:length(r)
        idex = [2*i-1,2*i];
        xCur = x(idex);
        vCur = v(idex);
        rCur = r(i);
        for j=1:length(bD)
            jdex = [2*j-1,2*j];
            bNCur = bN(jdex);
            bDCur = bD(j);

            % TODO: Vectorize!
            [flag,t] = intersect_plane( xCur, vCur, rCur, stepLeft, bNCur, bDCur );
            if( flag )
                if( t - tMin < -1e-6 && t > 0 )
                    % New minimum
                    tMin = t;
                    iMins = { [0;i;j] };
                elseif( t - tMin < 1e-6 )
                    iMins{end+1} = [0;i;j];
                end
            end
        end
        
         for j=i+1:length(r)
            jdex = [2*j-1,2*j];
            xOther = x(jdex);
            vOther = v(jdex);
            rOther = r(j);

            % TODO: Vectorize!
            [flag,t] = intersect_particle( xCur, vCur, rCur, stepLeft, xOther, vOther, rOther );
            if( flag )
                if( t - tMin < -1e-6 && t > 0 )
                    % New minimum
                    tMin = t;
                    iMins = { [1;i;j] };
                elseif( t - tMin < 1e-6 )
                    iMins{end+1} = [1;i;j];
                end
            end
        end
    end
    
    if( ~isempty(iMins) && tMin < Inf )
        if( tMin < 1e-8 )
            tMin = 1e-8;
        end
        
        % Advance to the collision
        x = x + tMin*v;
        stepLeft = stepLeft - tMin;
        
        %debug_plot( x )
        
        J = zeros( length(iMins), length(v) );
        c = zeros( length(iMins), 1 );
        
        for k=1:length(iMins)
            iMin = iMins{k}(2);
            jMin = iMins{k}(3);
            idex = [2*iMin-1,2*iMin];
            jdex = [2*jMin-1,2*jMin];

            xCur = x(idex);
            vCur = v(idex);

            if( iMins{k}(1) == 0 )
                bNCur = bN(jdex);
                
                J( k, idex ) = bNCur;
                c( k ) = -kr*bNCur'*vCur;
            else
                xOther = x(jdex);
                vOther = v(jdex);
                n = xCur - xOther;
                n = n / sqrt(n'*n);
                
                J( k, idex ) = n;
                J( k, jdex ) = -n;
                c( k ) = -kr*n'*(vCur - vOther);
            end
        end
        
        k = J \ c; %Jk = c;
        s = Q*(k - v);
        %lambda = lsqnonneg( Qi'*J', s );
        lambda = SBB_NNLS(Qi'*J',zeros( size(J,1),1 ),s);

        v = v + Q' \ (Q \ (J' * lambda));
    else
        x = x + stepLeft*v;
        stepLeft = 0;
    end
end