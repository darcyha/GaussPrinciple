function ContactComparison
rng(1);
global x;
global timing;
 
timing = [];

for baseNum = 20:100
    h = .01;
    x = [0:2:2*(baseNum-1) ; ones(1,baseNum)];
    for i=1:baseNum-1
        x = [x, [(0:2:2*(baseNum-i-1))+i;ones(1,baseNum-i)+i*sqrt(3)]];
    end
    v = zeros( size(x) );
    r = ones( 1, size(x,2) )';
    x = reshape( x, prod(size(x)), 1 );
    v = reshape( v, prod(size(v)), 1 );
    m = 1;
    kr = 0.3; % restitution

    bNorms = [ 0; 1 ];
    bDists = [ 0 ];

    opts = solopt();
   
    n = length(x);
    %M = spdiags( m*ones(n,1), 0, n, n );
    Q = spdiags( sqrt(m)*ones(n,1), 0, n, n );
    Qi = spdiags( 1./(sqrt(m)*ones(n,1)), 0, n, n );

    % TODO: Enforce valid velocity for objects already in contact at t0.
    [a, vCur] = compute_accelerations( x, v, r, Q, Qi, bNorms, bDists, true );

    % Compute v1/2 from v0 + h/2 a0
    vCur = vCur + h/2*a;
    v = vCur;
end

plot( timing(2,:), timing(3,:) )

fprintf( '%s\n', mat2str( timing([2,3],:) ) );

path = 'D:\\particles_tk004_%04d.csv';

fileID = fopen( sprintf(path, 0),'w');
fprintf(fileID,'float32 Position[0],float32 Position[1],float32 Position[2],float32 Radius\r\n');
fprintf(fileID,'%g, %g, 0.0, %g\r\n',[reshape( x, 2, [] ); r']);
fclose(fileID);
