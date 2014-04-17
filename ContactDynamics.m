function ContactDynamics
rng(1);
global x;
x = [0.1;.8];
v = [0;0];
%x = [.051;1;.2;1;0.2;0.05;0.1;0.05];
%v = [.2;0.001;.01;0;0;0;-.1;0];
%x = [0.1;0.05];
%v = [-.025;0];
%x = [ 0.05;0.05;  0.15;0.05;  0.1;((1+sqrt(3))/2)*0.1;  0.25;0.05;  0.2;((1+sqrt(3))/2)*0.1;  0.15; ((1+2*sqrt(3))/2)*0.1];
%v = zeros( size(x) );
%x = [.051;.91;  0.05;0.55;  0.45;0.6];
%v = [.2;0.001;  0.5;0;  -.5;0];
%x = [x; rand(10,1).*repmat([.5;1],5,1)+repmat([0.05;.2],5,1)];
%v = [v; zeros(10,1)];
h = .01;
r = 0.01*ones( size(x,1)/2, 1 );
%r = 0.05 - 0.04*rand( size(x,1)/2, 1 );
m = 1;
%kr = 0.01; % restitution
kr = 0.3; % restitution

%n = length(x);
%M = spdiags( m*ones(n,1), 0, n, n );
%Q = spdiags( sqrt(m)*ones(n,1), 0, n, n );
%Qi = spdiags( 1./(sqrt(m)*ones(n,1)), 0, n, n );

bNorms = [ ...
    %cos(pi/2 + pi/90) ; sin(pi/2 + pi/90) ; ...
    0; 1 ; ...
    %cos(pi/2 + pi/6) ; sin(pi/2 + pi/6) ; ...
    -1 ; 0 ; ...
    1 ; 0 ];

bDists = [ 0 ; -5 ; 0 ];

global opts; 
opts = solopt();

n = length(x);
%M = spdiags( m*ones(n,1), 0, n, n );
Q = spdiags( sqrt(m)*ones(n,1), 0, n, n );
Qi = spdiags( 1./(sqrt(m)*ones(n,1)), 0, n, n );

% TODO: Enforce valid velocity for objects already in contact at t0.
[a, vCur] = compute_accelerations( x, v, r, Q, Qi, bNorms, bDists, false );

% Compute v1/2 from v0 + h/2 a0
vCur = vCur + h/2*a;
v = vCur;

path = 'D:\\particles_tk005_%04d.csv';
fileID = fopen( sprintf(path, 0),'w');
fprintf(fileID,'float32 Position[0],float32 Position[1],float32 Position[2],float32 Radius\r\n');
fprintf(fileID,'%g, %g, 0.0, %g\r\n',[reshape( x, 2, [] ); r']);
fclose(fileID);

for i=1:5000
    if( mod(i, 100) == 0 )
        fprintf( '\rFrame %04d', i );
    end
    if( mod(i, 2) == 0 )
        xNew = zeros(2,1);
        xNew(1:2:end-1) = rand(1,1)*4.9+0.05;
        xNew(2:2:end) = rand(1,1)*2+4;
        x = [x;xNew];
        v = [v;zeros(2,1)];
        %r = [r;0.05*ones(1,1)];
        r = [r;0.03+0.005*rand(1,1)];
    end
    
    n = length(x);
    %M = spdiags( m*ones(n,1), 0, n, n );
    m = ones( size(x) );
    %m(1:2:end) = r/0.05;
    %m(2:2:end) = r/0.05;
    Q = spdiags( sqrt(m), 0, n, n );
    Qi = spdiags( 1./(sqrt(m)), 0, n, n );
    
    xCur = x;
    vCur = v;
    
    [xCur, vCur] = process_collisions( h, kr, xCur, vCur, r, Q, Qi, bNorms, bDists );
    [a, vCur] = compute_accelerations( xCur, vCur, r, Q, Qi, bNorms, bDists, false );
    
    vCur = vCur + h*a;

    x = xCur;
    v = vCur;
    
    fileID = fopen( sprintf(path, i),'w');
    fprintf(fileID,'float32 Position[0],float32 Position[1],float32 Position[2],float32 Radius\r\n');
    fprintf(fileID,'%g, %g, 0.0, %g\r\n',[reshape( x, 2, [] ); r']);
    fclose(fileID);
end
