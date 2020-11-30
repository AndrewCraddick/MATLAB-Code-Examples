% Part 2: d) Plotting V and E of a circular ring of charge

%% Computing a symbolic expression for V for anywhere in space

syms x y z % phiprime is angle that an elemental dq of the circular charge is located at, x,y and z are arbitrary points in space outside the charge distribution
N = 300; % number of increments to sum, this cannot exceed the resolution of my meshgrid
R = 2; % radius of circle is 2 meters

dphi = 2*pi/N; % discretizing the circle of charge which spans 2pi radians
     
     integrand = 0;
   for phiprime = 0:dphi:2*pi %phiprime ranges from 0 to 2pi in increments of dphi
      
      integrand = integrand + dphi./(sqrt(((x - R.*cos(phiprime) )).^2 + ((y - R.*sin(phiprime) ).^2) + z.^2));
   
   end
   
intgrl = sum(integrand); % uncessary but harmless step that I leave to show that I am using the sum of the above expression for each dphi    

eps0 = 8.854e-12;
kC = 1/(4*pi*eps0);
rhol = 1*10^-9; % linear charge density

Vtot = kC*rhol*R.*intgrl; % symbolic expression for Vtot

%% Graphing V & E in plane perpedicular to ring & passing through center

vlevel = linspace(20, 65, 15);
[Y1, Z1] = meshgrid(-4:.5:4, -4:.5:4);
Vcont1 = subs(Vtot, [x,y,z], {0,Y1,Z1});
contour(Y1,Z1,Vcont1,vlevel)
hold on

xlabel('y - axis [m]')
ylabel('z - axis [m]')
title('Normalized E in a plane perpedicular to a ring of charge (N = 300)')
str = {'Red line is side view', 'of ring of charge'};
text(-1,2,str)

    circle1 =  rectangle('Position',[-2 0 4 .1],'Curvature',[1,1]); % visually displaying ring of charge on plot
    set(circle1,'FaceColor',[1, 0, 0],'EdgeColor',[1, 0, 0]);

g = gradient(-1.*(kC*rhol*R.*intgrl),[x,y,z]); % taking negative gradient of V and finding symbolic equations for Ex, Ey and Ez

 % substituting all the values of the 2D coordinate system for the symbolic x and y variables to get numeric values for Ex and Ey
Ey1 = subs(g(2), [x y z], {0,Y1,Z1});
Ez1 = subs(g(3), [x y z], {0,Y1,Z1});

E1 = sqrt(Ey1.^2 + Ez1.^2); % full numeric magnitude of E in y-z plane

Eynorm1 = Ey1./E1; % This normalizes the electric field lines
Eznorm1 = Ez1./E1;  

quiver(Y1,Z1,Eynorm1,Eznorm1);

hold off

%% Plotting V and E in the plane of ring

vlevel2 = linspace(20, 65, 12);
interval = linspace (-3,3,11);

[X2, Y2] = meshgrid(interval, interval);
Vcont2 = subs(Vtot, [x,y,z], {X2,Y2,0}); 

contour(X2,Y2,Vcont2,vlevel2)
xlabel('x - axis [m]')
ylabel('y - axis [m]')
title('Normalized E in the plane of a ring of charge (N = 300)')

hold on
    circle2 =  rectangle('Position',[-2 -2 4 4],'Curvature',[1,1]); % visually displaying line of charge on plot
    set(circle2,'EdgeColor',[1, 0, 0]);
    circle2. LineWidth = 4;
        
g = gradient(-1.*(kC*rhol*R.*intgrl),[x,y,z]); 
 
Ex2 = subs(g(1), [x y z], {X2,Y2,0});
Ey2 = subs(g(2), [x y z], {X2,Y2,0});

E2 = sqrt(Ex2.^2 + Ey2.^2);

Exnorm2 = Ex2./E2; 
Eynorm2 = Ey2./E2;

quiver(X2,Y2,Exnorm2,Eynorm2,'AutoScalefactor',.4);
