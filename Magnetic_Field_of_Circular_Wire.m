%% Lab 7:Numerical Simulation of Magnetic Coils: single circular loop of wire

%% Constants
clear variables % clears workspace
number_of_XYZ_divisions = 75; % used for XYZ
number_of_divisions = 150; % increase for greater accuracy of B field
R = 1; % radius of circular wire is 1 meter
lambda = 2*10^-6; % uniform current denisty
u0 = 1.26e-12;
dphi = 2*pi/number_of_divisions;
dtheta = pi/number_of_divisions;

%% Circular Wire component calculation (Physical Dipole)

phi = linspace(0,2*pi,number_of_divisions);

x = linspace(-3,3,number_of_XYZ_divisions);
y = x;
z = y;
[X,Y,Z] = meshgrid(x, y, z);

Bx = zeros(size(X));
By = Bx;
Bz = By;

for phi_increment = 1:number_of_divisions
    Bx_numerator = Z.*cos(phi(phi_increment));
    By_numerator = Z.*sin(phi(phi_increment));
    Bz_numerator = (sin(phi(phi_increment)).*(R*sin(phi(phi_increment)) - Y)) ...
                  +(cos(phi(phi_increment)).*(R*cos(phi(phi_increment)) - X));
    
    Rx = X - R*cos(phi(phi_increment));
    Ry = Y - R*sin(phi(phi_increment));
    Rz = Z;
    
    Rmag = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
    
    Bx = Bx + (dphi*R*Bx_numerator ./(Rmag.^2));
    By = By + (dphi*R*By_numerator ./(Rmag.^2));
    Bz = Bz + (dphi*R*Bz_numerator ./(Rmag.^2));
    
end

Bx = u0*lambda*.5*R.*Bx; % including scientific constants
By = u0*lambda*.5*R.*By;
Bz = u0*lambda*.5*R.*Bz;

%B = By + Bz;

%% 2D Vector Plot

index = 1:15:75;

pX = X(index, index, index);
pY = Y(index, index, index);
pZ = Z(index, index, index);
pBx = Bx(index, index, index);
pBy = By(index, index, index);
pBz = Bz(index, index, index);

%quiver3(pX, pY, pZ, pBx, pBy, pBz) % creates 3D vector plot
startx = x(index);
starty = y(index);
startz = -3*ones(size(pX)); % gives beginning/end/increments to sweep across z axis

streamline(X,Y,Z,Bx,By,Bz,pX,pY,startz) % creates contour plot
view(3)

%quiver( pY, pZ, pBy, pBz) % 2D vector plot

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Magnetic field lines of circular loop of wire')
xlim([-3 3])
ylim([-3 3])
zlim([-3 3])


% 2D circular loop side view
%circle1 =  rectangle('Position',[-1 0 2 .1],'Curvature',[1,1]);
%    set(circle1,'FaceColor',[1, 0, 0],'EdgeColor',[1, 0, 0]);


% 3D circular loop
circle2 =  rectangle('Position',[-1 -1 2 2],'Curvature',[1,1]);
    set(circle2,'EdgeColor',[1, 0, 0]);
    circle2. LineWidth = 4;
    
% str = {'Red line is current loop'};
% text(-3,2.6,str)

%% Magnetic field of a PURE Dipole

clear variables

u0 = 1.26e-12;
m = 9.11*10^-31;
sci_constants = u0*m/4*pi;
number_of_XYZ_divisions = 30;
x = linspace(-3,3,number_of_XYZ_divisions);
z = x;
y = zeros(size(x));

[X,Y,Z] = meshgrid(x, y, z);

prefactor = 1./(((X.^2) + (Y.^2) + (Z.^2)).^2);

Bxy_prefactor = ((2 + 1./sqrt((X.^2) + (Y.^2) + (Z.^2))));
Bx_x = X.*Z;
By_y = Y.*Z;
Bz_z = 2*(Z.^2) - (X.^2) - (Y.^2);
Bx = sci_constants*prefactor.*Bxy_prefactor.*Bx_x;
By = sci_constants*prefactor.*Bxy_prefactor.*By_y;
Bz = sci_constants*prefactor.*Bz_z;

Bmag = sqrt((Bx.^2)+(By.^2)+(Bz.^2));
Bx_norm = Bx./Bmag;
By_norm = By./Bmag;
Bz_norm = Bz./Bmag;

quiver( X, Z, Bx_norm, Bz_norm)
% quiver3(X, Y, Z, Bx_norm, By_norm, Bz_norm,'autoscalefactor', 9)
hold on

circle1 =  rectangle('Position',[-1 0 2 .1],'Curvature',[1,1]);
    set(circle1,'FaceColor',[1, 0, 0],'EdgeColor',[1, 0, 0]);

% circle2 =  rectangle('Position',[-1 -1 2 2],'Curvature',[1,1]);
%     set(circle2,'EdgeColor',[1, 0, 0]);
%     circle2. LineWidth = 4;

str = {'Line is loop for comparison'};
text(-2.5,2,str)
    
xlabel('x [m]')
ylabel('z [m]')
zlabel('z [m]')
title('Magnetic field lines of Pure Dipole')
xlim([-3 3])
ylim([-3 3])
zlim([-3 3])



