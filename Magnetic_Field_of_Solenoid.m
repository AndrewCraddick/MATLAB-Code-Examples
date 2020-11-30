%% Lab 7:Numerical Simulation of Magnetic Coils, SOLENOID

%% Constants
clear variables

R = 2; % radius of loops
b = .05; % proportional to the slope of helix as it spirals upwards
lambda = 2*10^-6; % uniform current denisty
I = 2*pi*R*lambda;
u0 = 1.26e-12;
number_of_loops = 100;
length_of_solenoid = 2*pi*number_of_loops*b;
number_of_t_divisions = 3000; % used for t
number_of_XYZ_divisions = 100; % used for XYZ
z_dimension = length_of_solenoid + 5;
t = linspace(0,2*pi*number_of_loops,number_of_t_divisions);
dt = 2*pi*number_of_loops/number_of_t_divisions; % total radians/total increments


%% Solenoid (with N number_of_loops) component calculation

x = linspace(-4,4,number_of_XYZ_divisions);
y = zeros(size(x));
% y = linspace(-4,4,number_of_XYZ_divisions);
z = linspace(-5,z_dimension,number_of_XYZ_divisions);
[X,Y,Z] = meshgrid(x, y, z);
normal_factor = sqrt((R^2) + b^2); % used to normalize vector tangent to solenoid

Bx = zeros(size(X));
By = Bx;
Bz = By;

for t_increment = 1:number_of_t_divisions
    
    Bx_numerator = ((R/normal_factor) .* cos(t(t_increment)) .* (Z - b*t(t_increment)))...
        - ((b/normal_factor) .* (Y - R*sin(t(t_increment))));
    
    By_numerator = ((-R/normal_factor) .* sin(t(t_increment)) .*(Z - b*t(t_increment)))...
        - ((b/normal_factor) .* (X - R*cos(t(t_increment))));
    
    Bz_numerator = ((-R/normal_factor) .* sin(t(t_increment)) .*(Y - R*sin(t(t_increment)))) ...
        - ((R/normal_factor) .* cos(t(t_increment)) .* (X - R*cos(t(t_increment))));
    
    Rmag = sqrt(((X - R*cos(t(t_increment))).^2) + ((Y - R*sin(t(t_increment))).^2) ...
        + ((Z - b*t(t_increment)).^2));
    
    Bx = Bx + (dt*Bx_numerator ./(Rmag.^2));
    By = By + (dt*By_numerator ./(Rmag.^2));
    Bz = Bz + (dt*Bz_numerator ./(Rmag.^2));
    
end

Bx = u0*I.*Bx; % including scientific constants
By = u0*I.*By;
Bz = u0*I.*Bz;


%% 2D Vector Plot

Bmag = sqrt((Bx.^2)+(By.^2)+(Bz.^2));
Bx_norm = Bx./Bmag;
% By_norm = By./Bmag;
By_norm = zeros(size(Bx));
Bz_norm = Bz./Bmag;

index = 1:3:100;

pBx = Bx_norm(index, index, index);
pBy = By_norm(index, index, index);
pBz = Bz_norm(index, index, index);
% pBx = Bx(index, index, index);
% pBy = By(index, index, index);
% pBz = Bz(index, index, index);

pX = X(index, index, index);
pY = Y(index, index, index);
pZ = Z(index, index, index);

% quiver( pY, pZ, pBy, pBz) % 2D vector plot
quiver3(pX, pY, pZ, pBx, pBy, pBz, 'autoscalefactor', 11) % creates 3D vector plot
hold on

% t_plot = linspace(0,2*pi*number_of_loops,3000);
% x_plot = R*cos(t_plot);
% y_plot = R*sin(t_plot);
% z_plot = b*t_plot;
% plot3(y_plot,x_plot,z_plot)
% hold on
% 
% startx = x(index);
% starty = y(index);
% startz = -4*ones(size(pX)); % gives beginning/end/increments to sweep across z axis
% 
% streamline(X,Y,Z,Bx,By,Bz,pX,pY,startz) % creates contour plot
% view(3)

title('Magnetic field in plane perpendicular to solenoid')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xlim([-4 4])
ylim([-4 4])
zlim_plot = z_dimension + 9;

zlim([-4 35])
hold off



