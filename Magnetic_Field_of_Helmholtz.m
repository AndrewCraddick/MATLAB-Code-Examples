%% Lab 7:Numerical Simulation of Magnetic Coils, HELMHOLTZ & ANTI-HELMHOTZ

%% Constants   -------------------------------------------------
clear variables

R = 2; % radius of loops
b = .001; % proportional to the slope of helix as it spirals upwards
lambda = 2*10^-6; % uniform current denisty
I = 2*pi*R*lambda;
u0 = 1.26e-12;

number_of_loops = 15;
number_of_t_divisions = 2000; % used for t
number_of_XYZ_divisions = 100; % used for XYZ
total_radians = 2*pi*number_of_loops;
length_of_solenoid = 2*pi*number_of_loops*b;
distance_from_xy_plane = R/2; % ENDS R/2 units away from z=0
z__below_offset_term = distance_from_xy_plane + length_of_solenoid;

t = linspace(0,total_radians,number_of_t_divisions);
dt = total_radians/number_of_t_divisions;

%% Component calculation   --------------------------------------

% Equivalent to two short neighboring solenoids (with N number_of_loops) 
x = linspace(-4,4,number_of_XYZ_divisions);
%y = linspace(-4,4,number_of_XYZ_divisions); % uncomment to graph in full
% xzy space
y = zeros(size(x)); % creates graph in xz plane
z = linspace(-5,5,number_of_XYZ_divisions);
[X,Y,Z] = meshgrid(x, y, z);
normal_factor = sqrt((R^2) + b^2); % used to normalize tangent vector

Bx_below = zeros(size(X));
By_below = Bx_below;
Bz_below = By_below;

Bx_above = zeros(size(X));
By_above = Bx_above;
Bz_above = By_above;

for t_increment = 1:number_of_t_divisions
    
    % calculation for coil below  ---------------------------
    
    Bx_numerator_below = ((R/normal_factor) .* cos(t(t_increment)) .* (Z - b*t(t_increment)))...
        - ((b/normal_factor) .* (Y - R*sin(t(t_increment))));
    
    By_numerator_below = ((-R/normal_factor) .* sin(t(t_increment)) .*(Z - b*t(t_increment)))...
        + ((b/normal_factor) .* (X - R*cos(t(t_increment))));
    
    Bz_numerator = ((-R/normal_factor) .* sin(t(t_increment)) .*(Y - R*sin(t(t_increment)))) ...
        - ((R/normal_factor) .* cos(t(t_increment)) .* (X - R*cos(t(t_increment))));
    
    Rmag_below = sqrt(((X - R*cos(t(t_increment))).^2) + ((Y - R*sin(t(t_increment))).^2) ...
        + ((Z - b*t(t_increment) - z__below_offset_term).^2));
    
    Bx_below = Bx_below + (dt*Bx_numerator_below ./(Rmag_below.^2));
    By_below = By_below + (dt*By_numerator_below ./(Rmag_below.^2));
    Bz_below = Bz_below + (dt*Bz_numerator ./(Rmag_below.^2));
    
    % calculation for coil above  -----------------------------
    
    Bx_numerator_above = ((R/normal_factor) .* cos(t(t_increment)) .* (Z - b*t(t_increment)))...
        - ((b/normal_factor) .* (Y - R*sin(t(t_increment))));
    
    By_numerator_above = ((-R/normal_factor) .* sin(t(t_increment)) .*(Z - b*t(t_increment)))...
        - ((b/normal_factor) .* (X - R*cos(t(t_increment))));

    Rmag_above = sqrt(((X - R*cos(t(t_increment))).^2) + ((Y - R*sin(t(t_increment))).^2) ...
        + ((Z - b*t(t_increment) + distance_from_xy_plane).^2));
    
    Bx_above = Bx_above + (dt*Bx_numerator_above ./(Rmag_above.^2));
    By_above = By_above + (dt*By_numerator_above ./(Rmag_above.^2));
    Bz_above = Bz_above + (dt*Bz_numerator ./(Rmag_above.^2));
    
end

Bx = u0*I.*(Bx_below + Bx_above); % net components contributed by each coil
By = u0*I.*(By_below + By_above); % and including scientific constants
Bz = u0*I.*(Bz_below + Bz_above);


%% Plotting   ----------------------------------------

index = 10:2:90;

Bmag = sqrt((Bx.^2)+(By.^2)+(Bz.^2));
Bx_norm = Bx./Bmag;
% By_norm = By./Bmag;
By_norm = zeros(size(Bx));
Bz_norm = Bz./Bmag; % NOTE: this will make the vectors the SAME length but
                    %       NOT a length of ONE when using y = zeros...
                    %       because Bmag is missing Y component

pBx = Bx_norm(index, index, index); 
pBy = By_norm(index, index, index);
pBz = Bz_norm(index, index, index);
% pBx = Bx(index, index, index); % uncomment to use non-normalized vectors
% pBy = By(index, index, index);
% pBz = Bz(index, index, index);

pX = X(index, index, index);
pY = Y(index, index, index);
pZ = Z(index, index, index);

% quiver( pX, pZ, pBx, pBz) % 2D vector plot
quiver3(pX, pY, pZ, pBx, pBy, pBz, 'autoscalefactor', 19) % 3D vector plot
hold on

% Visualizing coils  ------------------------------
t_plot_below = linspace(0,total_radians,2000);
x_plot_below = R*cos(t_plot_below);
y_plot_below = R*sin(t_plot_below);
z_plot_below = b*t_plot_below - z__below_offset_term;
plot3(y_plot_below,x_plot_below,z_plot_below)
% plot(y_plot_below,z_plot_below) % 2D visual of helmholtz

t_plot_above = linspace(0,total_radians,2000);
x_plot_above = R*cos(t_plot_above);
y_plot_above = R*sin(t_plot_above);
z_plot_above = b*t_plot_above + distance_from_xy_plane;
plot3(y_plot_above,x_plot_above,z_plot_above);
% plot(y_plot_above,z_plot_above);


% Streamline (contour) plot -----------------------------------
% startx = x(index);
% starty = y(index);
% startz = -5*ones(size(pX)); % gives beginning/end/increments to sweep 
%                             % across z axis 
% streamline(X,Y,Z,Bx,By,Bz,pX,pY,startz) % creates contour plot
% view(3)

axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
hold off
