%% Plotting Electric field and Potential of non-uniformly charged sphere

% Constants
clear variables
number_of_divisions = 75;
R = 1; % radius of sphere is 1 meter
sigma = 2*10^-6; % uniform surface charge denisty
eps0 = 8.854e-12;
k = 1/(4*pi*eps0);
dphi = 2*pi/number_of_divisions; % discretizing the interval over which I sum
dtheta = pi/number_of_divisions;

%% Potential V Computation

phi = linspace(0,2*pi,number_of_divisions);
theta = linspace(0,pi,number_of_divisions);
y = linspace(-8,8,number_of_divisions/2);
z = y;
x = linspace(0,8,number_of_divisions/2); %% space extends equally in all directions

[X,Y,Z] = meshgrid(x, y, z);
V = zeros(size(X));
for phi_increment = 1:number_of_divisions
    for theta_increment = 1:number_of_divisions
        Rx = X - R*sin(theta(theta_increment)).*cos(phi(phi_increment));
        Ry = Y - R*sin(theta(theta_increment)).*sin(phi(phi_increment));
        Rz = Z - R*cos(theta(theta_increment));
        
            cos_theta_sin_theta = cos(theta(theta_increment)).*sin(theta(theta_increment));
            
            Rmag = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
            V = V + (cos_theta_sin_theta * dtheta * dphi/Rmag);
    end
end

V = eps0*k*sigma*(R^2)*V;

%% Gradient Computation

dx = x(2)-x(1);
dy = dx;
dz = dx;
[Ex, Ey, Ez] = gradient(-V, dx, dy, dz);

%% Visualization

levels = 7;
isovalues = [-70 -40 -20 0 20 40 70];
for contour = 1:levels
    [faces1, verts1, colors1] = isosurface(X, Y, Z, V/10^-9, isovalues(contour), V/10^-9); 
    p_sphere = patch('Vertices', verts1, 'Faces', faces1, 'FaceVertexCData',... 
              colors1, 'FaceColor','interp','EdgeColor','interp');
    isonormals(X, Y, Z, V/10^-9, p_sphere)
    p_sphere.EdgeColor = 'none';
    p_sphere.FaceAlpha = 0.3;
end
daspect([1 1 1])
view(3); 
axis vis3d 
camlight 
lighting gouraud 
h = colorbar;
h.Label.String = '[nV]';
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Spherical conductor with non-uniform charge distribution')

%% Plot of a physical electric dipole for comparison
% Here is the plot of a dipole and the sphere overlayed at the approximate
% value p should be to make the two equivalent

d = 1; % charges are 1 meter apart
q_dipole = 1.6*10^-(19);
p = 6.1*q_dipole*d;
V_dipole = (p/4*pi*eps0).*Z./((X.^2 + Y.^2 + Z.^2).^(3/2));

levels_dipole = 5;
isovalues_dipole = [-60 -30 0 30 60]; 

for contour = 1:levels_dipole
    [faces2, verts2, colors2] = isosurface(X, Y, Z, V_dipole/10^-32, isovalues_dipole(contour), V_dipole/10^-32); 
    p_dipole = patch('Vertices', verts2, 'Faces', faces2, 'FaceVertexCData',... 
              colors2, 'FaceColor','red','EdgeColor','red');
    
    isonormals(X, Y, Z, V_dipole, p_dipole)  
    p_dipole.EdgeColor = 'none';
    p_dipole.FaceAlpha = 0.5; 
end
hold on

levels = 5;
isovalues = [-60 -30 0 30 60];
for contour = 1:levels
    [faces3, verts3, colors3] = isosurface(X, Y, Z, V/10^-9, isovalues(contour), V/10^-9); 
    p_sphere = patch('Vertices', verts3, 'Faces', faces3, 'FaceVertexCData',... 
              colors3, 'FaceColor','blue','EdgeColor','blue');
    isonormals(X, Y, Z, V/10^-9, p_sphere)
    p_sphere.EdgeColor = 'none';
    p_sphere.FaceAlpha = 0.5;
end


daspect([1 1 1]) 
view(3); 
axis vis3d 
camlight 
lighting gouraud 
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Potentials are approximately equal at p = 9.76 \cdot 10^{4}')

hold off
legend([p_sphere(1) p_dipole(1)],'Potential of spherical shell','Potential of dipole')
