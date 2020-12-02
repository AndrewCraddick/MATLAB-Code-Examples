%% Plotting Electric field and Potential of uniformly charged sphere

% Constants -----------------------------------------------
clear variables % clears prexisting data from workspace
number_of_divisions = 50;
R = 1; % radius of sphere is 1 meter, not realistic but who cares
sigma = 2*10^-6; % surface charge denisty
eps0 = 8.854e-12;
k = 1/(4*pi*eps0);
dphi = 2*pi/number_of_divisions; % discretizing the interval over which I sum
dtheta = pi/number_of_divisions;

%% Potential V Computation

phi = linspace(0,2*pi,number_of_divisions);
theta = linspace(0,pi,number_of_divisions);
x = linspace(0,4,number_of_divisions/2);%cut sphere in half, better viewing
y = linspace(-4,4,number_of_divisions);
z = y;
[X,Y,Z] = meshgrid(x, y, z);

V_sphere = zeros(size(X)); % preallocate space for faster processing
for phi_increment = 1:number_of_divisions
    for theta_increment = 1:number_of_divisions
        Rx = X - R*sin(theta(theta_increment)).*cos(phi(phi_increment));
        Ry = Y - R*sin(theta(theta_increment)).*sin(phi(phi_increment));
        Rz = Z - R*cos(theta(theta_increment));
        
            cos_theta_sin_theta = sin(theta(theta_increment)); % Numerator: cos(theta(theta_increment)).*
            
            Rmag = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
            V_sphere = V_sphere + (cos_theta_sin_theta * dtheta * dphi/Rmag);
    end
end

V_sphere = k*sigma*(R^2)*V_sphere; % include scientific constants

%% Gradient Computation & Plot of E field of uniformly charged Spherical conductor

dx = x(2)-x(1);
dy = dx;
dz = dx;
[Ex, Ey, Ez] = gradient(-V_sphere, dx, dy, dz);

% Now creating indices to part select the specific coordinates 
% (and the components that correspond to those coordinates) that I want to 
% be graphed, otherwise the graph will put a vector at every X,Y,Z
index = 1:7:50;
index2 =1:7:25;
psphere = X(index, index2, index); 
p2 = Y(index, index2, index);
p3 = Z(index, index2, index);
p4 = Ex(index, index2, index);
p5 = Ey(index, index2, index);
p6 = Ez(index, index2, index);

%Emag = sqrt(Ex.^2 + Ey.^2 + Ez.^2); % not used
quiver3(psphere, p2, p3, p4, p5, p6, 'autoscalefactor', .8)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Spherical conductor with non-uniform charge, electric field')

%% Visualization of Potential Isosurfaces of uniformly charged Spherical conductor

levels_sphere = 4;
isovalues_sphere = [15 12 9 6]; % pick out specific, spread out isosurfaces to plot

for contour = 1:levels_sphere % for loop done 7 times
    [faces1, verts1, colors1] = isosurface(X, Y, Z, V_sphere/10^4, isovalues_sphere(contour), V_sphere/10^4); 
    % 7 isovalues used to create 7 isosurfaces, returns the number of faces 
    % and vertices of the isosurfaces, the last slot of the parentheses 
    % (V in this case) determines what variable the color is indicating the 
    % magnitude of, I divided by 10,000 to get "nicer" units
    
    p_sphere = patch('Vertices', verts1, 'Faces', faces1, 'FaceVertexCData',...
              colors1, 'FaceColor','interp','EdgeColor','interp'); 
          % creates polygons with verts specifying vertex values, 
          % and faces specifying which vertices to connect, interp means to 
          % "interpolate" or match the color of each face and edge from the 
          % adjacent face/edge to make the color continuous.
          % p_sphere is patch object that contains all the data for the polygons
          
    isonormals(X, Y, Z, V_sphere, p_sphere) % Calculates normals from the 
                                            % vertices of the patch identified 
                                            % by the handle p, this is used to 
                                            % create a smoother surface 
    p_sphere.EdgeColor = 'none';
    p_sphere.FaceAlpha = 0.3; % transparency of isosurface, higher is less transparent
end

daspect([1 1 1]) % sets aspect ratio for the axes
view(3); % this command is used to view the shapes created by the patch command
axis vis3d % freezes asoect ratio properties of axes
camlight % creates light which is to the right of the object by default unless specified otherwise
lighting gouraud % produces continuous shading between the faces of the 
                 % polygons that make up the constant potential surfaces
h1 = colorbar; % displays color spectrum to show relationship between color and magnitude
h1.Label.String = 'V \div 10^4'; % name of the color spectrum image
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title({'Planes of equipotential produced by', 'Spherical conductor with uniform charge distribution'})
% legend('Potential of spherical shell')

%% Plotting potential of a point charge for comparison
% This is just meant to show how the electric potential of a sphere and a
% point charge appear the same

q = 1.5*1.6*10^-19; % charge of a proton

V_point_charge = k*q./sqrt(X.^2 + Y.^2 + Z.^2);

levels_dipole = 5;
isovalues_dipole = [-10 -5 0 5 10]; 

for contour = 1:levels_dipole
    [faces2, verts2, colors2] = isosurface(X, Y, Z, V_dipole*10^32, isovalues_dipole(contour), V_dipole*10^32); 
    p_dipole = patch('Vertices', verts2, 'Faces', faces2, 'FaceVertexCData',... 
              colors2, 'FaceColor','interp','EdgeColor','interp');
    
    isonormals(X, Y, Z, V_dipole, p_dipole)  
    p_dipole.EdgeColor = 'none';
    p_dipole.FaceAlpha = 0.3; 
end
daspect([1 1 1]) 
view(3); 
axis vis3d 
camlight
lighting gouraud
h2 = colorbar;
h2.Label.String = 'V \cdot 10^{31}';
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')








