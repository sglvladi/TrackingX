
target_position = [sqrt(8);sqrt(8);0];

sensor_position = [2, -1, 0]';
sensor_orientation = [pi/4, 0, 0]';

Rx = @(t) [1, 0, 0; 0, cos(t), -sin(t); 0, sin(t), cos(t)];
Ry = @(t) [cos(t), 0, sin(t); 0, 1, 0; -sin(t), 0, cos(t)];
Rz = @(t) [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];

rot_mat = Rz(-sensor_orientation(1))*Ry(-sensor_orientation(2))*Rx(-sensor_orientation(3));

xyz = target_position - sensor_position;

xyz_rot = rot_mat*xyz;

x = xyz_rot(1);
y = xyz_rot(2);
z = xyz_rot(3);

[theta,rho,z]  = cart2pol(x,y,z);
angle = rad2deg(theta);
% theta = theta - sensor_orientation(1);

figure;
plot([0,target_position(1)],[0,target_position(2)]);
hold on
plot([sensor_position(1),target_position(1)],[sensor_position(2),target_position(2)]);



