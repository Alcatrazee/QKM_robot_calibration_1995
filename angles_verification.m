angles = importdata('test_angles.txt');
RDK = Robolink();
robot = RDK.Item('QKM HL6');
angles = deg2rad(angles)
angles(:,3) = angles(:,3) + ones(10,1)*pi/2;
%good_angles = deg2rad(importdata('verified_angles.txt'));
for i=1:10
    robot.MoveJ(rad2deg(angles(i,:)));
end