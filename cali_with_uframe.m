%% brief of this program
% This program is for calibration simulation of robot arm QKM HL6 with
% RoboDK. This program can calibrate the robot automatically, with randomly
% generated joint angles, this program can reduce the composition to less
% than 1e-10.
% How to use it: open  RoboDK file(calibration_QKM_uframe.rdk),
% then run this script.

%% close all unecessary windows
close all
clear;
clc

%% RoboDK relative variables
global RDK;
global robot;
global tracker;
global robot_base;
RDK = Robolink();
robot = RDK.Item('QKM HL6');
robot_base = RDK.Item('QKM HL6 Base');
tracker = RDK.Item('Tracker Frame');

%% referenrce frame establisment
T_tracker_ref = Get_ref(RDK,tracker,robot_base);

%% parameters decleration
% norminal length of links
length_of_links = [491 350 350 84];                              

% norminal q 
q_vec_0 = [ -450 450 300;                         
            -450 450 length_of_links(1);
            -450 450 length_of_links(1)+length_of_links(2);
           length_of_links(3)-450 450 length_of_links(1)+length_of_links(2);
           length_of_links(3)-450 450 length_of_links(1)+length_of_links(2);
           length_of_links(3)-450 450 length_of_links(1)+length_of_links(2)]';

% norminal w
w_vec_0 = [0 0 1;
           0 1 0;
           0 1 0;
           1 0 0;
           0 1 0;
           1 0 0]';
       
%% g_st0 calculation(or M matrix from Park's paper)
HomePosition = [0 0 90 0 0 0]';
robot.MoveJ(HomePosition);
gst_0_in_tracker = Get_Calibration_plate_posture();
g_st0 = T_tracker_ref\gst_0_in_tracker;                                     % represents in reference frame

% calculate log(g_st0)
theta_M = norm(log_my(g_st0));
normed_M = log_my(g_st0)/theta_M;       
normed_M(1:3) = normed_M(1:3) - normed_M(1:3)'*normed_M(4:6)*normed_M(4:6);

%% norminal twist_matrix_0; size = 6 x 7
twist_matrix_0 = [[cross(q_vec_0,w_vec_0);w_vec_0],normed_M];               % nominal twist
twist_matrix_copy = twist_matrix_0;

%% read SMR positions and joint angles from files
num_of_pts = 10;
% theta_random_vec = GetRandomAngles(num_of_pts);
theta_random_vec = importdata('nice_angles.txt');
theta_random_vec_deg = rad2deg(theta_random_vec);
theta_random_vec(:,7) = ones(num_of_pts,1)*theta_M;                         % the 7th column shall be set to thetaM
samples = MeasurePosture(theta_random_vec_deg,num_of_pts);
theta_random_vec(:,3) = theta_random_vec(:,3) - ones(num_of_pts,1)*pi/2;

% transform SMR postures to reference frame.
for i=1:num_of_pts
    samples(:,:,i) = T_tracker_ref\samples(:,:,i);
end

%% variables declaration
df_f_inv = zeros(num_of_pts*6,1);                                           % df*f^-1
dp = zeros(num_of_pts*6,1);                                                 % deviation of configuration parameters
A = zeros(num_of_pts*6,42);                                                 % A matrix
j = 0;                                                                      % iteration time
norm_dp=[];                                                                 % for norm of dp visialization

%% parameters identification and composition
while j<20
    %% calculate A matrix and df*f^-1 (parameters identification)
    for i=1:num_of_pts                                                      % repeat num_of_pts times
        [T_n,~,~] = FK_new(twist_matrix_0,theta_random_vec(i,:));           % Tn calculation
        T_a = samples(:,:,i);                                               % in this case,we have ee postures represent in reference frame
        A(1+i*6-6:i*6,:) = A_matrix(twist_matrix_0,theta_random_vec(i,:));  % A matrix calculation
        df_f_inv(1+i*6-6:i*6) = log_my(T_a/T_n);                            % solve for log(df_f_inv)
    end
    dp = A\df_f_inv;                                                        % solve for dp(derive of twist)
    %% composition
    for i=1:6
        twist_matrix_0(:,i) = twist_matrix_0(:,i) + dp(1+i*6-6:i*6,1);                                                          % composition
        twist_matrix_0(:,i) = twist_matrix_0(:,i)/norm(twist_matrix_0(4:6,i));                                                  % normalization
        twist_matrix_0(1:3,i) = twist_matrix_0(1:3,i) - twist_matrix_0(1:3,i)'*twist_matrix_0(4:6,i)*twist_matrix_0(4:6,i);     % make v perpendicular to w
    end
    twist_matrix_0(:,7) = twist_matrix_0(:,7)+dp(37:42,1);                  % alternate the last rows
    %% data prepration of visialization
    j=j+1;                                                                  % counter plus 1 
    norm_dp = [norm_dp norm(dp)];                                           % minimization target value calculation
    disp ([num2str(j) 'th iteration'])                                                 % show number of iteration
    disp (norm(dp))                                                         % show value of norm of dp
    if norm(dp) < 1e-10                                                  	% quit the for loop if deviation is less than 1e-5
        break;
    end
    %% plot
    clf;                                                                    % clear plot
    draw_manipulator_my(twist_matrix_0,theta_M,'b');                        % draw nominal axis
    drawnow;
end
%% plot again
fig2 = figure(2);                                                           % create another window
bar3(norm_dp)                                                               % plot discrete data

%% verify the new twist
num_of_test_points = 10;
%test_angles = GetRandomAngles(num_of_test_points);                          % generate new points
test_angles = importdata('test_angles.txt');
test_angles_deg = rad2deg(test_angles);                                     % you know what it is
test_angles(:,7) = ones(num_of_test_points,1)*theta_M;                      % the 7th column shall be set to thetaM
truth_poses = MeasurePosture(test_angles_deg,num_of_test_points);           % launch sim on RoboDK for data
test_angles(:,3) = test_angles(:,3) - ones(num_of_test_points,1)*pi/2;      % for FK calculation
for i=1:num_of_pts
    truth_poses(:,:,i) = T_tracker_ref\truth_poses(:,:,i);
end

% variables decleration
estimated_poses = zeros(4,4,num_of_test_points);
deviation = zeros(4,4,num_of_test_points);                                  % quotient of measured value and estimated value
delta_eps = zeros(6,num_of_test_points);
norm_delta_eps = delta_eps;
delta_theta = zeros(1,num_of_test_points);
dis_deviation_sum = 0;

% start testing
for i=1:num_of_test_points
    %truth_poses(:,:,i) = T_tracker_ref\test_poses(:,:,i);                  % transform reference frame   
    [estimated_poses(:,:,i),~,~] = FK_new(twist_matrix_0,test_angles(i,:)); % estimate poses of EE with the same angles
    deviation(:,:,i) = truth_poses(:,:,i)/estimated_poses(:,:,i);           % calculate deviation of two poses
    dis_deviation_sum = dis_deviation_sum + norm(deviation(1:3,4,i));
    delta_eps(:,i) = log_my(deviation(:,:,i));   % norm of deviation of gesture
    delta_theta(i) = norm(delta_eps(:,i));
    norm_delta_eps(:,i) = delta_eps(:,i)/norm(delta_eps(:,i));
end
disp 'average delta theta (rad)'
norm(delta_theta)
disp 'average omega (mm)'
norm(norm_delta_eps(4:6,:))

q_origin = zeros(3,6);
q_after = zeros(3,6);
for i = 1:6
    q_after(:,i) = cross(twist_matrix_0(4:6,i),twist_matrix_0(1:3,i));
    q_origin(:,i) = cross(twist_matrix_copy(4:6,i),twist_matrix_copy(1:3,i));
end
dis_deviation_sum/num_of_test_points

%% some useful functions
% get reference frame using 3 SMRs
function posture = Get_ref(RDK,tracker,robot_base)
    Ref1 = RDK.Item('Ref1');
    Ref1.setParentStatic(tracker)
    posture1 = Ref1.Pose();
    Ref1.setParentStatic(robot_base);
    
    Ref2 = RDK.Item('Ref2');
    Ref2.setParentStatic(tracker)
    posture2 = Ref2.Pose();
    Ref2.setParentStatic(robot_base);
    
    Ref3 = RDK.Item('Ref3');
    Ref3.setParentStatic(tracker)
    posture3 = Ref3.Pose();
    Ref3.setParentStatic(robot_base);
    
    poses = [posture1(1:3,4) posture2(1:3,4) posture3(1:3,4)];
    posture = getReferenceFrame(poses,1);
end

% get end effector's posture
function posture = Get_Calibration_plate_posture()
    SMR1_position =  get_SMR_pos(1);
    SMR2_position =  get_SMR_pos(2);
    SMR3_position =  get_SMR_pos(3);
    SMRs_position = [SMR1_position SMR2_position SMR3_position];
    posture = getReferenceFrame(SMRs_position,3);
end

% calculate SMR poses with respect to reference frame
function pos = get_SMR_pos(SMR_n)
    global RDK;
    global robot;
    global tracker;
    if(SMR_n == 1)
        SMR = RDK.Item('SMR1');
        Pose = [     1.000000,    -0.000000,    -0.000000,   -70.000000 ;
                     -0.000000,     1.000000,     0.000000,   -70.000000 ;
                     -0.000000,     0.000000,     1.000000,    42.500000 ;
                      0.000000,     0.000000,     0.000000,     1.000000 ];
    elseif(SMR_n == 2)
        SMR = RDK.Item('SMR2');
        Pose = [     1.000000,    -0.000000,    -0.000000,   -70.000000 ;
                     -0.000000,     1.000000,     0.000000,    70.000000 ;
                     -0.000000,     0.000000,     1.000000,    42.500000 ;
                      0.000000,     0.000000,     0.000000,     1.000000 ];
    elseif(SMR_n == 3)
        SMR = RDK.Item('SMR3');
        Pose = [     1.000000,    -0.000000,    -0.000000,    70.000000 ;
                     -0.000000,     1.000000,     0.000000,   -70.000000 ;
                     -0.000000,     0.000000,     1.000000,    42.500000 ;
                      0.000000,     0.000000,     0.000000,     1.000000 ];
    end
    SMR.setParentStatic(tracker)
    posture = SMR.Pose();
    SMR.setParentStatic(robot);
    SMR.setPose(Pose);
    pos = posture(1:3,4);
end

% calculate random angles within joint limit
function theta_random_vec = GetRandomAngles(NumOfPts)
    theta_random_vec = randn(NumOfPts,7);
    Joint_limit = deg2rad([170 110 136 185 120 360]);
    for i = 1:NumOfPts
        for j=1:6
            while(abs(theta_random_vec(i,j))>Joint_limit(j))
                    theta_random_vec(i,j) = randn();
            end
        end
    end
end

% drive the robot in simulator to certain posture and return posture of
% plate's frame
function poses = MeasurePosture(angles,NumOfPts)
    global robot;
    poses = zeros(4,4,NumOfPts);
    for i = 1:NumOfPts
        robot.MoveJ(angles(i,:));
        poses(:,:,i) = Get_Calibration_plate_posture();
    end
end