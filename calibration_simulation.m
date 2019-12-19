%% close all unecessary windows
close all
clear;
clc
%% parameters definition
% length of links
length_of_links = [399 448 451 100];                              

% q(point on the rotation axis) vector without offset
q_vec_0 = [ 0 0 length_of_links(1);                         
            0 0 length_of_links(1);
            0 0 length_of_links(1)+length_of_links(2);
           length_of_links(3) 0 length_of_links(1)+length_of_links(2)+42;
           length_of_links(3) 0 length_of_links(1)+length_of_links(2)+42;
           length_of_links(3) 0 length_of_links(1)+length_of_links(2)+42]';
              

% omega matrix(axis of rotation)
w_vec_0 = [0 0 1;
           0 1 0;
           0 1 0;
           1 0 0;
           0 1 0;
           1 0 0]';
        
% g_st0
% g_st0 = [eye(3) [length_of_links(3)+82 0 length_of_links(1)+length_of_links(2)+42]';0 0 0 1];     
SMR_initial_file_name = 'initial_SMR_poses.txt';
angles_initial_file_name = 'initial_angles.txt';
[g_st0,~] = retrive_data(SMR_initial_file_name,angles_initial_file_name);
% mess up with axis 5 and 6  
q5 = [450.999998914246,-3.86878701577952e-06,889]';
w5 = [-0.0276250956103992,0.999618354219507,0]';
q6 =[532.968703662278;2.26525467221662;889];
w6 =[0.999618355813239;0.0276250379409021;0];
% calculate twist of M
theta_M = norm(log_my(g_st0));
normed_M = log_my(g_st0)/theta_M;       
normed_M(1:3) = normed_M(1:3) - normed_M(1:3)'*normed_M(4:6)*normed_M(4:6);

%% norminal twist_matrix_0; size = 6 x 7
twist_matrix_0 = [[cross(q_vec_0,w_vec_0);w_vec_0],normed_M];                                     % nominal twist
twist_matrix = twist_matrix_0;%+randn(6,7)/10;                                                    % actual twist
twist_matrix(:,5) = [-cross(w5,q5);w5];                                                           % 
twist_matrix(:,6) = [-cross(w6,q6);w6];
% [T_a,~,~] = FK_new(twist_matrix,[0 0 0 0 0 0 theta_M]);
twist_matrix(:,1:6) = twist_matrix(:,1:6)./vecnorm(twist_matrix(4:6,1:6)); % normalize omega(this line is for matlab before 2019a)

%% calculate normalized omega
vec_norm = vecnorm(twist_matrix(4:6,1:6));
for i = 1:6
    twist_matrix(:,i) = twist_matrix(:,i)./vec_norm(i);
end

%% calculate v again(by def,v is always perpendicular to omega)
for i=1:6
    twist_matrix(1:3,i) = twist_matrix(1:3,i) - twist_matrix(1:3,i)'*twist_matrix(4:6,i)*twist_matrix(4:6,i);       % make v perpendicular to w
end

%% define number of points,the more the better for calculation,but for sampling,the more the worse.

%% the last joint angle shall be set to 0


%% define the maximum and minimum of joint angle.
max_angle_vec = deg2rad(ones(1,6)*130);
min_angle_vec = deg2rad(-ones(1,6)*130);

%% define the random angle value for each random point
%theta_random_vec(:,1:6) = (max_angle_vec - min_angle_vec) .* rand(num_of_pts,6) + min_angle_vec;
% assign random angles
% front = max_angle_vec - min_angle_vec;
% for i=1:num_of_pts
%    for j=1:6
%        theta_random_vec(i,j) =  front(j)*rand()+min_angle_vec(j);
%    end
% end

%% read data from files
data_file_name = 'SMR_poses.txt';
angle_file_name = 'old_angles.txt';
[samples,theta_random_vec] = retrive_data(data_file_name,angle_file_name);
theta_random_vec = deg2rad(theta_random_vec);
num_of_pts = length(theta_random_vec);  
theta_random_vec(:,7) = ones(num_of_pts,1)*theta_M;         % the 7th column shall be set to thteaM
%% variables declaration
df_f_inv = zeros(num_of_pts*6,1);                   % df*f^-1
dp = zeros(num_of_pts*6,1);                         % deviation of configuration parameters
A = zeros(num_of_pts*6,42);                         % A matrix
j = 0;                                              % iteration time
norm_dp=[]; % for norm of dp visialization
%% parameters identification and composition
while j<1000
    %% calculate A matrix and df*f^-1 (parameters identification)
    for i=1:num_of_pts                                                      % repeat num_of_pts times
        [T_n,~,~] = FK_new(twist_matrix_0,theta_random_vec(i,:));           % Tn calculation
        %[T_a,~,~] = FK_new(twist_matrix,theta_random_vec(i,:));             % Ta calculation
        T_a = samples(:,:,i);
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
    twist_matrix_0(:,7) = twist_matrix_0(:,7)+dp(37:42,1);                          % alternate the last rows
    %% data prepration of visialization
    j=j+1;                                                                  % counter plus 1 
    norm_dp = [norm_dp norm(dp)];                                           % minimization target value calculation
    disp (norm(dp))                                                         % show value of norm of dp
    disp (j)                                                                % show number of iteration
    if norm(dp) < 10e-6                                         % quit the for loop if deviation is less than 1e-5
        break;
    end
    %% plot
    clf;                                                                    % clear plot
%     draw_manipulator(twist_matrix,'r');                                     % draw actual axis 
    draw_manipulator_my(twist_matrix_0,theta_M,'b');                                   % draw nominal axis
    drawnow;
end
%% plot again
fig2 = figure(2);                                                           % create another window
bar3(norm_dp)                                                               % plot discrete data
