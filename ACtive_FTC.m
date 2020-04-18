%Controller Design for Ground Target tracking
%Written by Shen Qiang
clear all; close all; clc; format long
global DistTorq  MOI  MOI_Unc  U_RW  MA  rate_desire  S  PHI  mu  U_RW_Command
global Est_L Est_G
%------------- satellite moment of inertia modelling ------------------
% MOI = [ 130,  6.5,    6;...
%         6.5,  140,    5.5;...
%         6,    5.5,    135];

% MOI = [ 350,        3,      4;...
%     3,          270,     10;...
%     4,       10,      190];

MOI = [  10,         1.2,      0.5;...
    1.2,        19,       1.5;...
    0.5,       1.5,      25 ];

MOI_Nom = MOI ;

MOI_Err = 0*diag([7, 8, -5]) ;  %no MOI uncertainty

MOI_Unc = MOI_Nom + MOI_Err ;

MA = 1/sqrt(3)*[-1,     -1,     1,      1;...  %Mounting matrix of four RWs
    1,     -1,   -1,      1;...
    1,       1,    1,      1];

% specify the simulation parameters
dT = 0.1; %sampling time
Tsimu = 200; %Simulation time, seconds
m = round(Tsimu/dT+1);
XTime = linspace(0,Tsimu,m)';

% High-level controller selection parameter
ContrIndex = 0;
% if =0, use PD controller
% if =1, use Adaptive Robust controller

% Disturbance type
DistIndex = 1;
% if = 0, Slow-varying disturbance
% if = 1; Fast-varying disturbance

% Fault estimation switch
Fault_EST_Index = 0;
% if = 0, do not activate fault estimation
% if = 1, activate fault estimation

% ---------------------- High lever controller ----------------------------
%Case 0: PD controller
%specify the controller parameters
damp = 1;  %damping ratio
ts = 30;   %settling time
wn = 8/ts;  % natural frequency
%wn=0.45
k = 2*wn^2;
d = 2*damp*wn;
K = 1*k*MOI; %control gain
D = 1*d*MOI; %control gain

%Case 1: Robust adaptive controller
%specify the controller parameters
RAC_alpha = 0.2;
RAC_beta = 1.8;
RAC_epsilon0 = norm( pinv(MA) );
RAC_k = 100;
RAC_epsilon1 = 0.1;

RAC_h0 = 0.01;
mu = 0.01;

f_hat = zeros(3,1);

% ---------------------- End of  High lever controller --------------------

% %Maneuver angle in orbital frame
% Roll_tilt=0*pi/180;

%sensor noise and mounting misalignment
%Star Tracker measurment noise
STS_noise = 1*0.005*pi/180;
% STS_noise = 0*0.005*pi/180; %no star tracking sensor noise

%Star Tracker Mounting misalignment
STSMis_Yaw = 1*0.005*pi/180;
STSMis_Pitch = 1*0.8*0.005*pi/180;
STSMis_Roll = -1*0.005*pi/180;
% STSMis_Yaw = 0*0.005*pi/180; %no Star Tracker Mounting misalignment
% STSMis_Pitch = 0*0.8*0.005*pi/180; %no Star Tracker Mounting misalignment
% STSMis_Roll = 0*0.005*pi/180; %no Star Tracker Mounting misalignment

%FOG measurement noise
Rate_noise = 1*0.0005*pi/180;

% Rate_noise = 0*0.0005*pi/180; %no FOG measurement noise

%Mounting matrix of FOGs
FOGMis = [  cos(0.1*pi/180), cos(89.9*pi/180),cos(90*pi/180);...
    cos(89.9*pi/180),cos(0.1*pi/180),cos(89.9*pi/180);...
    cos(90*pi/180),cos(90*pi/180),cos(0.1*pi/180)];
% FOGMis = eye(3); %no FOG Mounting misalignment

%FOG bias
bias = 1*1/3600*pi/180*[1;1;1]; %1 deg/hour
% bias = 0*1/3600*pi/180*[1;1;1]; %no FOG bias

%Reaction wheel specification
%Maximum control torque
u_max = 0.2;

%Fault Detection Threshold
threshold = 2.0e-3;

tic; %Timer

j = 0;
%set initiall satellite attitude
q1_initial = -0.5;
q2_initial = 0.3;
q3_initial = -0.4;
q0_initial = sqrt(1-(q1_initial^2+q2_initial^2+q3_initial^2));

%initial measured attitude and rate
Q = [q1_initial, q2_initial, q3_initial, q0_initial]';
rate = [0.005 0.006 0.004]';
%initial fault detection
Omega_Detec_hat_0 = rate;
Omega_Est_hat_0 = rate;
Aux_hat_0 = zeros(3,1);
%Initial desired quaternion
Qc = [0 0 0 1]';
%initial RW momentum
H_RW = [1;1;-1;-1];

% LMI Computation
lmi_solve_L;
%%%%%%%%%%%


for i=0:dT:Tsimu %time instant
    j=j+1; %iteration indicator
    
    %----------------------------disturbance modelling------------------
    %Total disturbance torque
    if DistIndex == 0
        %  DistTorq = [0.01; 0.02; -0.01];
        DistTorq = zeros(3,1);
    elseif DistIndex == 1
        %        omega_0 = 0.01;
        %         DistTorq = 1e-3*[ 3*cos(10*omega_0*j*dT)+4*sin(3*omega_0*j*dT)-1;
        %             -1.5*sin(2*omega_0*j*dT)-3*cos(5*omega_0*j*dT)+1.5;
        %             2*sin(10*omega_0*j*dT)-1.5*cos(4*omega_0*j*dT)+1;];
        
        DistTorq = [ -0.005 *sin(j*dT); 0.005 *sin(j*dT); -0.005 *sin(j*dT) ];
        %  DistTorq = zeros(3,1);
    end
    %-------------------End of disturbance modelling------------------
    
    % ----------- Fault simulation information ---------------
    % four RWs happen fault at time: t1_LOE_Start, t1_Bias_Start ;
    %                                t2_LOE_Start, t2_Bias_Start ;
    %                                t3_LOE_Start, t3_Bias_Start ;
    %                                t4_LOE_Start, t4_Bias_Start ; respectively
    %
    % multiplicative fault --> loss of actuator effectiveness
    % 1-e1, 1-e2, 1-e3, 1-e4
    % additive fault --> constant bias
    % f1, f2, f3, f4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t1_LOE_Start = 5;
    t1_Bias_Start = 0;
    
    e1 = 0.6 ;
    % e1 = 0 ;
    f1 = 0 ;
    %%%%%%%%%%%%
    t2_LOE_Start = 5 ;
    t2_Bias_Start = 100 ;
    
    e2 = 0.2 ;
    f2 = -0.03 ;
    %     e2 = 0 ;
    %     f2 = 0 ;
    %%%%%%%%%%%%
    t3_LOE_Start = 0 ;
    t3_Bias_Start = 100 ;
    
    e3 = 0 ;
    f3 =  -0.06 ;
    %     e3 = 0 ;
    %     f3 = 0 ;
    
    %%%%%%%%%%%%
    t4_LOE_Start = 0;
    t4_Bias_Start = 0 ;
    
    e4 = 0;
    f4 = 0;
    %-------------- END of fault information --------------------
    
    Q_meas = Q;
    rate_meas = rate;
    
    %-------------- The desired rate ----------------
    rate_desire = zeros(3,1);
    rate_d(:,j) = rate_desire;
    %-------------End of The desired rate -------------
    
    %--------- calculate the derivative of deisre rate -----------
    % numerical calculation / (can also be analytical)
    if (j==1)||(j==2)
        rate_dDot(:,j) = [0 0 0]';
    end
    if j>2
        rate_dDot(:,j) = (rate_d(:,j)-rate_d(:,j-1))/dT;
    end
    %--------- End of calculate the derivative of deisre rate -----------
    
    %---------------Attitude command generation-------------------------
    Y = ode4(@Quat_Desire_Propagator,[(j-1)*dT j*dT],Qc);
    Qc = Y(end,1:4)';
    qc1 = Qc(1);
    qc2 = Qc(2);
    qc3 = Qc(3);
    qc0 = Qc(4);
    Quat_d(:,j) = Qc;
    %---------------End of Attitude command generation-------------------------
    
    %------------------ calculate error quaternion -------------
    Qc_inv = [-qc1, -qc2, -qc3, qc0]';
    Qe = QuatMultiplication(Qc_inv, Q_meas);
    Qe_Vec = [Qe(1);Qe(2);Qe(3)];% error for controller
    Qe_Scal = Qe(4);
    Qe_Vec_x = SkewSymmetric(Qe_Vec);
    Quat_error(:,j) = Qe;
    %------------------ End of calculate error quaternion -------------
    
    %------- calculate rate error -------------
    RotMat = (Qe_Scal^2-Qe_Vec'*Qe_Vec)*eye(3)+2*Qe_Vec*Qe_Vec'-2*Qe_Scal*Qe_Vec_x;
    Rate_error(:,j) = rate_meas-RotMat*rate_d(:,j);    %error for controller
    %------------- End of calculate rate error -------------
    
    EulerAngleErr(j,:)=SpinCalc('QtoEA321',Qe',1e-6,1);
    for i=1:3
        if EulerAngleErr(j,i)>100
            EulerAngleErr(j,i)=EulerAngleErr(j,i)-360;
        end
    end
    
    %*********************************************************************
    %-----------------------Controller Design---------------------------------
    %************************************************************************
    
    if ContrIndex == 0;
        %---------------- Implementation of PD-like controller ------------------
        %Control input for attitude control
        rate_meas_x = SkewSymmetric(rate_meas);
        U_PD = rate_meas_x*MOI_Nom*RotMat*rate_desire+MOI_Nom*RotMat*rate_dDot(:,j) ...
            -K*Qe_Vec-D*Rate_error(:,j) ;  % H_RW = Jrw*omega_RWs
        
        U_RW_Command = pinv(MA)*U_PD;
        
        % RW Saturation
        if norm(U_RW_Command, inf) >= u_max
            U_RW_Command = u_max*U_RW_Command/norm(U_RW_Command, inf);
        end
        %----------------End of Implementation of PD-like controller --------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif ContrIndex == 1;
        %---------------- Implementation of Reconfigurable controller ------------------
        %Control input for attitude control
        S = Rate_error(:,j) + RAC_alpha*atan(RAC_beta*Qe_Vec);
        % S =  Rate_error(:,j) + RAC_alpha*Qe_Vec;
        PHI = 1+norm(rate_meas)+norm(rate_meas)^2;
        %---------------parameter adaptive law-------------------------
        Y = ode4(@AdapLaw,[(j-1)*dT j*dT], RAC_h0);
        h_hat = Y(end);
        RAC_h0 = h_hat;
        h_hat_out(:,j) = h_hat;
        
        RAC_epsilon2 = mu / PHI;
        
        YITA = RAC_k + S'*f_hat/(norm(S)^2+RAC_epsilon1^2) + h_hat*PHI/(norm(S)+RAC_epsilon2);
        
        if ( norm(S) >= u_max/(RAC_epsilon0*YITA) )
            Sat_YITA = S/norm(S);
        else
            Sat_YITA = RAC_epsilon0*YITA*S/u_max;
        end
        
        U_RAC = - u_max/RAC_epsilon0 * pinv(MA)*Sat_YITA  ;
        
        U_RW_Command = U_RAC;
        %----------------End of Implementation of Robust adaptive controller --------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %------------ Begin of RWs Fault Modelling --------------
    % four RWs happen fault at time t1,t2,t3,t4 respectively
    % multiplicative fault --> loss of actuator effectiveness
    %
    % additive fault --> constant bias
    %
    % The following is corresponding to no falut case
    %     t1_LOE_Start = Tsimu;
    %     t1_Bias_Start = Tsimu;
    %     e1_RW_Eff = 1 ;
    %     u1_RW_Bias = 0 ;
    %
    %     t2_LOE_Start = Tsimu;
    %     t2_Bias_Start = Tsimu;
    %     e2_RW_Eff = 1 ;
    %     u2_RW_Bias = 0 ;
    %
    %     t3_LOE_Start = Tsimu;
    %     t3_Bias_Start = Tsimu;
    %     e3_RW_Eff = 1 ;
    %     u3_RW_Bias = 0 ;
    %
    %     t4_LOE_Start = Tsimu;
    %     t4_Bias_Start = Tsimu;
    %     e4_RW_Eff = 1 ;
    %     u4_RW_Bias = 0 ;
    
    %---------- 1st RW ------------
    if j >= round(t1_LOE_Start/dT+1)
        e1_RW_Eff = e1 ;
    else
        e1_RW_Eff = 0 ;
    end
    
    if j >= round(t1_Bias_Start/dT+1)
        u1_RW_Bias = f1 ;
    else
        u1_RW_Bias = 0;
    end
    
    %---------- 2ed RW ------------
    if j >= round(t2_LOE_Start/dT+1)
        e2_RW_Eff = e2 ;
    else
        e2_RW_Eff = 0 ;
    end
    
    if j >= round(t2_Bias_Start/dT+1)
        u2_RW_Bias = f2 ;
    else
        u2_RW_Bias = 0 ;
    end
    
    %---------- 3rd RW ------------
    if j >= round(t3_LOE_Start/dT+1)
        e3_RW_Eff = e3 ;
    else
        e3_RW_Eff = 0 ;
    end
    
    if j >= round(t3_Bias_Start/dT+1)
        u3_RW_Bias = f3 ;
    else
        u3_RW_Bias = 0 ;
    end
    
    %---------- 4th RW ------------
    if j >= round(t4_LOE_Start/dT+1)
        e4_RW_Eff = e4 ;
    else
        e4_RW_Eff = 0 ;
    end
    
    if j >= round(t4_Bias_Start/dT+1)
        u4_RW_Bias = f4 ;
    else
        u4_RW_Bias = 0 ;
    end
    
    % Effectiveness matrix of four RWs
    e_RW_Eff = diag([e1_RW_Eff, e2_RW_Eff, e3_RW_Eff, e4_RW_Eff]);
    
    % Bias vector of four RWs
    u_RW_Bias = [u1_RW_Bias; u2_RW_Bias; u3_RW_Bias; u4_RW_Bias] ;
    
    % Faulty RWs output control torque
    
    U_RW = ( eye(4) - e_RW_Eff ) * U_RW_Command + u_RW_Bias ;
    
    F_Total = MA * ( -e_RW_Eff  * U_RW_Command + u_RW_Bias ) ;
    f_total_out(:,j) = F_Total;
    %---------- End of RWs Fault Modelling -------------
    
    e_RW_Eff_out(:,j) = [ e_RW_Eff(1,1); e_RW_Eff(2,2); e_RW_Eff(3,3); e_RW_Eff(4,4) ] ;
    u_RW_Bias_out(:,j) = u_RW_Bias ;
    U_RW_Command_out(:,j) = U_RW_Command ;
    RWU(:,j) = U_RW; % Commanded reaction wheel control toque
    H_RW_out(:,j) = H_RW;
    
    %------------Attitude Propagator-------------------------------
    Y = ode4(@attitudedynamics,[(j-1)*dT j*dT],[rate; Q; H_RW]');
    rate = Y(end,1:3)'; %rate is the nominal angular velocity without sensor noise
    Q = Y(end,4:7)'; %Q is nominal quaternion without sensor noise
    H_RW = Y(end,8:11)'; %RW momentum
    %------------End of Attitude Propagator-------------------------------
    
    
    %-----------Sensor Modelling------------------------------
    %Rate sensor noise and bias
    w_noise(:,j) = Rate_noise*[randn;randn;randn];
    rate_meas = FOGMis*(rate+w_noise(:,j)+bias);
    % rate_meas = rate;
    
    %Star tracker measurement noise and mountig misalignment
    xnoise = randn;
    DCM_Yaw_noise = [cos(STS_noise*xnoise+STSMis_Yaw), sin(STS_noise*xnoise+STSMis_Yaw),  0;...
        -sin(STS_noise*xnoise+STSMis_Yaw),cos(STS_noise*xnoise+STSMis_Yaw),       0;...
        0 ,               0,               1];
    ynoise = randn;
    DCM_Pitch_noise=[cos(STS_noise*ynoise+STSMis_Pitch), 0,             -sin(STS_noise*ynoise+STSMis_Pitch);...
        0,                      1,             0;...
        sin(STS_noise*ynoise+STSMis_Pitch) ,      0,               cos(STS_noise*ynoise+STSMis_Pitch)];
    znoise = randn;
    DCM_Roll_noise = [1,                   0,                       0;...
        0,                     cos(STS_noise*znoise+STSMis_Roll),       sin(STS_noise*znoise+STSMis_Roll);...
        0 ,                   -sin(STS_noise*znoise+STSMis_Roll),         cos(STS_noise*znoise+STSMis_Roll)];
    
    DCM_noise = DCM_Roll_noise*DCM_Pitch_noise*DCM_Yaw_noise;
    
    Qnoise = SpinCalc('DCMtoQ',DCM_noise,1e-3,0);
    EulerAngleNoise(:,j) = SpinCalc('QtoEA321',Qnoise,1e-6,1);
    
    for k = 1:3
        if EulerAngleNoise(k,j)>100
            EulerAngleNoise(k,j) = EulerAngleNoise(k,j)-360;
        end
    end
    
    QnoiseKron = [Qnoise(4), Qnoise(3), -Qnoise(2), Qnoise(1);...
        -Qnoise(3), Qnoise(4), Qnoise(1), Qnoise(2);...
        Qnoise(2), -Qnoise(1), Qnoise(4), Qnoise(3);...
        -Qnoise(1), -Qnoise(2), -Qnoise(3), Qnoise(4)];
    
    Q_meas = QnoiseKron*Q;
    Q_meas= Q_meas/norm(Q_meas); %norm(Q_means should be almost euqal to 1)
    %-----------End of sensor modelling------------------------------
    
    %--------- Fault Detection ---------------------------
    Y = ode4(@Detection,[(j-1)*dT j*dT], Omega_Detec_hat_0', rate_meas);
    
    Omega_Detec_hat = Y(end,1:3)';
    Omega_Detec_hat_0 = Omega_Detec_hat;
    Omega_Detec_hat_out(:,j) =  Omega_Detec_hat;
    
    Omega_Detec_error = norm(Omega_Detec_hat - rate_meas);
    Omega_Detec_error_out(:,j) = Omega_Detec_error;
    
    %     if ( Omega_Detec_error  > threshold   &&  Fault_EST_Index == 0)
    %         Fault_EST_Index = 1;
    %     elseif ( Omega_Detec_error <  threshold   &&  Fault_EST_Index == 1)
    %         Fault_EST_Index = 0;
    %     end
    
    if ( Omega_Detec_error  >  threshold && Fault_EST_Index == 0&& ContrIndex == 0 )
        Fault_EST_Index = 1;
        iter_start_est = j;
        Omega_Est_hat_0 = rate_meas;
    end
    %---------- End of Fault Detection -------------------------
    
    %---------- Fault Estimation ---------------------------------
    if (Fault_EST_Index == 1)
        
        % state and  fault estimation
        Est_additional_para = [f_hat; rate_meas]';
        Y = ode4(@Estimation,[(j-1)*dT j*dT], [ Omega_Est_hat_0;  Aux_hat_0]', Est_additional_para);
        
        Omega_Est_hat = Y(end,1:3)';
        Omega_Est_hat_0 = Omega_Est_hat;
        Omega_Est_out(:,j) = Omega_Est_hat;
        Omega_Est_error_out(:,j) = Omega_Est_hat - rate_meas;
        
        Aux_hat = Y(end,4:6)';
        Aux_hat_0 = Aux_hat;
        Aux_hat_out(:,j) = Aux_hat;
        
        f_hat = Aux_hat + Est_G*MOI_Unc*Omega_Est_hat;
        f_hat_out(:,j) = f_hat;
        f_error_out(:,j) = F_Total - f_hat;
        
        if ( j - iter_start_est >= 1)
            f_hat_prev = f_hat_out(:,j-1);
            if ( norm(Omega_Est_hat - rate_meas) + norm(f_hat - f_hat_prev) < 0.002  && ContrIndex == 0 ) %
                ContrIndex = 1;
                iter_start_reconf = j;
            end
        end
        
    else
        f_hat = zeros(3,1);
        f_hat_out(:,j) = f_hat;
    end
    
    % controller reconfiguration
    
    % ---------------End of Fault Estimation ------------------
    
    % disp(ContrIndex);
    
end

toc

% figure(101)
% plot (XTime,rate_d(1,:),'r','LineWidth',1.5);hold on
% plot (XTime,rate_d(2,:),'b','LineWidth',1.5);hold on
% plot (XTime,rate_d(3,:),'g','LineWidth',1.5);
% legend('\omega_{d1}','\omega_{d2}','\omega_{d3}')
% title('Desired angular velocity');
% figure(102)
% plot (XTime,rate_dDot(1,:),'r','LineWidth',1.5);hold on
% plot (XTime,rate_dDot(2,:),'b','LineWidth',1.5);hold on
% plot (XTime,rate_dDot(3,:),'g','LineWidth',1.5);
% l1 = legend('$$\dot{\omega}_{d1}$$','$$\dot{\omega}_{d2}$$','$$\dot{\omega}_{d3}$$');
% set(l1,'Interpreter','latex','FontSize',10);
% title('Desired angular velocity derivative');

figure(9)
plot(XTime,EulerAngleErr(:,3),'r','LineWidth',1.5);hold on
plot(XTime,EulerAngleErr(:,2),'b--','LineWidth',1.5);hold on
plot(XTime,EulerAngleErr(:,1),':','Color',[0,0.4,0],'LineWidth',1.5);hold on
xlabel('Time (s)','FontSize',12,'FontWeight','bold');%,'FontWeight','bold'
l9 = legend('Roll','Pitch','Yaw');
set(l9,'Interpreter','latex','FontSize',16,'Location', 'best'); % , 'Position',[0.785,0.755,0.045,0.15]
ylabel('Attitude (deg)','FontSize',12,'FontWeight','bold'); %
hold off

figure(10)
plot (XTime,Rate_error(1,:),'r','LineWidth',1.5);hold on
plot (XTime,Rate_error(2,:),'b--','LineWidth',1.5);hold on
plot (XTime,Rate_error(3,:),':','Color',[0,0.4,0],'LineWidth',1.5);
l10 = legend('$$\omega_{1}$$','$$\omega_{2}$$','$$\omega_{3}$$');
set(l10,'Interpreter','latex','FontSize',16,'Location','Best');
%title('Angular velocity error');
xlabel('Time (s)','FontSize',12,'FontWeight','bold');
ylabel('Angular velocity (rad/s)','FontSize',12,'FontWeight','bold');

% figure(111)
% plot (XTime,Quat_d(1,:),'r','LineWidth',1.5);hold on
% plot (XTime,Quat_d(2,:),'b','LineWidth',1.5);hold on
% plot (XTime,Quat_d(3,:),'g','LineWidth',1.5);hold on
% plot (XTime,Quat_d(4,:),'y','LineWidth',1.5);
% legend('q_{d1}','q_{d2}','q_{d3}','q_{d0}')
% title('Desired quaternion');

% figure(11)
% plot (XTime,Quat_error(4,:),'k','LineWidth',1.5);hold on
% plot (XTime,Quat_error(1,:),'r','LineWidth',1.5);hold on
% plot (XTime,Quat_error(2,:),'b','LineWidth',1.5);hold on
% plot (XTime,Quat_error(3,:),'g','LineWidth',1.5);
% legend('q_{e0}','q_{e1}','q_{e2}','q_{e3}')
% %title('Quaternion error');
% xlabel('Time (s)','FontSize',12,'FontWeight','bold');
% ylabel('Attitude tracking error','FontSize',12,'FontWeight','bold');


% figure(19);
% subplot(2,2,1);plot(XTime,RWU(1,:),'LineWidth',1.5);%title('Torque of RW1, Nm')
% ylabel('Reaction output Torque u_{1}, Nm');
% subplot(2,2,2);plot(XTime,RWU(2,:),'LineWidth',1.5);%title('Torque of RW2, Nm')
% ylabel('Reaction output Torque u_{2}, Nm');
% subplot(2,2,3);plot(XTime,RWU(3,:),'LineWidth',1.5);%xlabel('Time, sec');title('Torque of RW3, Nm')
% xlabel('Time (s)','FontSize',12,'FontWeight','bold');
% ylabel('Reaction output Torque u_{3}, Nm');
% subplot(2,2,4);plot(XTime,RWU(4,:),'LineWidth',1.5);%xlabel('Time, sec');title('Torque of RW4, Nm')
% xlabel('Time (s)','FontSize',12,'FontWeight','bold');
% ylabel('Reaction output Torque u_{4}, Nm');

% figure(24);
% plot(XTime,e_RW_Eff_out,'LineWidth',1.5);
% title('Effectiveness')
% xlabel('Time, sec');
%
% figure(25);
% plot(XTime,u_RW_Bias_out,'LineWidth',1.5);
% title('Bias')
% xlabel('Time, sec');


figure(21)
plot (XTime,U_RW_Command_out(1,:),'r','LineWidth',1.5);hold on
plot (XTime,U_RW_Command_out(2,:),'b--','LineWidth',1.5);hold on
plot (XTime,U_RW_Command_out(3,:),':','Color',[0,0.4,0],'LineWidth',1.5);hold on
plot (XTime,U_RW_Command_out(4,:),'k-.','LineWidth',1.5);
l21 = legend('$$u_{c1}$$','$$u_{c2}$$','$$u_{c3}$$','$$u_{c4}$$');
set(l21,'Interpreter','latex','FontSize',16, 'Location', 'Best');  %    ,'Position',[0.79,0.13,0.045,0.2]
ylabel('Command control torque (Nm)','FontSize',12,'FontWeight','bold')
xlabel('Time (s)','FontSize',12,'FontWeight','bold');


% figure(26);
% subplot(2,2,1);plot(XTime,U_RW_Command_out(1,:),'LineWidth',1.5);%title('Command Torque of RW1, Nm')
% ylabel('Command Torque u_{c1}, Nm');
% subplot(2,2,2);plot(XTime,U_RW_Command_out(2,:),'LineWidth',1.5);%title('Command Torque of RW2, Nm')
% ylabel('Command Torque u_{c2}, Nm');
% subplot(2,2,3);plot(XTime,U_RW_Command_out(3,:),'LineWidth',1.5);%xlabel('Time, sec');title('Command Torque of RW3, Nm')
% xlabel('Time (s)','FontSize',12,'FontWeight','bold');
% ylabel('Command Torque u_{c3}, Nm');
% subplot(2,2,4);plot(XTime,U_RW_Command_out(4,:),'LineWidth',1.5);%xlabel('Time, sec');title('Command Torque of RW4, Nm')
% xlabel('Time (s)','FontSize',12,'FontWeight','bold');
% ylabel('Command Torque u_{c4}, Nm');
% set(gcf,'Position',[100,300,700,500]);


% figure(27);
% plot(XTime,H_RW_out,'LineWidth',1.5);
% title('RW momentum')
% xlabel('Time, sec');


figure(50);
plot(XTime, f_total_out(1,:)','r', 'LineWidth',1.5);
hold on;
plot(XTime, f_total_out(2,:)', 'b--', 'LineWidth',1.5);
hold on;
plot(XTime, f_total_out(3,:)', ':','Color',[0,0.4,0], 'LineWidth',1.5);
hold off;
l50 = legend('$${f_1}$$','$$f_2$$','$$f_3$$' );
set(l50,'Interpreter','latex','FontSize',16,'Location', 'best');  %
ylabel('Total fault effects','FontSize',12,'FontWeight','bold');
xlabel('Time (s)','FontSize',12,'FontWeight','bold');

figure(51)
plot(XTime, Omega_Detec_error_out, 'b', 'LineWidth',1.5);
hold on;
plot(XTime, threshold*ones(length(XTime)), '-- r', 'LineWidth',1.5 );
l51 = legend('$${\|  {  \bf \tilde{\mathbf{\omega} } }_{b,d}\|}$$ in fault detection', 'Threshold');
set(l51,'Interpreter','latex','FontSize',16,'Location', 'best','FontWeight','bold' );  %
hold off;
ylabel('Fault detection residual','FontSize',12,'FontWeight','bold')
xlabel('Time (s)','FontSize',12,'FontWeight','bold');


figure(52);
plot(XTime, f_hat_out(1,:)', 'r', 'LineWidth',1.5);
hold on;
plot(XTime, f_hat_out(2,:)', 'b--', 'LineWidth',1.5);
hold on;
plot(XTime, f_hat_out(3,:)', ':','Color',[0,0.4,0], 'LineWidth',1.5);
hold off;
l52 = legend('$$\hat{f}_1$$','$$\hat{f}_2$$','$$\hat{f}_3$$' );
set(l52,'Interpreter','latex','FontSize',16, 'Location', 'best' );  %
ylabel('Estimated total fault','FontSize',12, 'FontWeight','bold')
xlabel('Time (s)','FontSize',12,'FontWeight','bold');

figure(53);
plot(XTime, f_error_out(1,:)', 'r', 'LineWidth',1.5);
hold on;
plot(XTime, f_error_out(2,:)', 'b--', 'LineWidth',1.5);
hold on;
plot(XTime, f_error_out(3,:)', ':','Color',[0,0.4,0], 'LineWidth',1.5);
hold off;
l53 = legend('$$\tilde{f}_1$$','$$\tilde{f}_2$$','$$\tilde{f}_3$$' );
set(l53,'Interpreter','latex','FontSize',16, 'Location', 'best' );  %
ylabel('Fault estimation error','FontSize',12,'FontWeight','bold')
xlabel('Time (s)','FontSize',12,'FontWeight','bold');

figure(54);
plot(XTime, Omega_Est_error_out(1,:)', 'r', 'LineWidth',1.5);
hold on;
plot(XTime, Omega_Est_error_out(2,:)', 'b--', 'LineWidth',1.5);
hold on;
plot(XTime, Omega_Est_error_out(3,:)', ':','Color',[0,0.4,0], 'LineWidth',1.5);
hold off;
l54 = legend('$${ \tilde{\omega}_{{b,i1}} }$$','$$\tilde{\omega}_{b,i2}$$','$$\tilde{\omega}_{b,i3}$$' );
set(l54,'Interpreter','latex','FontSize',16, 'Location', 'Best');  %
ylabel('Angular velocity estimation error in fault estimation','FontSize',16,'FontWeight','bold')
xlabel('Time (s)','FontSize',16,'FontWeight','bold');




% figure(41);
% plot(XTime,d_err_out,'LineWidth',1.5);
% l41 = legend('$$\tilde{d}_{1}$$','$$\tilde{d}_{2}$$','$$\tilde{d}_{3}$$');
% set(l41,'Interpreter','latex','FontSize',14,'Position',[0.795,0.73,0.045,0.16]);
% ylabel({'$\tilde{d}$'},'Interpreter','latex','FontSize',18,'FontWeight','bold');
% xlabel('Time (s)','FontSize',12,'FontWeight','bold');
%
% figure(42);
% plot(XTime,s_obs_error_out,'LineWidth',1.5);
% l42 = legend('$$\tilde{s}_{1}$$','$$\tilde{s}_{2}$$','$$\tilde{s}_{3}$$');
% set(l42,'Interpreter','latex','FontSize',14,'Position',[0.795,0.75,0.045,0.16]);
% ylabel({'$\tilde{s}$'},'Interpreter','latex','FontSize',18,'FontWeight','bold');
% xlabel('Time (s)','FontSize',12,'FontWeight','bold');

% figure(27);
% plot(XTime,b_hat_out,'LineWidth',1.5);
% title('b_{hat}')
% xlabel('Time, sec');

% figure(28);
% plot(XTime,beta_out,'LineWidth',1.5);
% title('beta')
% xlabel('Time, sec');


% Energy consumption
PowerConsmp = 0;

for i = 1:length(U_RW_Command_out)
    PowerConsmp = PowerConsmp + 0.5*norm(U_RW_Command_out(:,i))^2*dT;
    %PowerConsmp_BFMag = PowerConsmp_BFMag + 0.5*norm(U_RW_actual(:,i))^2*dT;
end

fprintf('Power Consumption of active FTC: %6.3f \n', PowerConsmp);

% figure(31);
% plot(XTime,h_lump_Err_out,'LineWidth',1.5);
% title('Disturbance Tracking Error')
% xlabel('Time, sec');
%
% figure(32);
% plot(XTime,h_lump_out,XTime,h_lump_hat_out,'LineWidth',1.5);
% title('Disturbance Tracking Error')
% xlabel('Time, sec');

