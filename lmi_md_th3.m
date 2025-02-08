%% %%%%%%%%%%%%%%%%%%%%%%%%%% 数值求解 %%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% pdt-switching md method th3 + 2023-10-4 %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
system_parameters
%% 循环定义决策变量
gamma = sdpvar(1, 1,'symmetric'); % l2-gain 性能指标
epsilon = sdpvar(1, 1,'symmetric'); % lemma1 待求解参数   
Q = sdpvar(2, 2,'symmetric');  
R1 = sdpvar(1, 1,'symmetric');  R2 = sdpvar(1, 1,'symmetric');  R = blkdiag(R1, R2);  
for i = 1 : mode
        eval(['P_i_', int2str(i), '= sdpvar(2, 2,''symmetric'');']);
        eval(['X_i_', int2str(i), '= sdpvar(2, 2,''full'');']); 
        eval(['Y_i_', int2str(i), '= sdpvar(2, 2,''full'');']); 
end
%% 循环定义 LMI (30)
for i = 1 : mode
    Ai = eval(['A', int2str(i)]);  Wi = eval(['W', int2str(i)]);  tilde_Wi = eval(['tilde_W', int2str(i)]);
    Ci = eval(['C', int2str(i)]);  Hi = eval(['H', int2str(i)]);  Ji = eval(['J', int2str(i)]);
    Bi = eval(['B', int2str(i)]);  Ei = eval(['E', int2str(i)]);  Di = eval(['D', int2str(i)]);  
    Li = eval(['L', int2str(i)]);
    P_i = eval(['P_i_', int2str(i)]);  Xi = eval(['X_i_', int2str(i)]);  Yi = eval(['Y_i_', int2str(i)]);
    Xi_i_11 = alpha*P_i + (P_i*Ai + Bi*Yi*Ci) + (P_i*Ai + Bi*Yi*Ci)' + epsilon*Ci'*Ji'*Ji*Ci + Di'*Di + G*R*G;
    Xi_i_14 = P_i*Ei + Di'*Li;
    Xi_i_16 = P_i*Bi - Bi*Xi + rho*(Yi*Ci)';
    Xi_i_44 = Li'*Li - gamma*eye(1);
    Xi_i_66 = -rho*(Xi + Xi');
    LMI = [Xi_i_11           P_i*Wi         P_i*tilde_Wi        Xi_i_14       P_i*Bi*Hi        Xi_i_16;
          (P_i*Wi)'          Q-R            zeros(2,2)          zeros(2,1)    zeros(2,1)       zeros(2,2);
          (P_i*tilde_Wi)'    zeros(2,2)    -exp(-alpha*nu)*Q    zeros(2,1)    zeros(2,1)       zeros(2,2);
          (Xi_i_14)'         zeros(1,2)     zeros(1,2)          Xi_i_44       zeros(1,1)       zeros(1,2);
          (P_i*Bi*Hi)'       zeros(1,2)     zeros(1,2)          zeros(1,1)   -epsilon*eye(1)   zeros(1,2);
          (Xi_i_16)'         zeros(2,2)     zeros(2,2)          zeros(2,1)    zeros(2,1)       Xi_i_66];
    eval(['LMIA_', int2str(i),'= [LMI <= delta*eye(size(LMI)), P_i >= delta*eye(size(P_i))];']); 
end
%% 循环定义 LMIB
for i = 1 : mode
   P_i = eval(['P_i_', int2str(i)]);
   P_l = eval(['P_i_', int2str(i+(-1)^(i-1))]);
   LMI_P = P_i - mu*P_l;
   eval(['LMIB_', int2str(i),'= LMI_P <= delta*eye(size(LMI_P));']);
end
%% 线性约束条件，需要同时满足
Constraint = [gamma >= 0, Q >= 0, R >= 0, epsilon >= 0, LMIA_1, LMIA_2, LMIB_1, LMIB_2];
sdpsettings('solver', 'mosek', 'verbos', 0); % 设置求解器为 mosek，并打印少量信息
reuslt = optimize(Constraint, gamma); % 自动求解的优化函数
if reuslt.problem == 0 % problem = 0 代表求解成功
    [primal,~] = check(Constraint); % 检查约束
    if min(primal) >= 0 && all(primal([1:2,4]) > 0) % 判断残差值
        disp('***********************************');
        disp('****Constraints are guaranteed*****');
        disp('**Primal residual is/are positive**');
        disp('***********************************');
        gamma = sqrt(value(gamma))
    else
        disp('********************************************');
        disp('**Warning: Primal residual is/are negative**');
        disp('********************************************');
        gamma = sqrt(value(gamma))
    end
else
    disp('**************************************');
    disp('****Constraints are not guaranteed****');
    disp('**************************************'); 
    reuslt.info
    yalmiperror(reuslt.problem) % 打印出错信息
    check(Constraint); % 检查残余值
end
%% 求解控制增益
for i = 1 : mode
    Xi = eval(['X_i_', int2str(i)]);  Yi = eval(['Y_i_', int2str(i)]);
    eval(['K', int2str(i), '= inv(value(Xi))*value(Yi);']);
end
disp('**********************************');
disp('********** 控制增益 Ki ***********');
disp('**********************************'); 
% K1, K2
%% gamma_sup
disp('**********************************');
disp('********* th3 gamma_max **********');
disp('**********************************'); 
fengzi = alpha*mu^(delta_T);
fengmu = alpha - (delta_T*log(mu))/(tau+sigma);
gamma_max_th3 = gamma*sqrt(fengzi/fengmu)
span_counter = span_counter
disp(['运行时间: ', num2str(toc)]); 