%% %%%%%%%%%%%%%%%%%%%%%%%%%% 数值求解 %%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% pdt-switching qtd method th2 + 2023-10-4 %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
system_parameters
%% 循环定义决策变量
gamma = sdpvar(1, 1,'symmetric'); % l2-gain 性能指标
epsilon = sdpvar(1, 1,'symmetric'); % lemma1 待求解参数 
Q = sdpvar(2, 2,'symmetric');  
R11 = sdpvar(1, 1,'symmetric');   R22 = sdpvar(1, 1,'symmetric');  R = blkdiag(R11, R22);  
for i = 1 : mode
    for j = 0 : span_counter
        eval(['P_ij_', int2str(i), int2str(j), '= sdpvar(2, 2,''symmetric'');']);
        eval(['X_ij_', int2str(i), int2str(j), '= sdpvar(2, 2,''full'');']); 
        eval(['Y_ij_', int2str(i), int2str(j), '= sdpvar(2, 2,''full'');']); 
    end
end
%% 循环定义 LMIs (19)-(20)
for i = 1 : mode
    Ai = eval(['A', int2str(i)]);  Wi = eval(['W', int2str(i)]);  tilde_Wi = eval(['tilde_W', int2str(i)]);
    Ci = eval(['C', int2str(i)]);  Hi = eval(['H', int2str(i)]);  Ji = eval(['J', int2str(i)]);
    Bi = eval(['B', int2str(i)]);  Ei = eval(['E', int2str(i)]);  Di = eval(['D', int2str(i)]);  
    Li = eval(['L', int2str(i)]);
    for j = 0 : span_counter-1
        X_ij = eval(['X_ij_', int2str(i), int2str(j)]);
        Y_ij = eval(['Y_ij_', int2str(i), int2str(j)]);  
        P_ij = eval(['P_ij_', int2str(i), int2str(j)]); 
        bar_P_ij = eval(['P_ij_', int2str(i), int2str(j+1)]); 
        
        Xi_1ij_11 = alpha*P_ij + (P_ij*Ai + Bi*Y_ij*Ci) + (P_ij*Ai + Bi*Y_ij*Ci)' + epsilon*Ci'*Ji'*Ji*Ci + (bar_P_ij - P_ij)/ht + Di'*Di + G*R*G;
        Xi_1ij_14 = P_ij*Ei + Di'*Li;
        Xi_1ij_15 = P_ij*Bi*Hi;
        Xi_1ij_16 = P_ij*Bi - Bi*X_ij + rho*(Y_ij*Ci)';
        
        Xi_2ij_11 = alpha*bar_P_ij + (bar_P_ij*Ai + Bi*Y_ij*Ci) + (bar_P_ij*Ai + Bi*Y_ij*Ci)' + epsilon*Ci'*Ji'*Ji*Ci + (bar_P_ij - P_ij)/ht +  Di'*Di + G*R*G;
        Xi_2ij_14 = bar_P_ij*Ei + Di'*Li;
        Xi_2ij_15 = bar_P_ij*Bi*Hi;
        Xi_2ij_16 = bar_P_ij*Bi - Bi*X_ij + rho*(Y_ij*Ci)';
        Xi_ij_44 = Li'*Li - gamma*eye(1);
        Xi_ij_66 = -rho*(X_ij + X_ij');
        
        LMI_1 = [Xi_1ij_11          P_ij*Wi         P_ij*tilde_Wi       Xi_1ij_14      Xi_1ij_15        Xi_1ij_16;
                (P_ij*Wi)'          Q-R             zeros(2,2)          zeros(2,1)     zeros(2,1)       zeros(2,2);
                (P_ij*tilde_Wi)'    zeros(2,2)     -exp(-alpha*nu)*Q    zeros(2,1)     zeros(2,1)       zeros(2,2);
                (Xi_1ij_14)'        zeros(1,2)      zeros(1,2)          Xi_ij_44       zeros(1,1)       zeros(1,2);
                (Xi_1ij_15)'        zeros(1,2)      zeros(1,2)          zeros(1,1)    -epsilon*eye(1)   zeros(1,2);
                (Xi_1ij_16)'        zeros(2,2)      zeros(2,2)          zeros(2,1)     zeros(2,1)       Xi_ij_66];
            
        LMI_2 = [Xi_2ij_11              bar_P_ij*Wi     bar_P_ij*tilde_Wi    Xi_2ij_14      Xi_2ij_15         Xi_2ij_16;
                (bar_P_ij*Wi)'          Q-R             zeros(2,2)           zeros(2,1)     zeros(2,1)        zeros(2,2);
                (bar_P_ij*tilde_Wi)'    zeros(2,2)     -exp(-alpha*nu)*Q     zeros(2,1)     zeros(2,1)        zeros(2,2);
                (Xi_2ij_14)'            zeros(1,2)      zeros(1,2)           Xi_ij_44       zeros(1,1)        zeros(1,2);
                (Xi_2ij_15)'            zeros(1,2)      zeros(1,2)           zeros(1,1)    -epsilon*eye(1)    zeros(1,2);
                (Xi_2ij_16)'            zeros(2,2)      zeros(2,2)           zeros(2,1)     zeros(2,1)        Xi_ij_66];
        eval(['LMIA_', int2str(i), int2str(j), '= [LMI_1 <= delta*eye(size(LMI_1)), LMI_2 <= delta*eye(size(LMI_2)), P_ij >= delta*eye(size(P_ij)), bar_P_ij >= delta*eye(size(bar_P_ij))];']); 
    end
end
%% 循环定义 LMIs (21)-(22)
for i = 1 : mode
    Ai = eval(['A', int2str(i)]);  Wi = eval(['W', int2str(i)]);  tilde_Wi = eval(['tilde_W', int2str(i)]);
    Ci = eval(['C', int2str(i)]);  Hi = eval(['H', int2str(i)]);  Ji = eval(['J', int2str(i)]);
    Bi = eval(['B', int2str(i)]);  Ei = eval(['E', int2str(i)]);  Di = eval(['D', int2str(i)]);  
    Li = eval(['L', int2str(i)]);  P_i0 = eval(['P_ij_', int2str(i), int2str(0)]);
    for j = 1 : span_counter
        X_ij = eval(['X_ij_', int2str(i), int2str(j)]);
        Y_ij = eval(['Y_ij_', int2str(i), int2str(j)]);
        P_lj = eval(['P_ij_', int2str(i+(-1)^(i-1)), int2str(j)]);
        bar_Pim = eval(['P_ij_', int2str(i), int2str(j)]);
        
        Xi_3ij_11 = alpha*bar_Pim + (bar_Pim*Ai + Bi*Y_ij*Ci) + (bar_Pim*Ai + Bi*Y_ij*Ci)' + epsilon*Ci'*Ji'*Ji*Ci + Di'*Di + G*R*G;
        Xi_3ij_14 = bar_Pim*Ei + Di'*Li;
        Xi_3ij_15 = bar_Pim*Bi*Hi;
        Xi_3ij_16 = bar_Pim*Bi - Bi*X_ij + rho*(Y_ij*Ci)';
        Xi_ij_44 = Li'*Li - gamma*eye(1);
        Xi_ij_66 = -rho*(X_ij + X_ij');
        
        LMI_3 = [Xi_3ij_11            bar_Pim*Wi     bar_Pim*tilde_Wi    Xi_3ij_14     Xi_3ij_15        Xi_3ij_16;
                (bar_Pim*Wi)'         Q-R            zeros(2,2)          zeros(2,1)    zeros(2,1)       zeros(2,2);
                (bar_Pim*tilde_Wi)'   zeros(2,2)    -exp(-alpha*nu)*Q    zeros(2,1)    zeros(2,1)       zeros(2,2);
                (Xi_3ij_14)'          zeros(1,2)     zeros(1,2)          Xi_ij_44      zeros(1,1)       zeros(1,2);
                (Xi_3ij_15)'          zeros(1,2)     zeros(1,2)          zeros(1,1)   -epsilon*eye(1)   zeros(1,2);
                (Xi_3ij_16)'          zeros(2,2)     zeros(2,2)          zeros(2,1)    zeros(2,1)       Xi_ij_66];
        LMI_Pl = P_i0 - mu*P_lj;  
        eval(['LMIB_', int2str(i), int2str(j), '= [LMI_3 <= delta*eye(size(LMI_3)), bar_Pim >= delta*eye(size(bar_Pim)), LMI_Pl <= delta*eye(size(LMI_Pl))];']); 
    end
end
%% 线性矩阵不等式拼接
LMIA_Sets = ''; % 初始化空字符，用于拼接 LMIs (19)-(20)
for i = 1 : mode
    for j = 0 : span_counter-1
        strA_temp = eval(['LMIA_', int2str(i), int2str(j), ',']);
        LMIA_Sets = [LMIA_Sets, strA_temp]; % 线性矩阵不等式拼接
    end
end
LMIB_Sets = ''; % 初始化空字符，用于拼接  LMIs (21)-(22)
for i = 1 : mode
    for j = 1 : span_counter
        strB_temp = eval(['LMIB_', int2str(i), int2str(j), ',']);
        LMIB_Sets = [LMIB_Sets, strB_temp]; % 线性矩阵不等式拼接
    end
end
%% 线性约束条件，需要同时满足
Constraint = [gamma >= 0, Q >= 0, R >= 0, epsilon >= 0, LMIA_Sets, LMIB_Sets];
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
    for j = 0 : span_counter
        X_ij = eval(['X_ij_', int2str(i), int2str(j)]);  Y_ij = eval(['Y_ij_', int2str(i), int2str(j)]);
        eval(['K', int2str(i), int2str(j), '= inv(value(X_ij))*value(Y_ij);']);
    end
end
disp('**********************************');
disp('********** 控制增益 Kij ***********');
disp('**********************************'); 
K10, K11, K12, K13, K14
K20, K21, K22, K23, K24
%% gamma_sup
disp('**********************************');
disp('********* th2 gamma_max **********');
disp('**********************************'); 
fengzi = alpha*mu^(delta_T);
fengmu = alpha - (delta_T*log(mu))/(tau+sigma);
gamma_max_th2 = gamma*sqrt(fengzi/fengmu)
span_counter = span_counter
disp(['运行时间: ', num2str(toc)]); 