%% %%%%%%%%%%%% 循环求解不同摄动幅度和时滞下的 L2-gain %%%%%%%%%%%%
%% %%%%%%%%% pdt-switching qtd method th2 + 2023-10-4 %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
%% 区间切分
xox = 6;  dt = 1;  Nt = xox/dt; % 坐标轴1
yoy = 6;  dx = 1;  Nx = yoy/dx; % 坐标轴2
varphi_step = 0.5;  down_varphi = 0.5;  up_varphi = 3; 
nu_step = 0.25;  down_nu = 1;  up_nu = 2.25;  
varphi_slice = down_varphi:varphi_step:up_varphi;  nu_slice = down_nu:nu_step:up_nu;
varphi_array = zeros(Nt,Nx); nu_array = zeros(Nt,Nx);
gamma_max_th2 = zeros(Nt,Nx); % 用于存储所计算出的性能指标Gamma值的矩阵 便于绘图
for i = 1:Nt
    for j = 1:Nx
    varphi_array(i,j) = varphi_slice(i); 
    nu_array(i,j) = nu_slice(j); 
    end
end
%% 重复迭代计算gamma的过程（二重循环迭代）
for m = 1:Nt
    for n = 1:Nx
       %% pdt-switching 参数
        mode = 2; % 系统模态数
        rho = 0.33; % lemmma 引入的参数
        mu = 1.1; % 切换时刻 李雅普诺夫函数的上升率
        alpha = 0.4; % 指数衰减项
        sigma = 3; % 快切换区间长度a
        ht = 0.2; % 最小跨度
        delta_T = sigma/ht + 1;
        %tau = (delta_T*log(mu))/alpha - sigma % 最小下界
        tau = 0.9; % tau > (delta_T*log(mu))/alpha - sigma 
        span_counter = fix(tau/ht); % 跨度计数器
        
        varphi = varphi_array(m,n); % 非脆弱摄动
        nu = nu_array(m,n); % 非脆弱摄动
       %% 2-neuron CSTNN 参数
        A1 = [-1    0; 
               0   -1];  
        A2 = A1;  
        W1 = [2  -0.1; 
             -5   3.0];  
        W2 = [2  -0.1; 
             -5   4.5];  
        tilde_W1 = [-1.5  -0.1;
                    -0.2  -2.5];   
        tilde_W2 = [-1.5  -0.1;
                    -0.2  -4.0];  
        B1 = [-2 1; 0 -2];  B2 = [-2 0.1; 0 -2];  % Changed
        E1 = [1;    2];  E2 = [0.12; 0.1];
        C1 = [1 0; 0 1]; C2 = C1;
        D1 = [0.1 -0.1];  D2 = [0.6 -0.1]; % Changed
        L1 = 0.1;  L2 = L1;  
        Lipschitz_condition = 1;  
        d1 = Lipschitz_condition;  
        d2 = Lipschitz_condition; 
        G = blkdiag(d1, d2);
        H1 = [1.5; 0.6];  H2 = H1; % Changed
        J1 = varphi*[0.1  0.1];  J2 = J1;
        F1 = 0.6;  F2 = 0.8; % 非脆弱幅度
        %% 病态求解
        delta = 0;
        %% 循环定义决策变量
        gamma_v = sdpvar(1, 1,'symmetric'); % l2-gain 性能指标
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
                Xi_ij_44 = Li'*Li - gamma_v*eye(1);
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
                Xi_ij_44 = Li'*Li - gamma_v*eye(1);
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
        Constraint = [gamma_v >= 0, Q >= 0, R >= 0, epsilon >= 0, LMIA_Sets, LMIB_Sets];
        sdpsettings('solver', 'mosek', 'verbos', 0); % 设置求解器为 mosek，并打印少量信息
        reuslt = optimize(Constraint, gamma_v); % 自动求解的优化函数
        if reuslt.problem == 0 % problem = 0 代表求解成功
            gamma_v = sqrt(double(gamma_v));  % 将计算出的Gamma值存储在矩阵gammaset中，方便绘图
        else
        end
        %% gamma_sup
        disp('**********************************');
        disp('******* th3 gamma_max set ********');
        disp('**********************************'); 
        fengzi = alpha*mu^(delta_T);
        fengmu = alpha - (delta_T*log(mu))/(tau+sigma);
        gamma_max_th2(m,n) = gamma_v*sqrt(fengzi/fengmu)
    end
end
%% 性能指标gamma 三维图
figure('name', '性能指标 gamma_max_th2 三维网格图');
[x,y] = meshgrid(1:dt:xox,1:dx:yoy);
z_th2 = gamma_max_th2(1:Nt, 1:Nx);
h_th2 = surf(x, y, z_th2, 'LineWidth', 1.1, 'marker', '.', 'markersize', 12);
h_th2.FaceAlpha = 0.2;
h_th2.LineStyle = '-';
h_th2.FaceColor = 'interp';
xdata_th2 = get(h_th2,'xdata'); ydata_th2 = get(h_th2,'ydata'); 
set(h_th2,'xdata', down_varphi:varphi_step:up_varphi,'ydata', down_nu:nu_step:up_nu);  
xlabel('$\varphi$','Interpreter','latex','Fontsize',17);
ylabel('$\nu$','Interpreter','latex','Fontsize',17);
zlabel('$\gamma_{\diamond}$','Interpreter','latex','Fontsize',17);
colormap(hsv);
grid on;  grid minor;  box on;
axis([down_varphi up_varphi, down_nu up_nu, 2 4]);  
disp(['运行时间: ', num2str(toc)]); 