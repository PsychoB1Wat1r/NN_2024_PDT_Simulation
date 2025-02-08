%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%% pdt-switching 系统参数 + 2023-10-4 %%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode = 2; % 系统模态数
rho = 0.33; % lemmma 引入的参数
mu = 1.1; % 切换时刻 李雅普诺夫函数的上升率
alpha = 0.4; % 指数衰减项

sigma = 3; % 快切换区间长度
ht = 0.2; % 最小跨度
delta_T = sigma/ht + 1;
%tau = (delta_T*log(mu))/alpha - sigma % 最小下界
tau = 0.9; % tau > (delta_T*log(mu))/alpha - sigma 
span_counter = fix(tau/ht); % 跨度计数器

nu = 1; % 神经元时滞  % changed in [1.0, 1.2]
varphi = 1.5; % 非脆弱摄动  % changed in [0.5, 2.5]
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
F1 = 0.6;  F2 = 0.8; 
%% 病态求解
delta = 0;