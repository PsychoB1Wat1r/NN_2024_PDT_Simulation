%% %%%%%%%%%%%%%%%%%%%%%%%%% 混沌行为 %%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%% 开环系统相平面图 + 2023-10-4 %%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
%% 仿真绘图参数
time = 400; % 仿真时长
dt = 0.01; % 仿真步长
N = time/dt; % % 系统迭代次数
sigma = 3; % 快切换区间长度
mode = 2; % 系统模态数
alpha = 0.4; % 指数衰减项
mu = 1.1; % 切换时刻 李雅普诺夫函数的上升率
f = 5;
nu = 1; % 神经元时滞
Ntau = nu/dt;
varphi = 1.5; % 非脆弱摄动
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
%% 控制增益 开环控制增益置零
K1 = [0   0
      0   0];
K2 = [0   0
      0   0];
%% 迭代过程
x_array = zeros(2,N+Ntau); % 系统状态存储矩阵
g_x = zeros(2,N+Ntau); % 激活函数存储矩阵
g_delay_x = zeros(2,N+Ntau); % 延迟 激活函数存储矩阵
w_array = zeros(1,N+Ntau); % 噪声信号存储矩阵
%% 系统初值
for i = 1 : Ntau
    x_array(:,i) = [0.3; -0.3];
end
x_array(:,Ntau+1) = [0.3; -0.3];
Mode_array = [];  Phi_array = [];
[Mode,Phi] = ssp_pdt_qtd_mode(time+nu,dt,mode,sigma,alpha,mu,f);
for n = Ntau+1 : N+Ntau
      temp_mode = Mode(n);
      if temp_mode == 1
        A = A1;  B = B1;  C = C1;  W = W1;  tilde_W = tilde_W1;
        H = H1;  E = E1;  D = D1;  L = L1;  J = J1;
        K = K1;
      else 
        A = A2;  B = B2;  C = C2;  W = W2;  tilde_W = tilde_W2;
        H = H2;  E = E2;  D = D2;  L = L2;  J = J2;
        K = K2;
      end 
      g_x(:,n) = [tanh(x_array(1,n));tanh(x_array(2,n))];
      g_delay_x(:,n) = [tanh(x_array(1,n-Ntau));tanh(x_array(2,n-Ntau))];
      w_array(:,n) = 2.5*exp(-0.25*n*dt)*sin(2.5*pi*n*dt);
      x_array(:,n+1) = x_array(:,n) + dt*( A*x_array(:,n) + W*g_x(:,n) + tilde_W*g_delay_x(:,n) + E*w_array(:,n) );
end
%% 相平面图
figure('name', '相平面图');
plot(x_array(2, 1000:time/dt), x_array(1, 1000:time/dt), 'color', '#11009E', 'linewidth', 1.1);
set(gca, 'FontSize', 11, 'Linewidth', 1);
xlabel('$x_{2}(t)$', 'Interpreter', 'latex', 'fontsize', 15);
ylabel('$x_{1}(t)$', 'Interpreter', 'latex', 'fontsize', 15);
grid on;  grid minor;  box on;
% axis([-6 6, -1.5 1]);
%% 开环系统演化轨迹绘制
figure('name', '开环系统演化轨迹绘制');
plot(0:dt:time-dt, x_array(1, 1:time/dt),'-b','linewidth',1);
hold on
plot(0:dt:time-dt, x_array(2, 1:time/dt),'-r','linewidth',1);
set(gca, 'FontSize', 11, 'Linewidth', 1);
xlabel('$t$', 'Interpreter', 'latex', 'fontsize', 15);
ylabel('$x(t)$', 'Interpreter', 'latex', 'fontsize', 15);
h2 = legend('$x_{1}(t)$', '$x_{2}(t)$', 'Interpreter', 'latex', 'fontsize', 15);
set(h2,'Orientation', 'horizon', 'Box', 'on');
grid on;  grid minor;  box on;
disp(['运行时间: ', num2str(toc)]); 