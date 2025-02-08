%% %%%%%%%%%%%%%%%%%%%%%% 三维比较分析网格图 %%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%% L2-gain 的比较 + 2023-10-4 %%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
% axis off; % 清除坐标轴属性
%% 区间切分
xox = 6;  dt = 1;  Nt = xox/dt; % 坐标轴1
yoy = 6;  dx = 1;  Nx = yoy/dx; % 坐标轴2
varphi_step = 0.5;  down_varphi = 0.5;  up_varphi = 3; 
nu_step = 0.25;  down_nu = 1;  up_nu = 2.25;   
varphi_slice = down_varphi:varphi_step:up_varphi;  nu_slice = down_nu:nu_step:up_nu;
varphi = zeros(Nt,Nx); nu = zeros(Nt,Nx);
gamma_max_th2 = zeros(Nt,Nx); % 用于存储所计算出的性能指标Gamma值的矩阵 便于绘图
gamma_max_th3 = zeros(Nt,Nx); % 用于存储所计算出的性能指标Gamma值的矩阵 便于绘图
load('gamma_max_th2.mat'); load('gamma_max_th3.mat'); % 导入数据
% gamma_max_th2 = rot90(gamma_max_th2, 1);  gamma_max_th3 = rot90(gamma_max_th3, 1);
%% 性能指标gamma 三维网格图
figure('name', '性能指标 gamma_max 三维网格图');
[x,y] = meshgrid(1:dt:xox,1:dx:yoy);
z_th3 = gamma_max_th3(1:Nt, 1:Nx);
h_th3 = mesh(x, y, z_th3,  'LineWidth', 1.2, 'marker', 'd', 'markersize', 7, 'MarkerFaceColor', 'r');
h_th3.FaceAlpha = 0.5; % 透明度设置 0~1
h_th3.LineStyle = '-'; % 曲面线条样式
h_th3.EdgeColor = 'r'; % 曲面 线颜色
h_th3.FaceColor = 'none'; % 曲面 面颜色
xdata_th3 = get(h_th3,'xdata'); ydata_th3 = get(h_th3,'ydata'); 
set(h_th3,'xdata', down_varphi:varphi_step:up_varphi, 'ydata', down_nu:nu_step:up_nu);
hold on
z_th2 = gamma_max_th2(1:Nt, 1:Nx);
h_th2 = mesh(x, y, z_th2, 'LineWidth', 1.2, 'marker', 'o', 'markersize', 7);
h_th2.FaceAlpha = 1; % 透明度设置 0~1
h_th2.LineStyle = '-.'; % 曲面线条样式
h_th2.EdgeColor = '#87CEEB'; % 曲面 线颜色 
h_th2.FaceColor = 'none'; % 曲面 面颜色 
xdata_th2 = get(h_th2,'xdata'); ydat_th2 = get(h_th2,'ydata'); 
set(h_th2,'xdata', down_varphi:varphi_step:up_varphi, 'ydata', down_nu:nu_step:up_nu);
set(gca, 'FontSize', 11); % 坐标轴字体大小
xlabel('$\varphi$','Interpreter','latex','Fontsize',15);
ylabel('$\nu$','Interpreter','latex','Fontsize',15);
zlabel('$\gamma_{\diamond}$','Interpreter','latex','Fontsize',16);
axis([down_varphi up_varphi, down_nu up_nu-0.25, 2 4]);
axis ij;  grid on;  grid minor;  box on; 
% grid off;  box off;  
disp(['运行时间: ', num2str(toc)]); 