function [Mode_array,Phi_array] = ssp_pdt_qtd_mode(time,dt,mode,T,alpha,mu,f)
rand('state',20);

sum = 0;
stage = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ɵ������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time=15;
% dt=0.02;
% mode=2;
% 
% T=3;%���л����������ʱ��
% alpha=0.9;
% mu=1.01;
% f=10;%����л�Ƶ��

ht = 0.2;
f = 1/ht;
tau = ((T*f+1)*log(mu))/alpha-T; % ���л�����С����ʱ��

T_floor = 0.4*T; % ���л��������С�½�
tau_up = 2.2*tau; % ���л����������Ͻ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MSP = 1/f; % MSP��ʾ��С�л����ڣ�
MAX_Phi = fix(tau*f); % QTD������Phi


while sum < time/dt

    tau_temp = ceil((tau + (tau_up-tau)*rand())/dt);
    T_temp = fix((T - (T-T_floor)*rand())/dt);
    if sum + tau_temp > time/dt
        tau_temp = time/dt - sum;
        T_temp = 0;
    elseif sum + tau_temp + T_temp > time/dt
        T_temp = time/dt - (sum+tau_temp);         
    end  
    sum = sum + tau_temp + T_temp;
    tau_stage(stage) = tau_temp;
    T_stage(stage)=T_temp;
    stage=stage+1;
end 
    
stage=size(tau_stage,2);



for i=1:stage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���л���ģ̬ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mode_tau_i=round(rand()*(mode-1) + 1);
    
    if i > 1
        while mode_tau_i == peri_switch_mode
            mode_tau_i=round(rand()*(mode-1) + 1);
        end
    end
    tau_i=tau_stage(i);
    T_i=T_stage(i);
    peri_switch_mode=mode_tau_i;
    
    eval(strcat('tau_mode_array',int2str(i),'=ones(1,tau_i);'));
    eval(strcat('tau_Phi_array',int2str(i),'=ones(1,tau_i);'));
    eval(strcat('T_mode_array',int2str(i),'=zeros(1,T_i);'));
    eval(strcat('T_Phi_array',int2str(i),'=ones(1,T_i);'));
    
    for j=1:tau_i          
        eval(strcat('tau_mode_array',int2str(i),'(j)=mode_tau_i;'));% mode_array(j)=mode_tau_i;
    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���л���QTD���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for j=1:tau_i
        if  j <= MAX_Phi/(f*dt)
            for k=1:MAX_Phi
                blow_bound=(MSP/dt)*(k-1)+1;
                up_bound=(MSP/dt)*k;
                
                if blow_bound <=j && j <=up_bound 
                 
                    eval(strcat('tau_Phi_array',int2str(i),'(j)=k-1;'));%Phi_array_stage(j)=k-1;
                end
            
            end
        else
                    eval(strcat('tau_Phi_array',int2str(i),'(j)=MAX_Phi;'));% Phi_array(j)=MAX_Phi
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���л����� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sum_T_ij=0;
    j=0;
    number1=0;
    number2=0;
%     peri_switch_mode=mode_tau_i;
    while j < T_i
        
        T_ij=fix((1/f+0.3*rand()*(tau-1/f))/dt);%�ϸ���˵�ǣ�T_ij=fix((1/f+rand()*(tau-1/f))/dt)���ĵ���Ϊ�˷�ֹT_ij��tau̫�ӽ�
        if T_i - (sum_T_ij + T_ij) < (MSP/dt)
            T_ij=T_i-sum_T_ij;     
        end
        sum_T_ij=sum_T_ij+T_ij;
%         MAX_Phi_T_ij=fix(T_ij*f);
        
        mode_T_ij=round(rand()*(mode-1) + 1);
        
        while mode_T_ij ==peri_switch_mode
            mode_T_ij=round(rand()*(mode-1) + 1);
        end
        
        for m=(number1+1):(number1+T_ij) 
            
            eval(strcat('T_mode_array',int2str(i),'(m)=mode_T_ij;'));
            
            if m == number1+T_ij
                number1=number1+T_ij;
            end
            
        end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���л������ڵ�QTD���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        for m=(number2+1):(number2+T_ij)
     
                for k=1:MAX_Phi
                    blow_bound=(MSP/dt)*(k-1)+1;
                    up_bound=(MSP/dt)*k;
                
                    if blow_bound <=(m-number2) && (m-number2) <=up_bound 
                 
                        eval(strcat('T_Phi_array',int2str(i),'(m)=k-1;'));%Phi_array_stage(j)=k-1;
                    end
                    
                end
                
                if m == number2+T_ij
                        number2=number2+T_ij;
                end
            
        end
        
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        j=number1;
        peri_switch_mode=mode_T_ij;
        
    end

    
end

Mode_array=[];
Phi_array=[];


for i=1:stage

    eval(['Mode_array=[Mode_array,tau_mode_array',num2str(i),'];']); 
    eval(['Mode_array=[Mode_array,T_mode_array',num2str(i),'];']);
   
    eval(['Phi_array=[Phi_array,tau_Phi_array',num2str(i),'];']);
    eval(['Phi_array=[Phi_array,T_Phi_array',num2str(i),'];']);
end


end

