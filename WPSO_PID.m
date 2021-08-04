%带惯性权重的PSO-PID，本文所用代码源于教学代码修改而来
clc;
clear;
close all;
m=60;%粒子数
%dim=3;%粒子维数
dim=3;%粒子维数
vmax=[2,1,1];%速度最大值 kp,ki,kd
vmin=[-2,-1,-1];%速度最小值kp,ki,kd
xmax=[10,5,5];%位置最大值kp,ki,kd
xmin=[0,0,0];%位置最小值
tmax=100;%最大迭代次数
c1=2; 
c2=2; %学习因子
c3=2;
x=rand(m,dim);
v=rand(m,dim);
pfit=rand(1,m);
wmax=0.9; %惯性权重最大值 wpso
wmin=0.4; %惯性权重最小值 wpso
countermax=1;
gbest2=zeros(1,tmax);
pgt2=zeros(3,tmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%惯性权重线性递减调整
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%粒子群初始化
for counter=1:countermax
    for i=1:m    %位置初始化
        for j=1:dim
            x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j);       
        end
    end
    pb=x;%粒子最佳位置
    for i=1:m   %速度初始化
        for j=1:dim
             v(i,j)=rand()*(vmax(j)-vmin(j))+vmin(j);
        end
    end

    for i=1:m    %计算粒子初始位置适应度
        kp=x(i,1); %比例系数赋值
        ki=x(i,2); %积分系数赋值
        kd=x(i,3); %微分系数赋值
        [t,w1,u]=sim('shuixiang_pid',[0,10]);
       if min(u(:,2))<0
          pfit(i)=max(u(:,1))+c3*abs(min(u(:,2)));
       else
          pfit(i) =max(u(:,1));
       end
    end
    [gbest,index]=min(pfit);%gbest全局最小适应度
    pg=x(index,:);%全局最佳位置

    
    
    
    
    
    
    for t=1:tmax
         w=wmax-(wmax-wmin)*t/tmax;
        for i=1:m
            for j=1:dim             
                R1=c1*rand(1);
                R2=c2*rand(1);
                v(i,j)=w*v(i,j)+R1*(pb(i,j)-x(i,j))+R2*(pg(j)-x(i,j));
                if  v(i,j)>vmax(j)
                    v(i,j)=vmax(j);
                end
                if  v(i,j)<vmin(j)
                    v(i,j)=vmin(j);
                end        
                x(i,j)=x(i,j)+v(i,j);
                if x(i,j)>xmax(j)
                    x(i,j)=xmax(j);
                end
                if x(i,j)<xmin(j)
                    x(i,j)=xmin(j);
                end
            end
            kp=x(i,1);ki=x(i,2);kd=x(i,3); 
           [tt,w1,u]=sim('shuixiang_pid',[0,10]);      
%                pfitness=max(u(:,1));
          if min(u(:,2))<0
              pfitness=max(u(:,1))+c3*abs(min(u(:,2)));
           else
              pfitness=max(u(:,1));
           end
           if pfitness<pfit(i)
               pfit(i)=pfitness; %更新粒子最佳适应度
               pb(i,:)=x(i,:);%更新粒子最佳位置
           end
          if pfitness<gbest
               gbest=pfitness; %更全局最佳适应度
               pg=x(i,:);      %更新全局最佳位置
          end      
        end 
        gbest2(t)=gbest2(t)+gbest;    %全局最优位置适应度跟踪记录
        pgt2(1,t)=pgt2(1,t)+pg(1);      %参数跟踪记录
        pgt2(2,t)=pgt2(2,t)+pg(2);
        pgt2(3,t)=pgt2(3,t)+pg(3);
     
    end
end
 gbest2=gbest2/countermax;
 pgt2=pgt2/countermax;
 kp2=pgt2(1,tmax);
 ki2=pgt2(2,tmax);
  kd2=pgt2(3,tmax);%最优参数输出
 %%%%%%%%%%%%%%%%%%%%%%%最终优化参数统一仿真
% sim('PSO_PID',[0,10]);  
% chaotiao1=max(PID_out.signals.values(:,1));
% chaotiao1=100*(chaotiao1-1)/1;
% chaotiao2=max(PID_out.signals.values(:,2));
% chaotiao2=100*(chaotiao2-1)/1;
% chaotiao3=max(PID_out.signals.values(:,3));
% chaotiao3=100*(chaotiao3-1)/1;
% figure(1);
t=1:tmax;
plot(t,gbest2(t),'-x');
legend('标准PSO');
xlabel('迭代次数t');
ylabel('output目标函数值');、
