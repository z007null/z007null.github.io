%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%利用改进的粒子群算法完成LORENZ混沌系统的参数估计（三个参数未知）
%k1,k2,k3 误差系数 
%日期：2016.6.20 
%修订：2018.5.9
%地点：贵阳小河
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clc;
clc;
clear;
close all;
m=30;%粒子数
dim=3;%粒子维数
xmax=[11,30,3];%位置最大值
xmin=[9,20,2];%位置最小值
vmax=[0.2,1,0.1];%速度最大值 
vmin=[-0.2,-1,-0.1];%速度最小值
tmax=30;%最大迭代次数

nummax=1;
c1=2; 
c2=2; %学习因子

wmax=0.9; %惯性权重最大值 wpso
wmin=0.4; %惯性权重最小值 wpso


aa2=zeros(1,tmax);
bb2=zeros(1,tmax);
cc2=zeros(1,tmax);
g_fitness2=zeros(1,tmax);
pgnum2=zeros(nummax,3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSO
for num=1:nummax
    %%%位置初始化,计算粒子适应度
        for i=1:m
            for j=1:dim    
                    x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j); 
                   v(i,j)=rand()*(vmax(j)-vmin(j))+vmin(j);      
            end 
             k1=x(i,1);k2=x(i,2);k3=x(i,3);
             [n,w1,u]=sim('lorenz1',0.1);%无噪声
%              [n,w1,u]=sim('lorenz',5);%有噪声
             fit(i)=max(u(:,1));  
        end

        [gbest,index]=min(fit);%gbest全局最小适应度
        pg=x(index,:);%全局最佳位置
        pb=x;%个体最佳位置
    for t=1:tmax
        w=wmax-(wmax-wmin)*t/tmax;
     for i=1:m
        for j=1:dim
                R1=c1*rand();
                R2=c2*rand();  
                v(i,j)=w*v(i,j)+R1*(pb(i,j)-x(i,j))+R2*(pg(j)-x(i,j));  
                 if v(i,j)>vmax(j)
                             v(i,j)=vmax(j);
                end
                if v(i,j)<vmin(j)
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
         k1=x(i,1);k2=x(i,2);k3=x(i,3);
             [n,w1,u]=sim('lorenz1',0.1);%无噪声
%             [n,w1,u]=sim('lorenz',5);%有噪声
         pfitness=max(u(:,1));   
        if pfitness<fit(i)
            fit(i)=pfitness;
            pb(i,:)=x(i,:);
        end
       if pfitness<gbest
            gbest=pfitness;
            pg=x(i,:);
       end 
    end  
        aa2(t)=aa2(t)+pg(1);%估计参数a跟踪
        bb2(t)=bb2(t)+pg(2);%估计参数b跟踪
        cc2(t)= cc2(t)+pg(3);%估计参数c跟踪
        g_fitness2(t)=g_fitness2(t)+gbest;%统计适应度值进化曲线  
    end    
    pgnum2(num,:)=pg;    %统计每次试验的寻优结果
end


a_mean2=mean(pgnum2(:,1));%a参数的平均估计值
b_mean2=mean(pgnum2(:,2));%b参数的平均估计值
c_mean2=mean(pgnum2(:,3));%c参数的平均估计值


figure(1);
i=1:tmax;
plot(i,log10(g_fitness2(i)/nummax),'-r');

xlabel('t');
ylabel('v_f');

figure(2);
plot(i,aa2(i)/nummax,'-r');

xlabel('t');
ylabel('a_1'); 

figure(3);
plot(i,bb2(i)/nummax,'-r');

xlabel('t');
ylabel('b_1'); 

figure(4);
plot(i,cc2(i)/nummax,'-r');

xlabel('t');
ylabel('c_1'); 
  