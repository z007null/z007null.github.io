%利用改进的粒子群与惯性权重线性递减的粒子群分别对Griewank函数
%进行优化。粒子数：60；维数：10；最大迭代次数：500；分别测试50次取平均值。
%时间：2017.6.17
      
clc;
clear;
m=60;%粒子数
dim=10;%粒子维数
vmax=[60,60,60,60,60,60,60,60,60,60];%速度最大值 
vmin=[-60,-60,-60,-60,-60,-60,-60,-60,-60,-60];%速度最小值
xmax=[600,600,600,600,600,600,600,600,600,600];%位置最大值
xmin=[-600,-600,-600,-600,-600,-600,-600,-600,-600,-600];%位置最小值
tmax=300;%最大迭代次数
kmax=500;
c1=2; 
c2=2; %学习因子
x=rand(m,dim);
v=rand(m,dim);
pfit=rand(1,m);
flag=0;%判断进化停滞标志
wmax=0.9; %惯性权重最大值 wpso
wmin=0.4; %惯性权重最小值 wpso
gbest1= zeros(1,tmax);
gbest2= zeros(1,tmax);
position1=zeros(tmax,dim);%粒子位置跟踪
position2=zeros(tmax,dim);%粒子位置跟踪
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin
nummax=1; %测试次数
testbest1=zeros(1,nummax); %每次试验的最优值
testbest2=zeros(1,nummax); %每次试验的最优值
%%%%%%%%%%%%%%%%%%%%混沌粒子群
for number=1:nummax      
    for i=1:m    %位置混沌初始化  m1个粒子中选择最优的m个粒子为初始粒子
        result1=0;
        result2=1;
        for j=1:dim    
            x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j);
            result1=result1+(x(i,j)^2)/4000;
            result2=result2*cos(x(i,j)/j^0.5);              
        end     
        pfit(i)=1+result1-result2; 
    end    
           
    pb=x;%粒子初始最佳位置赋值
    for i=1:m   %速度初始化
        for j=1:dim
             v(i,j)=rand()*(vmax(j)-vmin(j))+vmin(j);
        end
    end 
   [gbest,index]=min(pfit);%gbest全局最小适应度
   gbestbefore=0;%上次全局最小值
   pg=x(index,:);%全局最佳位置
    for t=1:tmax
         w=wmax-(wmax-wmin)*t/tmax+(0.5-rand())*0.8;
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
                         x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j); 
                end
                if x(i,j)<xmin(j)
                        x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j); 
                end
            end
            result1=0;
            result2=1;   
            for j=1:dim
                result1=result1+(x(i,j)^2)/4000;
                result2=result2*cos(x(i,j)/j^0.5);          
            end    
            pfitness=1+result1-result2;
            if pfitness<pfit(i)
                pfit(i)=pfitness;
                pb(i,:)=x(i,:);
           end
            if pfitness<gbest
                gbest=pfitness;
                pg=x(i,:);
            end      
        end         
        cflag=abs((gbest-gbestbefore)/gbest);
        if cflag<0.01
            flag=flag+1;
            if flag>10
                    for k=1:kmax
                       for i=1:m  
                        u=pg;
                        if rand()<0.95
                            x(i,:)=u+(0.5-rand())*u*abs(log10((k/kmax)^2))/exp(50*(k/kmax)^2);
                            for j=1:dim
                                if x(i,j)>xmax(j)||x(i,j)<xmin(j)
                                     x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j);
                                end
                            end
                        else
                             x(i,:)=rand()*(xmax-xmin)+xmin;
                        end
                        result1=0;
                        result2=1;   
                        for j=1:dim
                            result1=result1+(x(i,j)^2)/4000;
                            result2=result2*cos(x(i,j)/j^0.5);          
                        end    
                        pfitness=1+result1-result2;
                        if pfitness<pfit(i)
                             pfit(i)=pfitness;
                             pb(i,:)= x(i,:);
                        end
                        if pfitness<gbest
                            flag=0;  
                            gbest=pfitness;
                            pg=x(i,:);                                             
                        end    
                       end                       
                    end
            end
        else
            flag=0;
        end      
        gbestbefore=gbest;
        gbest1(t)=gbest1(t)+ gbest;
    end   
    testbest1(number)=gbest;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%惯性权重线性递减方式
for number=1:nummax
   for i=1:m
        result1=0;
        result2=1;
        for j=1:dim    
            x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j);
            result1=result1+(x(i,j)^2)/4000;
            result2=result2*cos(x(i,j)/j^0.5);              
        end     
        pfit(i)=1+result1-result2;       
   end  
   pb=x;
   for i=1:m
       for j=1:dim
           v(i,j)=rand()*(vmax(j)-vmin(j))+vmin(j);
       end
   end
   [gbest,index]=min(pfit);%gbest全局最小适应度
    pg=x(index,:);%全局最佳位置
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
            result1=0;
            result2=1;   
            for j=1:dim
                result1=result1+(x(i,j)^2)/4000;
                result2=result2*cos(x(i,j)/j^0.5);          
            end    
            pfitness=1+result1-result2;
            if pfitness<pfit(i)
                pfit(i)=pfitness;
                pb(i,:)=x(i,:);
           end
            if pfitness<gbest
                gbest=pfitness;
                pg=x(i,:);
            end      
        end         
        gbest2(t)=gbest2(t)+gbest;
        position2(t,:)=position2(t,:)+pg;
    end  
    testbest2(number)=gbest;
end
testbestmax1=max(testbest1);%最大值
testbestmin1=min(testbest1);%最小值
testbestmean1=mean(testbest1);%平均值
%%%%%%%%%%%%%%%计算标准差
fangcha1=0;
for number=1:nummax
    fangcha1=fangcha1+(testbest1(number)-testbestmean1)^2;
end
fangcha1=sqrt(fangcha1/nummax);

testbestmax2=max(testbest2);%最大值
testbestmin2=min(testbest2);%最小值
testbestmean2=mean(testbest2);%平均值
%%%%%%%%%%%%%%%计算标准差
fangcha2=0;
for number=1:nummax
    fangcha2=fangcha2+(testbest2(number)-testbestmean2)^2;
end
fangcha2=sqrt(fangcha2/nummax);

t=1:tmax;
gbest1(t)=gbest1(t)/nummax;
gbest2(t)=gbest2(t)/nummax;
position1(t)=position1(t)/nummax;
position2(t)=position2(t)/nummax;
figure(1);
plot(t,log10(gbest1(t)),'-r',t,log10(gbest2(t)),'--b');
% plot(t,gbest1(t),'-r',t,gbest2(t),'--b');
figure(1);
legend('改进PSO','标准PSO');
xlabel('进化次数');
ylabel('适应值的对数');
figure(2);
plot(1:nummax,testbest1(1:nummax),'-d',1:nummax,testbest2(1:nummax),'-o');
legend('改进PSO','标准PSO');
xlabel('实验次数');
ylabel('适应值');
% figure(4);
% plot(t,position1(t,:));
% figure(5);
% plot(t,position2(t,:));
% grid on;






 

