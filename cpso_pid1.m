%PSO-PIDʵ��Աȣ���������Ⱥ��LORENZ������׼����Ⱥ
%ʱ�䣺2016.7.16
clc;
clear;
close all;
m=20;%������
m1=50;
dim=3;%����ά��
vmax=[2,1,1];%�ٶ����ֵ kp,ki,kd
vmin=[-2,-1,-1];%�ٶ���Сֵkp,ki,kd
xmax=[10,5,5];%λ�����ֵkp,ki,kd
xmin=[0,0,0];%λ����Сֵ
tmax=50;%����������
c1=2; 
c2=2; %ѧϰ����
c3=2;
%%%
x=rand(m,dim);
v=rand(m,dim);
pfit=rand(1,m);
flag=0;
wmax=0.9; %����Ȩ�����ֵ wpso
wmin=0.4; %����Ȩ����Сֵ wpso
%%%
countermax=1;
gbest1=zeros(1,tmax);
pgt1=zeros(3,tmax);
gbest2=zeros(1,tmax);
pgt2=zeros(3,tmax);
%%%%%%%%%%%Lorenz����

for counter=1:countermax
    k1=rand();
    k2=rand();
    k3=rand();
for i=1:m1    %λ�ó�ʼ��   
        if k1<=0.4       
           k1=k1/0.4;
        else
           k1=(1-k1)/(1-0.4);
        end  
        if k2<=0.4       
           k2=k2/0.4;
        else
           k2=(1-k2)/(1-0.4);
        end  
         if k3<=0.4       
           k3=k3/0.4;
        else
           k3=(1-k3)/(1-0.4);
         end  
        x1(i,1)=xmin(1)+(xmax(1)-xmin(1))*k1;
        x1(i,2)=xmin(2)+(xmax(2)-xmin(2))*k2;
        x1(i,3)=xmin(3)+(xmax(3)-xmin(3))*k3;
        kp=x1(i,1); %����ϵ����ֵ
        ki=x1(i,2); %����ϵ����ֵ
        kd=x1(i,3); %΢��ϵ��
        [t,w1,u]=sim('shuixiang_pid',[0,10]);
%         x1(i,4) =max(u(:,1));
        if min(u(:,2))<0
          x1(i,4)=max(u(:,1))+c3*abs(min(u(:,2)));
        else
          x1(i,4) =max(u(:,1));
       end
    end
      x1=sortrows(x1,4); 
       for i=1:m
           for j=1:dim
               x(i,j)=x1(i,j);         
           end
           pfit(i)=x1(i,4);
       end     
    pb=x;%�������λ��
    for i=1:m   %�ٶȳ�ʼ��
        for j=1:dim
             v(i,j)=rand()*(vmax(j)-vmin(j))+vmin(j);
        end
    end
    [gbest,index]=min(pfit);%gbestȫ����С��Ӧ��
    gbestbefore=gbest;
    pg=x(index,:);%ȫ�����λ��
    %%%%%%%%%%%%%%%%%����PSO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=rand();
    for t=1:tmax
        if k<=0.4       
           k=k/0.4;
        else
           k=(1-k)/(1-0.4);
        end     
        w=wmax-(wmax-wmin)*t/tmax+0.3*(0.5-k);

        for i=1:m
            for j=1:dim        
                R1=c1*rand();
                R2=c2*rand();
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
%            pfitness=max(u(:,1));
          if min(u(:,2))<0
              pfitness=max(u(:,1))+c3*abs(min(u(:,2)));
           else
              pfitness=max(u(:,1));
           end
           if pfitness<pfit(i)
               pfit(i)=pfitness; %�������������Ӧ��
               pb(i,:)=x(i,:);%�����������λ��
           end
          if pfitness<gbest
               gbest=pfitness; %��ȫ�������Ӧ��
               pg=x(i,:);      %����ȫ�����λ��
          end      
        end 
        if abs((gbest-gbestbefore)<0.0001) &&(t>1 )
            flag=flag+1;
            if flag>10
                k1=(pg(1)-xmin(1))/(xmax(1)-xmin(1));
                k2=(pg(2)-xmin(2))/(xmax(2)-xmin(2));
                k3=(pg(3)-xmin(3))/(xmax(3)-xmin(3));
                if k1<=0.4       
                   k1=k1/0.4;
                else
                   k1=(1-k1)/(1-0.4);
                end  
                if k2<=0.4       
                   k2=k2/0.4;
                else
                   k2=(1-k2)/(1-0.4);
                end  
                 if k3<=0.4       
                   k3=k3/0.4;
                else
                   k3=(1-k3)/(1-0.4);
                 end  
                for num=1:100       %��������
                        pg1(1)=pg(1)+0.5*(0.5-k1)*pg(1);
                        pg1(2)=pg(2)+0.5*(0.5-k2)*pg(2);
                        pg1(3)=pg(3)+0.5*(0.5-k3)*pg(3);                        
                        if pg1(1)>xmax(1)
                            pg1(1)=xmax(1);
                        end
                        if pg1(2)>xmax(2)
                            pg1(2)=xmax(2);
                        end    
                         if pg1(3)>xmax(3)
                            pg1(3)=xmax(3);
                         end                               
                        kp=pg1(1);ki=pg1(2);kd=pg1(3);                       
                        [tt,w1,u]=sim('shuixiang_pid',[0,10]);      
%                         pfitness=max(u(:,1));
                       if min(u(:,2))<0
                          pfitness=max(u(:,1))+c3*abs(min(u(:,2)));
                       else
                          pfitness=max(u(:,1));
                       end
                        if pfitness<gbest
                           gbest=pfitness; 
                           pg=pg1;
                           flag=0;
                        end                       
                end
                num=unidrnd(m);%������������
                pfit(num)=gbest;
                pb(num,:)=pg;
            end
        else
            flag=0;
        end
        gbestbefore=gbest;
        gbest1(t)=gbest1(t)+gbest;    %ȫ������λ����Ӧ�ȸ��ټ�¼
        pgt1(1,t)=pgt1(1,t)+pg(1);      %�������ټ�¼
        pgt1(2,t)=pgt1(2,t)+pg(2);
        pgt1(3,t)=pgt1(3,t)+pg(3);
    end 
end
 gbest1=gbest1/countermax;
 pgt1=pgt1/countermax;
 kp1=pgt1(1,tmax);
 ki1=pgt1(2,tmax);
 kd1=pgt1(3,tmax);%���Ų������

%%%%%%%%%%%%%%%%%%%%%%%%%%����Ȩ�����Եݼ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%����Ⱥ��ʼ��
for counter=1:countermax
    for i=1:m    %λ�ó�ʼ��
        for j=1:dim
            x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j);       
        end
    end
    pb=x;%�������λ��
    for i=1:m   %�ٶȳ�ʼ��
        for j=1:dim
             v(i,j)=rand()*(vmax(j)-vmin(j))+vmin(j);
        end
    end

    for i=1:m    %�������ӳ�ʼλ����Ӧ��
        kp=x(i,1); %����ϵ����ֵ
        ki=x(i,2); %����ϵ����ֵ
        kd=x(i,3); %΢��ϵ����ֵ
        [t,w1,u]=sim('shuixiang_pid',[0,10]);
       if min(u(:,2))<0
          pfit(i)=max(u(:,1))+c3*abs(min(u(:,2)));
       else
          pfit(i) =max(u(:,1));
       end
    end
    [gbest,index]=min(pfit);%gbestȫ����С��Ӧ��
    pg=x(index,:);%ȫ�����λ��

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
               pfit(i)=pfitness; %�������������Ӧ��
               pb(i,:)=x(i,:);%�����������λ��
           end
          if pfitness<gbest
               gbest=pfitness; %��ȫ�������Ӧ��
               pg=x(i,:);      %����ȫ�����λ��
          end      
        end 
        gbest2(t)=gbest2(t)+gbest;    %ȫ������λ����Ӧ�ȸ��ټ�¼
        pgt2(1,t)=pgt2(1,t)+pg(1);      %�������ټ�¼
        pgt2(2,t)=pgt2(2,t)+pg(2);
        pgt2(3,t)=pgt2(3,t)+pg(3);
     
    end
end
 gbest2=gbest2/countermax;
 pgt2=pgt2/countermax;
 kp2=pgt2(1,tmax);
 ki2=pgt2(2,tmax);
 kd2=pgt2(3,tmax);%���Ų������
%%%%%%%%%%%%%%%%%%%%%%%�����Ż�����ͳһ����Ƚ�

sim('PSO_PID',[0,10]);  
chaotiao1=max(PID_out.signals.values(:,1));
chaotiao1=100*(chaotiao1-1)/1;
chaotiao2=max(PID_out.signals.values(:,2));
chaotiao2=100*(chaotiao2-1)/1;
chaotiao3=max(PID_out.signals.values(:,3));
chaotiao3=100*(chaotiao3-1)/1;
figure(1);
t=1:tmax;
plot(t,gbest1(t),'-o',t,gbest2(t),'-x');
legend('����PSO','��׼PSO');