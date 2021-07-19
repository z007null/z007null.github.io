%������Ȩ�ص�PSO-PID���������ô���Դ�ڽ�ѧ�����޸Ķ���
clc;
clear;
close all;
m=60;%������
%dim=3;%����ά��
dim=3;%����ά��
vmax=[2,1,1];%�ٶ����ֵ kp,ki,kd
vmin=[-2,-1,-1];%�ٶ���Сֵkp,ki,kd
xmax=[10,5,5];%λ�����ֵkp,ki,kd
xmin=[0,0,0];%λ����Сֵ
tmax=100;%����������
c1=2; 
c2=2; %ѧϰ����
c3=2;
x=rand(m,dim);
v=rand(m,dim);
pfit=rand(1,m);
wmax=0.9; %����Ȩ�����ֵ wpso
wmin=0.4; %����Ȩ����Сֵ wpso
countermax=1;
gbest2=zeros(1,tmax);
pgt2=zeros(3,tmax);

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
 %%%%%%%%%%%%%%%%%%%%%%%�����Ż�����ͳһ����
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
legend('��׼PSO');
xlabel('��������t');
ylabel('outputĿ�꺯��ֵ');