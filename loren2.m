%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���øĽ�������Ⱥ�㷨���LORENZ����ϵͳ�Ĳ������ƣ���������δ֪��
%k1,k2,k3 ���ϵ�� 
%���ڣ�2016.6.20 
%�޶���2018.5.9
%�ص㣺����С��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clc;
clc;
clear;
close all;
m=30;%������
dim=3;%����ά��
xmax=[11,30,3];%λ�����ֵ
xmin=[9,20,2];%λ����Сֵ
vmax=[0.2,1,0.1];%�ٶ����ֵ 
vmin=[-0.2,-1,-0.1];%�ٶ���Сֵ
tmax=30;%����������

nummax=1;
c1=2; 
c2=2; %ѧϰ����

wmax=0.9; %����Ȩ�����ֵ wpso
wmin=0.4; %����Ȩ����Сֵ wpso


aa2=zeros(1,tmax);
bb2=zeros(1,tmax);
cc2=zeros(1,tmax);
g_fitness2=zeros(1,tmax);
pgnum2=zeros(nummax,3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSO
for num=1:nummax
    %%%λ�ó�ʼ��,����������Ӧ��
        for i=1:m
            for j=1:dim    
                    x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j); 
                   v(i,j)=rand()*(vmax(j)-vmin(j))+vmin(j);      
            end 
             k1=x(i,1);k2=x(i,2);k3=x(i,3);
             [n,w1,u]=sim('lorenz1',0.1);%������
%              [n,w1,u]=sim('lorenz',5);%������
             fit(i)=max(u(:,1));  
        end

        [gbest,index]=min(fit);%gbestȫ����С��Ӧ��
        pg=x(index,:);%ȫ�����λ��
        pb=x;%�������λ��
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
             [n,w1,u]=sim('lorenz1',0.1);%������
%             [n,w1,u]=sim('lorenz',5);%������
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
        aa2(t)=aa2(t)+pg(1);%���Ʋ���a����
        bb2(t)=bb2(t)+pg(2);%���Ʋ���b����
        cc2(t)= cc2(t)+pg(3);%���Ʋ���c����
        g_fitness2(t)=g_fitness2(t)+gbest;%ͳ����Ӧ��ֵ��������  
    end    
    pgnum2(num,:)=pg;    %ͳ��ÿ�������Ѱ�Ž��
end


a_mean2=mean(pgnum2(:,1));%a������ƽ������ֵ
b_mean2=mean(pgnum2(:,2));%b������ƽ������ֵ
c_mean2=mean(pgnum2(:,3));%c������ƽ������ֵ


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
  