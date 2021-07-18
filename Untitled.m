%%% ���������ܡ�180ҳϰ��1
%%%ʯ��ƽ
%%%2021.5.7
clc;
clear all;
close all;
dim=10;%����ά��
tmax=200;%����������
m=30;%
for j=1:dim
    xmax(j)=5.12;
    xmin(j)=-5.12;
    vmax(j)=0.2*xmax(j);
    vmin(j)=0.2*xmin(j);
end
c1=2; 
c2=2; %ѧϰ����
wmax=0.9; %����Ȩ�����ֵ wpso
wmin=0.4; %����Ȩ����Сֵ wpso
  gbest_bbpso=zeros(1,tmax);
  gbest_pso=zeros(1,tmax);
nummax=1; %���Դ���
best_bbpso=zeros(1,nummax); %ÿ�����������ֵ
best_pso=zeros(1,nummax); 

%%%%%%%%%%%%%%%%%%%BBPSO
for number=1:nummax   
    for i=1:m    %      
        for j=1:dim  
            x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j); 
        end   
        pfit(i)=0;    
        for j=1:dim
            pfit(i)=pfit(i)+(x(i,j)^2-10*cos(2*pi*x(i,j))+10);       
        end     
    end
    pb=x;%���ӳ�ʼ���λ�ø�ֵ(�������λ��)
   [gbest,index]=min(pfit);%gbestȫ����С��Ӧ��
   pg=x(index,:);%ȫ�����λ��
    for t=1:tmax    
        for i=1:m   
            for j=1:dim
                u=0.5*(pb(i,j)+pg(j));
                v=abs(pb(i,j)-pg(j));
                x(i,j)=normrnd(u,v);%��˹��̬�ֲ�
               if x(i,j)>xmax(j)||x(i,j)<xmin(j)
                    x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j); %������
                end                         
            end   
            pfitness=0;
            for j=1:dim
             pfitness=pfitness+(x(i,j)^2-10*cos(2*pi*x(i,j))+10);
            end 
            if pfitness<pfit(i)
                pfit(i)=pfitness;
                pb(i,:)=x(i,:);             
           end
           if pfitness<gbest
                gbest=pfitness;
                pg=x(i,:);
                index=i;
           end      
        end
      gbest_bbpso(t)=gbest_bbpso(t)+gbest;
    end
    best_bbpso(number)=gbest;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ȩ�����Եݼ���ʽ
for number=1:nummax
   for i=1:m
         pfit(i)=0;
         for j=1:dim 
             x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j); 
             pfit(i)=pfit(i)+(x(i,j)^2-10*cos(2*pi*x(i,j))+10);
         end         
   end  
   pb=x;
   for i=1:m
       for j=1:dim
           v(i,j)=rand()*(vmax(j)-vmin(j))+vmin(j);
       end
   end
   [gbest,index]=min(pfit);%gbestȫ����С��Ӧ��
    pg=x(index,:);%ȫ�����λ��
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
            pfitness=0;
            for j=1:dim
                pfitness=pfitness+(x(i,j)^2-10*cos(2*pi*x(i,j))+10);
            end   
            if pfitness<pfit(i)
                pfit(i)=pfitness;
                pb(i,:)=x(i,:);
           end
            if pfitness<gbest
                gbest=pfitness;
                pg=x(i,:);
            end      
        end         
        gbest_pso(t)=gbest_pso(t)+gbest;
    end    
     best_pso(number)=gbest;
end
figure(1);
   i=1:tmax;
   %plot(i,(gbest_bbpso(i)/nummax),'-r',i,(gbest_pso(i)/nummax),'--r');
 % plot(i,(gbest_bbpso(i)/nummax),'-r');
 plot(i,(gbest_pso(i)/nummax),'--r');
  %legend('BBPSO','PSO');
 % legend('BBPSO');
  legend('PSO');
xlabel('��������');
ylabel('������Ӧֵ');         

