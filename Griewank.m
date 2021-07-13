%���øĽ�������Ⱥ�����Ȩ�����Եݼ�������Ⱥ�ֱ��Griewank����
%�����Ż�����������60��ά����10��������������500���ֱ����50��ȡƽ��ֵ��
%ʱ�䣺2017.6.17

clc;
clear;
m=60;%������
dim=10;%����ά��
vmax=[60,60,60,60,60,60,60,60,60,60];%�ٶ����ֵ 
vmin=[-60,-60,-60,-60,-60,-60,-60,-60,-60,-60];%�ٶ���Сֵ
xmax=[600,600,600,600,600,600,600,600,600,600];%λ�����ֵ
xmin=[-600,-600,-600,-600,-600,-600,-600,-600,-600,-600];%λ����Сֵ
tmax=300;%����������
kmax=500;
c1=2; 
c2=2; %ѧϰ����
x=rand(m,dim);
v=rand(m,dim);
pfit=rand(1,m);
flag=0;%�жϽ���ͣ�ͱ�־
wmax=0.9; %����Ȩ�����ֵ wpso
wmin=0.4; %����Ȩ����Сֵ wpso
gbest1= zeros(1,tmax);
gbest2= zeros(1,tmax);
position1=zeros(tmax,dim);%����λ�ø���
position2=zeros(tmax,dim);%����λ�ø���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin
nummax=1; %���Դ���
testbest1=zeros(1,nummax); %ÿ�����������ֵ
testbest2=zeros(1,nummax); %ÿ�����������ֵ
%%%%%%%%%%%%%%%%%%%%��������Ⱥ
for number=1:nummax      
    for i=1:m    %λ�û����ʼ��  m1��������ѡ�����ŵ�m������Ϊ��ʼ����
        result1=0;
        result2=1;
        for j=1:dim    
            x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j);
            result1=result1+(x(i,j)^2)/4000;
            result2=result2*cos(x(i,j)/j^0.5);              
        end     
        pfit(i)=1+result1-result2; 
    end    
           
    pb=x;%���ӳ�ʼ���λ�ø�ֵ
    for i=1:m   %�ٶȳ�ʼ��
        for j=1:dim
             v(i,j)=rand()*(vmax(j)-vmin(j))+vmin(j);
        end
    end 
   [gbest,index]=min(pfit);%gbestȫ����С��Ӧ��
   gbestbefore=0;%�ϴ�ȫ����Сֵ
   pg=x(index,:);%ȫ�����λ��
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ȩ�����Եݼ���ʽ
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
testbestmax1=max(testbest1);%���ֵ
testbestmin1=min(testbest1);%��Сֵ
testbestmean1=mean(testbest1);%ƽ��ֵ
%%%%%%%%%%%%%%%�����׼��
fangcha1=0;
for number=1:nummax
    fangcha1=fangcha1+(testbest1(number)-testbestmean1)^2;
end
fangcha1=sqrt(fangcha1/nummax);

testbestmax2=max(testbest2);%���ֵ
testbestmin2=min(testbest2);%��Сֵ
testbestmean2=mean(testbest2);%ƽ��ֵ
%%%%%%%%%%%%%%%�����׼��
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
legend('�Ľ�PSO','��׼PSO');
xlabel('��������');
ylabel('��Ӧֵ�Ķ���');
figure(2);
plot(1:nummax,testbest1(1:nummax),'-d',1:nummax,testbest2(1:nummax),'-o');
legend('�Ľ�PSO','��׼PSO');
xlabel('ʵ�����');
ylabel('��Ӧֵ');
% figure(4);
% plot(t,position1(t,:));
% figure(5);
% plot(t,position2(t,:));
% grid on;






 

