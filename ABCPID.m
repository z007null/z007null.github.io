% �˹���Ⱥ�㷨
% ����˵����
% Foods [FoodNumber][D]; % ��ʼ����ʳ��Դ
% ObjVal[FoodNumber];    % Ŀ�꺯��
% Fitness[FoodNumber];   % ��Ӧ��ֵ��Ŀ�꺯��ֵ�ĵ���
% trial[FoodNumber];     % ��β����
% prob[FoodNumber];      % ����ĸ���ֵ
% solution [D];          % �������½⣬��ѡλ�� produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
% ObjValSol;             % �½��µ�Ŀ�꺯��ֵ
% FitnessSol;            % �½����Ӧ��ֵ
% neighbour, param2change; ��Ӧ�ڷ��� v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
% GlobalMin;             % Ŀ�꺯��ֵ��Сֵ
% GlobalParams[D];       % ÿһ�����и��㷨�õ������Ÿ���ֵ��δ֪���Ľ�
% GlobalMins[runtime];   % ѭ��������㷨�Ĵ�������¼�µ���С�⣬��֤�㷨��³���Ժ��ȶ���
clc,clear,close all
warning off
feature jit off
tic
% �㷨����
NP=20;           % ��Ⱥ��С
dim=3;%�۷�ά��
c3=2;
FoodNumber=NP/2; % ��Ⱥʳ��Դ������Ҳ���ǲ��� �� �ĸ���
limit=10;       % ������limit���β��۷�͹۲���ѭ������֮�󣬲��ܹ����Ľ�����ô��λ�ý�������
maxCycle=500;    % ������ѭ��
x=rand(NP,dim);

%/* Problem specific variables*/
for i=1:NP    %�������ӳ�ʼλ����Ӧ��
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
%objfun='Sphere';     % ���Ż�����
%objfun='pfit';
%D=10;               % δ֪��Ϊ100��
Xmax=[2,1,1];%�ٶ����ֵ kp,ki,kd
Xmin=[-2,-1,-1];%�ٶ���Сֵkp,ki,kd
D=Xmax-Xmin;
%ub=ones(1,D)*10;    % δ֪��ȡֵ�±߽�
%lb=ones(1,D)*(-10); % δ֪��ȡֵ�ϱ߽�
runtime=1;           % �㷨���д�����һ������1����

GlobalMins=zeros(1,runtime);   % ��Ӧ����Сֵ��ʼ��

for r=1:runtime
  
%��ʼ������ֵ
% Range = repmat((ub-lb),[FoodNumber 1]);       % ���ֵ
% Lower = repmat(lb, [FoodNumber 1]);           % ��Сֵ
% Foods = rand(FoodNumber,D) .* Range + Lower;  % ��ʼ������
 for i=1:NP    %λ�ó�ʼ��
        for j=1:dim
            Foods(i,j)=rand()*(Xmax(j)-Xmin(j))+Xmin(j);
           %x(i,j)=rand()*(xmax(j)-xmin(j))+xmin(j);       
        end
    end
%ObjVal=feval(objfun,Foods);       % Ŀ�꺯��ֵ
ObjVal=feval(pfit,Foods);       % Ŀ�꺯��ֵ
Fitness=calculateFitness(ObjVal); % ��Ӧ��ֵ��ȡ�䵼����Ϊ��Сֵ

% �趨��β���󣬳�ʼ��
trial=zeros(1,FoodNumber);

% �ҵ���õ�ʳ��Դ
BestInd=find(ObjVal==min(ObjVal));
BestInd=BestInd(end);
GlobalMin=ObjVal(BestInd);     % ����ֵ��С
GlobalParams=Foods(BestInd,:); % ��Ӧ��ʳ��Դ����

iter=1;
while ((iter <= maxCycle)) % ������ʼ

% ���۷�
    for i=1:(FoodNumber)
        % ��������ɱ�
        Param2Change=fix(rand*D)+1;
        % ���ѡ����������
        neighbour=fix(rand*(FoodNumber))+1;
        % ���ѡ��ĸ��岻����i
        while(neighbour==i)
            neighbour=fix(rand*(FoodNumber))+1;
        end
        
       sol=Foods(i,:);  % ����ѡ��
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       % ����ȡֵ��ΧԼ��
        ind=find(sol<Xmax); % ��СֵԼ��
        sol(ind)=Xmax(ind);
        ind=find(sol>Xmin); % ���ֵԼ��
        sol(ind)=Xmin(ind);
        
        % �����µ�Ŀ�꺯��ֵ����Ӧ��ֵ
        ObjValSol=feval(pfit,sol);
        FitnessSol=calculateFitness(ObjValSol);
        
       % �������Ÿ���ֵ
       if (FitnessSol>Fitness(i)) % ����²����ĸ���ֵ��Ӧ��ֵԽ�����������ֵԽС�����������
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; % /*if the solution i can not be improved, increase its trial counter*/
       end
    end
    
% �۲��
% �������
% �۲���������Դ��صĸ���ֵѡ����Դ������ֵ���㹫ʽ
% prob(i)=a*fitness(i)/max(fitness)+b*/
prob=(0.9.*Fitness./max(Fitness))+0.1;  
i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        % �������ѡ�����
        Param2Change=fix(rand*D)+1;
        % ���ѡ����������
        neighbour=fix(rand*(FoodNumber))+1;
        % ���ѡ��ĸ��岻����i      
        while(neighbour==i)
            neighbour=fix(rand*(FoodNumber))+1;
        end
       sol=Foods(i,:);  % ����ѡ��
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       % ����ȡֵ��ΧԼ��
       ind=find(sol<Xmax); % ��СֵԼ��
       sol(ind)=Xmax(ind);
       ind=find(sol>Xmin); % ���ֵԼ��
       sol(ind)=Xmin(ind);
        
       % �����µ�Ŀ�꺯��ֵ����Ӧ��ֵ
       ObjValSol=feval(pfit,sol);
       FitnessSol=calculateFitness(ObjValSol);
        
       % �������Ÿ���ֵ
       if (FitnessSol>Fitness(i)) % ����²����ĸ���ֵ��Ӧ��ֵԽ�����������ֵԽС�����������
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; % /*if the solution i can not be improved, increase its trial counter*/
       end
    end
    
    i=i+1;
    if (i==(FoodNumber)+1) 
        i=1;
    end   
end 

    % ��¼��õ�Ŀ�꺯��ֵ
    ind=find(ObjVal==min(ObjVal));
    ind=ind(end);
    if (ObjVal(ind)<GlobalMin)
        GlobalMin=ObjVal(ind);      % ����Ŀ�꺯��ֵ
        GlobalParams=Foods(ind,:);  % ���Ÿ���
    end
         
         
% ����
% ���ĳһ��ѭ����β���������趨limit�������¸��¸��壬���¼���
ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    Bas(ind)=0;
    sol=(Xmin-Xmax).*rand(1,D)+Xmax;
    ObjValSol=feval(pfit,sol);
    FitnessSol=calculateFitness(ObjValSol);
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
end

fprintf('iter=%d ObjVal=%g\n',iter,GlobalMin);
iter=iter+1;

end % End of ABC

GlobalMins(r)=GlobalMin;
end % end of runs
disp('���Ž�Ϊ��')
GlobalParams;
disp('����Ŀ�꺯��ֵΪ��')
GlobalMin;

toc
% save all