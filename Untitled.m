 clc;
close all;
clear all;
x = [7.1 	7.2 	7.3   7.4 	7.5 	7.6 	7.7  7.80 	7.90 	8.00 	8.10 	8.20 	8.30 	8.40  8.5 	8.6 	8.7 	8.8 	8.9 	9.0 	9.1 	9.2 	9.3 	9.4 	9.5 	9.6 	9.7 	9.8  9.9 	10.0 	10.1 	10.2 	10.3 	10.4 	10.5 	10.6 	10.7 	10.8 	10.9 	11.0 	11.1 	11.2  11.3 	11.4 	11.5 	11.6 	11.7 	11.8 
];
y = [ -0.07 	-0.12 	-0.19 	-0.27 	-0.35 	-0.44 	-0.52 	-0.63 	-0.75 	-0.86 	-0.99 	-1.11 	-1.26 	-1.42  -1.58 	-1.76 	-1.93 	-2.13 	-2.37 	-2.60 	-2.85 	-3.09 	-3.39 	-3.73 	-4.06 	-4.43 	-4.77 	-5.20  -5.69 	-6.13 	-6.66 	-7.14 	-7.73 	-8.42 	-9.00 	-9.65 	-10.35 	-11.05 	-11.62 	-11.94 	-12.07 	-12.12  -12.16 	-12.19 	-12.22 	-12.23 	-12.25 	-12.26 
];



plot(x,y,'r-*');
xlabel('X(mm)','FontSize', 16) 
ylabel('V(v)','FontSize', 16)
text(x,y,num2str([x;y].','(%.2f,%.2f)'))
title('��');
hold on;

A=49.20;
B=-5.55;
Y=A+B*x;
plot(x,Y,'b-o');
legend('V-X����','���ֱ��');
xlabel('X(mm)','FontSize', 16) 
ylabel('V(v)','FontSize', 16)
%text(x,Y,num2str([x;Y].','(%.2f,%.2f)'))
title('��');