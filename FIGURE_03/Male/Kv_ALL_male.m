set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 20)

data1 = load('KV_ALL_male.txt');
x1 = data1(:,1);
x2 = data1(:,2);
x3 = data1(:,3);
x4 = data1(:,4);
%Data 
y1 = [	-70	-60	-50	-40	-30	-20	-10	0	10	20	30	40	50 60];
%IV not normalized
y2 = [0	0	0	0	14.02788583	30.22651603	56.39344292	90.79726631	136.6200122	193.0895382	247.8605138	303.3609704	362.3538233	431.356572];
%Erros
y3 =[0	0	0	0	3.60383065	4.549188659	6.897883749	12.59313828	14.41411049	20.20465969	24.95499618	32.52994281	38.01065582	37.28359852];



data2 = load('kv15_male.txt');
x11 = data2(:,1);
x22 = data2(:,2);
x33 = data2(:,3);
x44 = data2(:,4);

%Data 
y11 = [	-70	-60	-50	-40	-30	-20	-10	0	10	20	30	40	50 60];
%IV not normalized
y22 = [0	0 0	6 12.02788583 24.10676054 44.39344292 67.79726631 99.62001223 138.8750147 182.5394671 234.5953481 287.6804212 354.030152];
%Erros
y33 =[3.268124431 2.928658665 3.08201154	2.371395139	3.60383065 4.549188659	6.897883749	12.59313828	14.41411049	20.20465969	24.95499618	32.52994281	38.01065582	37.28359852];

y44 =[0	0	0	0.0651	0.169931933	0.269629194	0.410922699	0.535269586	0.685679689	0.847248373	1	1	1	1];
y55 =[0	0	0	0	0.015955325	0.021015524	0.032878272	0.035953078	0.050949229	0.064259607	0.086789571	0.095819414	0.123251091	0.153047758];


data3 = load('kv21_male.txt');
f1 = data3(:,1);
f2 = data3(:,2);
f3 = data3(:,3);
f4 = data3(:,4);

%Data 
fy1 = [-70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60];
%IV not normalized
fy2 = [0	0 0	1 2	6.11975549 12 23 37	54.21452344	65.32104673	68.76562236	74.67340203	77.32641993];
%Erros
fy3 =[3.297159978	2.046862447	1.863660244	2.660742991	4.786143035	3.051819829	9.170201639	9.681767434	13.04046419	14.9026516	23.48465848	24.30543041	25.3348589	35.2928053];


fy4=[0	0	0	0	0.078962218	0.191278331	0.310403202	0.507448372	0.711672299	0.924285016	1	1	1	1];
%Activation
fy5 =[0	0	0	0	0.023818069	0.020748033	0.050044915	0.081049498	0.107732923	0.143777162	0.178384779	0.202678356	0.229117854	0.300970935];


%Tau activation
fy6 =[	-30	-20	-10	0	10	20	30	40	50 60];
fy7 =[43	36	29	25	21	17	15	12	10	8];
fy8 = [5	4.2	3	5	6	3.2	2.1	2.9	2	1];





figure (1)
set(gcf,'color','w');
errorbar(y1,y2,y3,'ok', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(x1,x2,'k', 'LineWidth',2)
axis([-60 40 0 700])
ylabel('Male I_{Kv_{TOT}} (pA)'),xlabel('Voltage (mV)')
legend('Data N = 7','Model')

figure (2)
set(gcf,'color','w');
errorbar(y11,y22,y33,'ok', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(x11,x22,'k', 'LineWidth',2)
axis([-60 40 0 350])
ylabel('Male I_{Kv15} (pA)'),xlabel('Voltage (mV)')
legend('Data N = 7','Model')

figure (3)
set(gcf,'color','w');
errorbar(y11,y44,y55,'ok', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(x11,x33,'k', 'LineWidth',2)
axis([-60 40 0 1])
ylabel('Male Kv1.5 G/G{max}'),xlabel('Voltage (mV)')
legend('Data N = 6','Model')


figure (4)
set(gcf,'color','w');
plot(x11,x44,'k', 'LineWidth',2)
axis([-60 40 0 50])
ylabel('Male Kv1.5 \tau_{activation}'),xlabel('Voltage (mV)')
legend('Model')

figure (5)
set(gcf,'color','w');
errorbar(fy1,fy2,fy3,'ok', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(f1,f2,'k', 'LineWidth',2)
axis([-60 40 0 350])
ylabel('Male I_{Kv2.1} (pA)'),xlabel('Voltage (mV)')
legend('Data N = 10','Model')

figure (6)
set(gcf,'color','w');
errorbar(fy1,fy4,fy5,'ok', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(f1,f3,'k', 'LineWidth',2)
axis([-60 40 0 1])
ylabel('Male Kv2.1 G/Gmax'),xlabel('Voltage (mV)')
legend('Data N = 7','Model')

figure (7)
set(gcf,'color','w');
errorbar(fy6,fy7,fy8,'ok', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(f1,f4,'k', 'LineWidth',2)
axis([-60 40 0 80])
ylabel('Male Kv2.1 \tau_{Activation (ms)}'),xlabel('Voltage (mV)')
legend('Data N = 7','Model')

