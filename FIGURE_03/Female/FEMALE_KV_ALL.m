set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 20)

data1 = load('female_all.txt');
x1 = data1(:,1);
x2 = data1(:,2);
x3 = data1(:,3);
x4 = data1(:,4);

%Data 
y1 = [	-70	-60	-50	-40	-30	-20	-10	0	10	20	30	40	50 60];
%IV 
y2 =[0	0	0	8	14.6774	24.8173	68.4717	123.5526	195.8348	260.1197	337.1074	406.3862	492.2415	568.1691];

%Error
y3 =[0 0 0 0 4.254753266 5.604139861 8.767539316 9.587487346 13.58646111 17.13589524 23.14388554 25.55184375 32.86695764 40.81273552];


data2 = load('kv15_female.txt');
cx1 = data2(:,1);
cx2 = data2(:,2);
cx3 = data2(:,3);
cx4 = data2(:,4);

%Data 
f1 = [	-70	-60	-50	-40	-30	-20	-10	0	10	20	30	40	50 60];
%IV 
f2 = [0	0	0	0	1.677403286	6.395005191	28.01841109	56.68694544	81.59176891	113.7719204	149.976736	179.959441	225.3015465	275.9825858];
%Error
f3 =[0 0 0 0 4.254753266 5.604139861 8.767539316 9.587487346 13.58646111 17.13589524 23.14388554 25.55184375 32.86695764 40.81273552];
f4 = [0	0	0	0	0.051586948	0.087056628	0.315658645	0.544723762	0.683524246	0.844801341	1	1	1	1];
f5 = [0	0	0	0	0.015955325	0.021015524	0.032878272	0.035953078	0.050949229	0.064259607	0.086789571	0.095819414	0.123251091	0.153047758];



data3 = load('kv21_female.txt');
vx1 = data3(:,1);
vx2 = data3(:,2);
vx3 = data3(:,3);
vx4 = data3(:,4);

%Data 
vy1 = [	-70	-60	-50	-40	-30	-20	-10	0	10	20	30	40	50 60];
%IV not normalized
vy2 =[0 0 0 0 10	18.42231578	40.45334159	66.8656923 114.2430168 146.3478395 187.1307253 226.4268017 266.939973 292.1865241];
%Error y2
vy3 = [4.613731845	5.585343992	5.329290873	5.752428129	9.317391009	5.756191789	9.368621419	10.14399681	14.95189252	18.78395898	24.10726415	24.60194821	28.90006495	29.58321458];


%Activation
vy4 =[0	0	0	0.05	0.137815282	0.200994401	0.365264431	0.514962266	0.767036508	0.870932742	1	1	1	1];
vy5 =[0	0	0	0	0.03180688	0.037480334	0.0349146	0.044982	0.077329227	0.037771624	0.05847791	0.0847791	0.038477914	0.13847791];


%Tau activation 
vy6 =[	-30	-20	-10	0	10	20	30	40	50 60];
vy7 = [172	144	116	100	84	68	60	48	40	32];
vy8 = [25	21	25	26	15	14.6	10.5	9.5	7.5	5];

figure (100)
set(gcf,'color','w');
errorbar(y1,y2,y3,'oc', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(x1,x2,'c', 'LineWidth',2)
axis([-60 40 0 700])
ylabel('Female I_{Kv_{ALL}} (pA)'),xlabel('Voltage (mV)')
legend('Data N = 10','Model')


figure (200)
set(gcf,'color','w');
errorbar(f1,f2,f3,'sc', 'LineWidth',1.2,'MarkerSize',10')
hold on
plot(cx1,cx2,'c', 'LineWidth',2)
axis([-60 40 0 350])
ylabel('Female I_{Kv1.5} (pA)'),xlabel('Voltage (mV)')
legend('Data N = 10','Model')


figure (300)
set(gcf,'color','w');
errorbar(f1,f4,f5,'sc', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(cx1,cx3,'c', 'LineWidth',2)
axis([-60 40 0 1])
ylabel('Female Kv1.5 G/G{max}'),xlabel('Voltage (mV)')
legend('Data N = 10','Model')

figure (400)
set(gcf,'color','w');
plot(cx1,cx4,'c', 'LineWidth',2)
axis([-60 50 0 50])
ylabel('Female Kv1.5 \tau_{activation}'),xlabel('Voltage (mV)')
legend('Data N = 10','Model')

figure (500)
set(gcf,'color','w');
errorbar(vy1,vy2,vy3,'sc', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(vx1,vx2,'c', 'LineWidth',2)
axis([-60 40 0 350])
ylabel('Female I_{Kv2.1} (pA)'),xlabel('Voltage (mV)')
legend('Data N = 20','Model')

figure (600)
set(gcf,'color','w');
errorbar(vy1,vy4,vy5,'sc', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(vx1,vx3,'c', 'LineWidth',2)
axis([-60 40 0 1])
ylabel('Female Kv2.1 G/Gmax'),xlabel('Voltage (mV)')
legend('Data N = 10','Model')

figure (700)
set(gcf,'color','w');
errorbar(vy6,vy7,vy8,'sc', 'LineWidth',1.2,'MarkerSize',10)
hold on
plot(vx1,vx4,'c', 'LineWidth',2)
axis([-60 40 0 200])
ylabel('Female Kv2.1 \tau_{Activation (ms)}'),xlabel('Voltage (mV)')
legend('Data N = 10','Model')


