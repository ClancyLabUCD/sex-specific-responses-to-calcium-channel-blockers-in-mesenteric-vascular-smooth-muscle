set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 2)
set(groot, 'defaultAxesFontSize', 20)

data1 = load('all_ODEs_Male.txt');
Time = data1(:,1);
V = data1(:,2);
dL= data1(:,3);
dF = data1(:,4);
p_K  = data1(:,5);
p_K15 = data1(:,6);
Na_in = data1(:,7);
K_in = data1(:,8);
Cl_in  = data1(:,9);
Ca_in = data1(:,10);
Ca_SRcen = data1(:,11);
Ca_jun= data1(:,12);
BUF_1= data1(:,13);
CSQ_SRcen1 = data1(:,14);
BUF_jun1 = data1(:,15);
R_10= data1(:,16);
xab = data1(:,17);
p_cc = data1(:,18);

data2 = load('Currents_Male.txt');
Time = data2(:,1);

I_CaL = data2(:,2);
I_Kv21= data2(:,3);
I_Kv15 = data2(:,4);
I_KvALL= data2(:,5);
I_BK = data2(:,6);

I_NSC_Na  = data2(:,7);
I_NSC_K = data2(:,8);
I_NaK =data2(:,9);
I_NCX   = data2(:,10);
I_PMCA  = data2(:,11);

I_Ca_leak  = data2(:,12);
I_Na_leak  = data2(:,13);
I_K_leak = data2(:,14);
I_Cl_leak = data2(:,15);
I_ALL_leak = data2(:,16);

J_SERCA = data2(:,17);
J_rel = data2(:,18);
I_Cl = data2(:,19);







figure(30)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 20)
subplot(4,4,1), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time*1000/1000,V, '-k'), ylabel('Vm (mV)')
subplot(4,5,2), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,dL, '-k'), ylabel('dL'), xlabel('Time (ms)')
subplot(4,5,3), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,dF, '-k'), ylabel('dF '), xlabel('Time (ms)')
subplot(4,5,4), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,p_K, '-k'), ylabel('p_K '), xlabel('Time (ms)')
subplot(4,5,5), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,p_K15, '-k'), ylabel('p_{K15}'), xlabel('Time (ms)')
subplot(4,5,6), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time,Na_in, '-k'), ylabel('Na_{in} (mM)'), xlabel('Time (ms)')
subplot(4,5,7), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,K_in, '-k'), ylabel('K_{in} (mM)'), xlabel('Time (ms)')
subplot(4,5,8), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,Cl_in, '-k'), ylabel('Cl_{in} (mM)'), xlabel('Time (ms)')
subplot(4,5,9), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,Ca_in, '-k'), ylabel('Ca_{Cyt} (mM)'), xlabel('Time (ms)')
subplot(4,5,10), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,Ca_SRcen, '-k'), ylabel('Ca_{SRcen} (mM)'), xlabel('Time (ms)')
subplot(4,5,11), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,Ca_jun, '-k'), ylabel('Ca_{jun} (mM)'), xlabel('Time (ms)')
 subplot(4,5,12), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,BUF_1, '-k'), ylabel('BUF_1'), xlabel('Time (ms)')  
subplot(4,5,13), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,CSQ_SRcen1, '-k'), ylabel('CSQ_{SRper1}'), xlabel('Time (ms)')
subplot(4,5,14), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,BUF_jun1, '-k'), ylabel('BUF_{jun1}'), xlabel('Time (ms)')
subplot(4,5,15), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,R_10, '-k'), ylabel('R_{10}'), xlabel('Time (ms)')
subplot(4,5,16), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,xab, '-k'), ylabel('xab'), xlabel('Time (ms)')
 subplot(4,5,17), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,p_cc, '-k'), ylabel('pcc'), xlabel('Time (ms)')

    
figure(31)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 20)

subplot(4,5,1), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_CaL, '-k'), ylabel('I_{CaL} (pA)')
subplot(4,5,2), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Kv21, '-k'), ylabel('I_{Kv_{2.1}} (pA)'), xlabel('Time (ms)')
subplot(4,5,3), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Kv15, '-k'), ylabel('I_{Kv_{1.5}} (pA)'), xlabel('Time (ms)')
subplot(4,5,4), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_KvALL, '-k'), ylabel('I_{Kv_{ALL}} (pA)'), xlabel('Time (ms)') 
subplot(4,5,5), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_BK, '-k'), ylabel('I_{BK} (pA)'), xlabel('Time (ms)')
    
subplot(4,5,6), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000, I_NSC_Na, '-k'), ylabel('I_{NSC_{Na}} (pA)'), xlabel('Time (ms)')
subplot(4,5,7), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_NSC_K, '-k'), ylabel('I_{NSC_{K}} '), xlabel('Time (ms)')
subplot(4,5,8), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_NaK, '-k'), ylabel('I_{NaK} (pA)'), xlabel('Time (ms)')
subplot(4,5,9), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_NCX, '-k'), ylabel('I_{NCX} (pA) '), xlabel('Time (ms)')
subplot(4,5,10), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000, I_PMCA, '-k'), ylabel('I_{PMCA} (pA)'), xlabel('Time (ms)')

subplot(4,5,11), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Ca_leak, '-k'), ylabel('I_{Ca_{leak}} (pA)'), xlabel('Time (ms)')
subplot(4,5,12), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Na_leak, '-k'), ylabel('I_{Na_{leak}} (pA)'), xlabel('Time (ms)')  
subplot(4,5,13), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_K_leak, '-k'), ylabel('I_{K_{leak}} (pA)'), xlabel('Time (ms)')
subplot(4,5,14), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Cl_leak, '-k'), ylabel('I_{Cl_{leak}} (pA)'), xlabel('Time (ms)')
    
subplot(4,5,15), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_ALL_leak, '-k'), ylabel('I_{ALL_{leak}} (pA)'), xlabel('Time (ms)')
subplot(4,5,16), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,J_SERCA, '-k'), ylabel('J_{SERCA} (pA)'), xlabel('Time (ms)')   
subplot(4,5,17), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,J_rel*2*96.4853415*0.005*1000 , '-b'), ylabel('J_{rel} (pA)'), xlabel('Time (ms)')
subplot(4,5,18), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Cl , '-k'), ylabel('J_{rel} (pA)'), xlabel('Time (ms)')
         


figure (32)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(4,1,1)      
    plot(Time/1000,V, '-k'), ylabel('Vm (mV)') , xlabel('Time (ms)')         
subplot(4,1,2)       
    plot(Time/1000,I_BK, '-k'), ylabel('I_{BK} (pA)'), xlabel('Time (ms)')      
subplot(4,1,3)      
    plot(Time/1000,I_CaL, '-k'), ylabel('I_{CaL} (pA)'), xlabel('Time (ms)')            
subplot(4,1,4)       
   plot(Time/1000,Ca_in*1E6, '-k'), ylabel('Ca_{Cyt} (nM)'), xlabel('Time (ms)') 

    
    figure (33)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(4,1,1)      
    plot(Time/1000,V, '-k'), ylabel('Vm (mV)') , xlabel('Time (ms)')         
subplot(4,1,2)       
    plot(Time/1000,I_KvALL, '-k'), ylabel('I_{Kv_{ALL}} (pA)'), xlabel('Time (ms)')      
subplot(4,1,3)      
    plot(Time/1000,I_Kv21, '-k'), ylabel('I_{Kv_{2.1}} (pA)'), xlabel('Time (ms)')            
subplot(4,1,4)       
    plot(Time/1000,I_Kv15, '-k'), ylabel('I_{Kv_{1.5}} (mM)'), xlabel('Time (ms)') 
    
    
       figure (34)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(3,1,1)      
       plot(Time/1000,Ca_in*1E6, '-k'), ylabel('Ca_{in} (nM)') , xlabel('Time (s)')  
    ylim([0 200])     
subplot(3,1,2)      
    plot(Time/1000,Ca_SRcen*1000, '-k'), ylabel('Ca_{SR} (\muM)'), xlabel('Time (s)')    
    ylim([0 150])
subplot(3,1,3)      
    plot(Time/1000,Ca_jun*1000, '-k'), ylabel('Ca_{Jun} (\muM)'), xlabel('Time (s)')            
ylim([0 20])  

   figure (35)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(3,1,1)      
       plot(Time/1000,Ca_in*1E6, '-k'), ylabel('Ca_{in} (nM)') , xlabel('Time (s)')  
    ylim([0 200])     
subplot(3,1,2)      
    plot(Time/1000,Na_in, '-k'), ylabel('Na_{in} (mM)'), xlabel('Time (s)')    
    ylim([0 20])
subplot(3,1,3)      
    plot(Time/1000,K_in, '-k'), ylabel('K_{in} (mM)'), xlabel('Time (s)')            
ylim([0 200])  


   figure (36)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(4,1,1)      
       plot(Time/1000,Ca_in*1E6, '-k'), ylabel('Ca_{in} (nM)') , xlabel('Time (s)')  
    ylim([0 200])     
subplot(4,1,2)      
    plot(Time/1000,Ca_SRcen*1000, '-k'), ylabel('Ca_{SR} (\muM)'), xlabel('Time (s)')    
    ylim([0 150])
    
    subplot(4,1,3)       
    plot(Time/1000,I_BK, '-k'), ylabel('STOCs (pA)'), xlabel('Time (s)')  
    subplot(4,1,4)      
    plot(Time/1000,Ca_jun*1000, '-k'), ylabel('Ca_{Jun} (\muM)'), xlabel('Time (s)')            
ylim([0 20]) 



    figure (37)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
    plot(Time/1000,V, 'k'), ylabel('Vm (mV)') , xlabel('Time (s)')  
ylim([-60 -20] )        
  xlim([475 500])


   


    
    
   figure (39)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(3,1,1) 
plot(Time/1000,V, 'k'), ylabel('Vm (mV)') , xlabel('Time (s)')  
ylim([-60 -20] )   
   xlim([485 500] ) 
subplot(3,1,2)      
    plot(Time/1000,I_CaL, '-k'), ylabel('I_{CaL} (pA)'), xlabel('Time (ms)')   
    ylim([-2 0] ) 
       xlim([485 500] ) 
subplot(3,1,3)      
       plot(Time/1000,Ca_in*1E6, '-k'), ylabel('Ca_{in} (nM)') , xlabel('Time (s)')  
    ylim([0 200]) 
       xlim([485 500] ) 




