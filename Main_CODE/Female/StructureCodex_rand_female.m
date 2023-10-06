set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 2)
set(groot, 'defaultAxesFontSize', 20)

data1 = load('all_ODEs_female.txt');
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


data2 = load('Currents_female.txt');
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




figure(300)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 20)
subplot(4,4,1), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,V, '-b'), ylabel('Vm (mV)')
subplot(4,4,2), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,dL, '-b'), ylabel('dL'), xlabel('Time (ms)')
subplot(4,4,3), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,dF, '-b'), ylabel('dF '), xlabel('Time (ms)')
subplot(4,4,4), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,p_K, '-b'), ylabel('p_K '), xlabel('Time (ms)')
subplot(4,4,5), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,p_K15, '-b'), ylabel('p_{K15}'), xlabel('Time (ms)')
subplot(4,4,6), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,Na_in, '-b'), ylabel('Na_{in} (mM)'), xlabel('Time (ms)')
subplot(4,4,7), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,K_in, '-b'), ylabel('K_{in} (mM)'), xlabel('Time (ms)')
subplot(4,4,8), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,Cl_in, '-b'), ylabel('Cl_{in} (mM)'), xlabel('Time (ms)')
subplot(4,4,9), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,Ca_in, '-b'), ylabel('Ca_{Cyt} (mM)'), xlabel('Time (ms)')
subplot(4,4,10), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,Ca_SRcen, '-b'), ylabel('Ca_{SRcen} (mM)'), xlabel('Time (ms)')
subplot(4,4,11), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,Ca_jun, '-b'), ylabel('Ca_{jun} (mM)'), xlabel('Time (ms)')
 subplot(4,4,12), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,BUF_1, '-b'), ylabel('BUF_1'), xlabel('Time (ms)')  
subplot(4,4,13), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,CSQ_SRcen1, '-b'), ylabel('CSQ_{SRper1}'), xlabel('Time (ms)')
subplot(4,4,14), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,BUF_jun1, '-b'), ylabel('BUF_{jun1}'), xlabel('Time (ms)')
subplot(4,4,15), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,R_10, '-b'), ylabel('R_{10}'), xlabel('Time (ms)')
subplot(4,4,16), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,xab, '-b'), ylabel('xab'), xlabel('Time (ms)')

    
figure(310)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 20)

subplot(4,5,1), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_CaL, '-b'), ylabel('I_{CaL} (pA)')
subplot(4,5,2), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Kv21, '-b'), ylabel('I_{Kv_{2.1}} (pA)'), xlabel('Time (ms)')
subplot(4,5,3), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Kv15, '-b'), ylabel('I_{Kv_{1.5}} (pA)'), xlabel('Time (ms)')
subplot(4,5,4), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_KvALL, '-b'), ylabel('I_{Kv_{ALL}} (pA)'), xlabel('Time (ms)') 
subplot(4,5,5), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_BK, '-b'), ylabel('I_{BK} (pA)'), xlabel('Time (ms)')
    
subplot(4,5,6), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000, I_NSC_Na, '-b'), ylabel('I_{NSC_{Na}} (pA)'), xlabel('Time (ms)')
subplot(4,5,7), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_NSC_K, '-b'), ylabel('I_{NSC_{K}} '), xlabel('Time (ms)')
subplot(4,5,8), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_NaK, '-b'), ylabel('I_{NaK} (pA)'), xlabel('Time (ms)')
subplot(4,5,9), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_NCX, '-b'), ylabel('I_{NCX} (pA) '), xlabel('Time (ms)')
subplot(4,5,10), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000, I_PMCA, '-b'), ylabel('I_{PMCA} (pA)'), xlabel('Time (ms)')

subplot(4,5,11), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Ca_leak, '-b'), ylabel('I_{Ca_{leak}} (pA)'), xlabel('Time (ms)')
subplot(4,5,12), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Na_leak, '-b'), ylabel('I_{Na_{leak}} (pA)'), xlabel('Time (ms)')  
subplot(4,5,13), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_K_leak, '-b'), ylabel('I_{K_{leak}} (pA)'), xlabel('Time (ms)')
subplot(4,5,14), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_Cl_leak, '-b'), ylabel('I_{Cl_{leak}} (pA)'), xlabel('Time (ms)')
    
subplot(4,5,15), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,I_ALL_leak, '-b'), ylabel('I_{ALL_{leak}} (pA)'), xlabel('Time (ms)')
subplot(4,5,16), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,J_SERCA, '-b'), ylabel('J_{SERCA} (pA)'), xlabel('Time (ms)')   
subplot(4,5,17), hold on, set(gca,'box','off','tickdir','out','fontsize',12) 
    plot(Time/1000,J_rel*2*96.4853415*0.005*1000 , '-b'), ylabel('J_{rel} (pA)'), xlabel('Time (ms)')

         


figure (320)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(4,1,1)      
    plot(Time/1000,V, '-b'), ylabel('Vm (mV)') , xlabel('Time (ms)')         
subplot(4,1,2)       
    plot(Time/1000,I_BK, '-b'), ylabel('I_{BK} (pA)'), xlabel('Time (ms)')      
subplot(4,1,3)      
    plot(Time/1000,I_CaL, '-b'), ylabel('I_{CaL} (pA)'), xlabel('Time (ms)')            
subplot(4,1,4)       
   plot(Time/1000,Ca_in*1E6, '-b'), ylabel('Ca_{Cyt} (nM)'), xlabel('Time (ms)') 

    
    figure (330)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(4,1,1)      
    plot(Time/1000,V, '-b'), ylabel('Vm (mV)') , xlabel('Time (ms)')         
subplot(4,1,2)       
    plot(Time/1000,I_KvALL, '-b'), ylabel('I_{Kv_{ALL}} (pA)'), xlabel('Time (ms)')      
subplot(4,1,3)      
    plot(Time/1000,I_Kv21, '-b'), ylabel('I_{Kv_{2.1}} (pA)'), xlabel('Time (ms)')            
subplot(4,1,4)       
    plot(Time/1000,I_Kv15, '-b'), ylabel('I_{Kv_{1.5}} (mM)'), xlabel('Time (ms)') 
    
    
       figure (340)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(3,1,1)      
       plot(Time/1000,Ca_in*1E6, '-b'), ylabel('Ca_{in} (nM)') , xlabel('Time (s)')  
    ylim([0 200])     
subplot(3,1,2)      
    plot(Time/1000,Ca_SRcen*1000, '-b'), ylabel('Ca_{SR} (\muM)'), xlabel('Time (s)')    
    ylim([0 150])
subplot(3,1,3)      
    plot(Time/1000,Ca_jun*1000, '-b'), ylabel('Ca_{Jun} (\muM)'), xlabel('Time (s)')            
ylim([0 20])  

   figure (350)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(3,1,1)      
       plot(Time/1000,Ca_in*1E6, '-b'), ylabel('Ca_{in} (nM)') , xlabel('Time (s)')  
    ylim([0 200])     
subplot(3,1,2)      
    plot(Time/1000,Na_in, '-b'), ylabel('Na_{in} (mM)'), xlabel('Time (s)')    
    ylim([0 20])
subplot(3,1,3)      
    plot(Time/1000,K_in, '-b'), ylabel('K_{in} (mM)'), xlabel('Time (s)')            
ylim([0 200])  


   figure (360)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(4,1,1)      
       plot(Time/1000,Ca_in*1E6, '-b'), ylabel('Ca_{in} (nM)') , xlabel('Time (s)')  
    ylim([0 200])     
subplot(4,1,2)      
    plot(Time/1000,Ca_SRcen*1000, '-b'), ylabel('Ca_{SR} (\muM)'), xlabel('Time (s)')    
    ylim([0 150])
    
    subplot(4,1,3)       
    plot(Time/1000,I_BK, '-b'), ylabel('STOCs (pA)'), xlabel('Time (s)')  
    subplot(4,1,4)      
    plot(Time/1000,Ca_jun*1000, '-b'), ylabel('Ca_{Jun} (\muM)'), xlabel('Time (s)')            
ylim([0 20]) 



    figure (370)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
    plot(Time/1000,V, '-b'), ylabel('Vm (mV)') , xlabel('Time (s)')  
ylim([-50 -20] )        
  xlim([485 500])


   

   figure (390)
set(gcf,'color','w');
set(groot, 'defaultLineLineWidth', 1.5)
set(groot, 'defaultAxesFontSize', 16)
subplot(3,1,1) 
plot(Time/1000,V, 'b'), ylabel('Vm (mV)') , xlabel('Time (s)')  
ylim([-60 -20] )   
xlim([485 500] ) 
subplot(3,1,2)      
    plot(Time/1000,I_CaL, '-b'), ylabel('I_{CaL} (pA)'), xlabel('Time (s)')   
    ylim([-2 0] ) 
    xlim([485 500] ) 
subplot(3,1,3)      
       plot(Time/1000,Ca_in*1E6, '-b'), ylabel('Ca_{in} (nM)') , xlabel('Time (s)')  
    ylim([0 300]) 
    xlim([485 500] ) 


