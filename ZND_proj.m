clear; clc;

mech = 'h2air_highT.cti';   
gas = Solution(mech); 
gas1 = Solution(mech);
nsp = nSpecies(gas);
P1 = [];
gas2 = Solution(mech); 
gas3 = Solution(mech);
nsp = nSpecies(gas2);
T1 = [];
gas4 = Solution(mech); 
gas5 = Solution(mech);
nsp = nSpecies(gas4);
phi = [];


% find hydrogen, nitrogen, and oxygen indices
ih2 = speciesIndex(gas,'H2');
io2  = speciesIndex(gas,'O2');
in2  = speciesIndex(gas,'N2');
x = zeros(nsp, 1);
x(ih2,1) = 2.0;
x(io2,1) = 1.0;
x(in2,1) = 3.76;

T0 = 300;
P0 = 100000;

% fig_num & fname are for 'znd' - use '0' for no output
% plots: 1 = make plots, otherwise no plots
fig_num = 0;
fname = 0;
plots = 1;

display('Initial Conditions')
display('   Pressure = 100000 Pa; Temperature = 300 K; Equivalence Ratio = 1')

% pressure loop
npoints=20;
   disp(['For ', num2str(npoints), ' values of P1'])
for i = 1:npoints
   P1(i) = 100000 + 900000/(npoints-1)*(i-1);
   disp([' ', num2str(i),  ': P1 = ', num2str(P1(i)), '; T1 = 300; phi = 1'])
   set(gas,'Temperature',T0,'Pressure',P1(i),'MoleFractions',x);
   
   %%%Constant Volume Explosion Data%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   [cj_speed(i)] = CJspeed(P1(i),T0,x,mech,0); 
   set(gas,'Temperature',T0,'Pressure',P1(i),'MoleFractions',x);
   [gas] = PostShock_fr(cj_speed(i), P1(i), T0, x, mech);
   

   %%%%%ZND Detonation Data%%%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   set(gas1, 'T', T0, 'P', P1(i), 'X', x);
   [gas] = PostShock_fr(cj_speed(i), P1(i), T0, x, mech);
   Ts(i) = temperature(gas); %frozen shock temperature   
   Ps(i) = pressure(gas); %frozen shock pressure
   % SOLVE ZND DETONATION ODES
   [ZNDout] = znd(gas,gas1,fig_num,cj_speed(i),fname);
   ind_time_ZND(i) = ZNDout.ind_time_ZND;
   ind_len_ZND(i) = ZNDout.ind_len_ZND;
   exo_time_ZND(i) = ZNDout.exo_time_ZND;
   exo_len_ZND(i) = ZNDout.exo_len_ZND;
   tsteps = size(ZNDout.T,2);
   Tf_ZND(i) = ZNDout.T(tsteps);
   
   %%Calculate CJstate Properties%%%
   [gas] = PostShock_eq(cj_speed(i),P1(i),T0,x, mech);
   T2(i) = temperature(gas);
   P2(i) = pressure(gas);
   rho2(i) = density(gas);
      
   
end

% temperature loop
npoints=20;
   disp(['For ', num2str(npoints), ' values of T1'])
for i = 1:npoints
   T1(i) = 300 + 1200/(npoints-1)*(i-1);
   disp([' ', num2str(i+20), ': P1 = 100000 Pa',  '; T1 = ', num2str(T1(i)), '; phi = 1'])
   set(gas2,'Temperature',T1(i),'Pressure',P0,'MoleFractions',x);
   
   %%%Constant Volume Explosion Data%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   [cj_speed_t(i)] = CJspeed(P0,T1(i),x,mech,0); 
   set(gas2,'Temperature',T1(i),'Pressure',P0,'MoleFractions',x);
   [gas2] = PostShock_fr(cj_speed_t(i), P0, T1(i), x, mech);
   

   %%%%%ZND Detonation Data%%%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   set(gas3, 'T', T1(i), 'P', P0, 'X', x);
   [gas2] = PostShock_fr(cj_speed_t(i), P0, T1(i), x, mech);
   Ts_t(i) = temperature(gas2); %frozen shock temperature   
   Ps_t(i) = pressure(gas2); %frozen shock pressure
   % SOLVE ZND DETONATION ODES
   [ZNDout] = znd(gas2,gas3,fig_num,cj_speed_t(i),fname);
   ind_time_ZND_t(i) = ZNDout.ind_time_ZND;
   ind_len_ZND_t(i) = ZNDout.ind_len_ZND;
   exo_time_ZND_t(i) = ZNDout.exo_time_ZND;
   exo_len_ZND_t(i) = ZNDout.exo_len_ZND;
   tsteps = size(ZNDout.T,2);
   Tf_ZND_t(i) = ZNDout.T(tsteps);
   
   %%Calculate CJstate Properties%%%
   [gas2] = PostShock_eq(cj_speed_t(i),P0,T1(i),x, mech);
   T2_t(i) = temperature(gas2);
   P2_t(i) = pressure(gas2);
   rho2_t(i) = density(gas2);
      
   
end

% phi loop
npoints=20;
   disp(['For ', num2str(npoints), ' values of phi'])
for i = 1:npoints
   phi(i) = 0.5 + 1.5/(npoints-1)*(i-1);
   disp([' ', num2str(i+40), ': P1 = 100000; T1 = 300; phi = ', num2str(phi(i))])
   x = zeros(nsp, 1);
   x(ih2,1) = 2*phi(i);
   x(io2,1) = 1.0;
   x(in2,1) = 3.76;
   
   set(gas4,'Temperature',T0,'Pressure',P0,'MoleFractions',x);
   
   %%%Constant Volume Explosion Data%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   [cj_speed_f(i)] = CJspeed(P0, T0, x, mech, 0);   
   [gas4] = PostShock_fr(cj_speed_f(i), P0, T0, x, mech);
   % SOLVE CONSTANT VOLUME EXPLOSION ODES
   [CVout] = explosion(gas4,fig_num);
   exo_time_CV_f(i) = CVout.exo_time;
   ind_time_CV_f(i) = CVout.ind_time;

   %%%%%ZND Detonation Data%%%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   set(gas5, 'T', T0, 'P', P0, 'X', x);
   [gas4] = PostShock_fr(cj_speed_f(i), P0, T0, x, mech);
   Ts_f(i) = temperature(gas4); %frozen shock temperature   
   Ps_f(i) = pressure(gas4); %frozen shock pressure
   % SOLVE ZND DETONATION ODES
   [ZNDout] = znd(gas4,gas5,fig_num,cj_speed_f(i),fname);   
   exo_time_ZND_f(i) = ZNDout.exo_time_ZND;
   exo_len_ZND_f(i) = ZNDout.exo_len_ZND;
   ind_time_ZND_f(i) = ZNDout.ind_time_ZND;
   ind_len_ZND_f(i) = ZNDout.ind_len_ZND;
   tsteps = size(ZNDout.T,2);
   Tf_ZND(i) = ZNDout.T(tsteps);
   
   %%Calculate CJstate Properties%%%
   [gas4] = PostShock_eq(cj_speed_f(i),P0, T0,x, mech);
   T2_f(i) = temperature(gas4);
   P2_f(i) = pressure(gas4);
   rho2_f(i) = density(gas4);
   
end

if(plots==1)
    % make plots
    close all;
    fontsize=12;
    
    % pressure plots
    figure;
    subplot(2,2,1);
    plot(P1(:),T2(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Temperature (K)','FontSize',fontsize);
    title('Post CJ State Temperature','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,2);
    plot(P1(:),P2(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Pressure (Pa)','FontSize',fontsize);
    title('Post CJ State Pressure','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,3);
    plot(P1(:),rho2(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Density (kg/m^3)','FontSize',fontsize);
    title('Post CJ State Density','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,4);
    plot(P1(:),cj_speed(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Ucj (m/s)','FontSize',fontsize);
    title('CJ Velocity','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    
    % temperature plots
    figure;
    subplot(2,2,1);
    plot(T1(:),T2_t(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('Temperature (K)','FontSize',fontsize);
    title('Post CJ State Temperature','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,2);
    plot(T1(:),P2_t(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('Pressure (Pa)','FontSize',fontsize);
    title('Post CJ State Pressure','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,3);
    plot(T1(:),rho2_t(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('Density (kg/m^3)','FontSize',fontsize);
    title('Post CJ State Density','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,4);
    plot(T1(:),cj_speed_t(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('Ucj (m/s)','FontSize',fontsize);
    title('CJ Velocity','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);
    
    
    % equivalence ratio plots
    figure;
    subplot(2,2,1);
    plot(phi(:),T2_f(:),'k');
    xlabel('Equivalence Ratio','FontSize',fontsize);
    ylabel('Temperature (K)','FontSize',fontsize);
    title('Post CJ State Temperature','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,2);
    plot(phi(:),P2_f(:),'k');
    xlabel('Equivalence Ratio','FontSize',fontsize);
    ylabel('Pressure (Pa)','FontSize',fontsize);
    title('Post CJ State Pressure','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,3);
    plot(phi(:),rho2_f(:),'k');
    xlabel('Equivalence Ratio','FontSize',fontsize);
    ylabel('Density (kg/m^3)','FontSize',fontsize);
    title('Post CJ State Density','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,4);
    plot(phi(:),cj_speed_f(:),'k');
    xlabel('Equivalence Ratio','FontSize',fontsize);
    ylabel('Ucj (m/s)','FontSize',fontsize);
    title('CJ Velocity','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);
    
    

    % Ts and Tf plots
    figure;
    subplot(3,2,1);
    plot(P1(:),Ts(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('T_s (K)','FontSize',fontsize);
    title('Frozen shock Temperature for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(3,2,2);
    plot(P1(:),Tf_ZND(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('T_f (K)','FontSize',fontsize);
    title('Final Temperature for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(3,2,3);
    plot(T1(:),Ts_t(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('T_s (K)','FontSize',fontsize);
    title('Frozen shock Temperature for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(3,2,4);
    plot(T1(:),Tf_ZND_t(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('T_f (K)','FontSize',fontsize);
    title('Final Temperature for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);
    
    subplot(3,2,5);
    plot(phi(:),Ts_f(:),'k');
    xlabel('Equivalence Ratio','FontSize',fontsize);
    ylabel('T_s (K)','FontSize',fontsize);
    title('Frozen shock Temperature for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(3,2,6);
    plot(phi(:),Tf_ZND(:),'k');
    xlabel('Equivalence Ratio','FontSize',fontsize);
    ylabel('T_f (K)','FontSize',fontsize);
    title('Final Temperature for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);


end
