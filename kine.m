%%calculation of heating rate density using the kinetic range for all the
%%given 10 minute event of one minute  magnetometer data
clear all
clc
close all
load All_events_start_date_for_10_min_window.mat


KE1=[];KE2=[];KE3=[];KE4=[]; KE5=[]; KE6=[]; max_heating_time_average=[]; all_time_period=[];
%%{
for i=2:169
    %This is the reading of the magnetometer data
    mag_data=get_magnetometer_data(i);
    %converting the date of the mag file into number of minutes
    dates=24*60*(datenum(mag_data(1,:),mag_data(2,:),mag_data(3,:),mag_data(4,:),mag_data(5,:),mag_data(6,:))...
        -datenum(2004,1,1));
    %take the all the events dates that is contained in one second mag_data
    %file
    time_period=events_start_date(events_start_date > dates(1)  & events_start_date < dates(end));
    %Taking indices of all the dates in of events value of 1 minutes data
    indices_event_date=find(events_start_date > dates(1)  & events_start_date < dates(end));
    
    ref_q_KE1=nan(1,length(time_period));ref_q_KE2=nan(1,length(time_period));
    ref_q_KE3=nan(1,length(time_period));ref_q_KE4=nan(1,length(time_period));
    ref_q_KE5=nan(1,length(time_period));ref_q_KE6=nan(1,length(time_period));
    ref_time_average_for_maxm_q=nan(1,length(time_period)); ref_period=nan(1,length(time_period));
    
    parfor j = 1:length(time_period)

        bx = 10^(-9).*mag_data(7,dates>=time_period(j) & dates <= (time_period(j)+9));
       
        by= 10^(-9).*mag_data(8,dates>=time_period(j) & dates <= (time_period(j)+9));
        
        bz=10^(-9).* mag_data(9,dates>=time_period(j) & dates <= (time_period(j)+9));
       
        %we will calculate only those MHD heating rate values only if one
        %second data are gerater than 100 values
        if length(bx)>100
        ref_indices=indices_event_date(j);
        lat=events_lat(ref_indices);
        R= events_R(ref_indices);
    
     ref_period(j)=time_period(j);
    ref_q_KE1(j) = kinetic_HR(bx',by',bz',R,lat,1);
    %5 second
    [ average_bx2,average_by2,average_bz2 ] = average_mag_data(bx,by,bz,5);
    
     ref_q_KE2(j) = kinetic_HR(average_bx2',average_by2',average_bz2',R,lat,5);
     %10 second
     
     [ average_bx3,average_by3,average_bz3 ] = average_mag_data(bx,by,bz,10);
     
    
    ref_q_KE3(j) = kinetic_HR(average_bx3',average_by3',average_bz3',R,lat,10);
    
    %15 second
    [ average_bx4,average_by4,average_bz4 ] = average_mag_data(bx,by,bz,15);
     
   
    ref_q_KE4(j) = kinetic_HR(average_bx4',average_by4',average_bz4',R,lat,15);
    
     %20 sec
    [ average_bx5,average_by5,average_bz5 ] = average_mag_data(bx,by,bz,20);
     
    
    ref_q_KE5(j) = kinetic_HR(average_bx5',average_by5',average_bz5',R,lat,20);
    
    %30 second 
    [ average_bx6,average_by6,average_bz6 ] = average_mag_data(bx,by,bz,30);
     
    
    ref_q_KE6(j) = kinetic_HR(average_bx6',average_by6',average_bz6',R,lat,30);
    
   list_q= [ref_q_KE1(j),ref_q_KE2(j),ref_q_KE3(j),ref_q_KE4(j),ref_q_KE5(j),ref_q_KE6(j)];
   list_of_time=[1 5 10 15 20 30];
   if isnan(list_q)
       max_q=nan;
       ref_time_average_for_maxm_q(j)=nan;
   else
    max_q=max(list_q);
   
    indices_for_max_q=find(list_q==max_q);
    max_q_time=list_of_time(indices_for_max_q);
    ref_time_average_for_maxm_q(j)=max_q_time;
   end
        end
    
    end
    KE1=[KE1 ref_q_KE1];KE3=[KE3 ref_q_KE3]; KE5=[KE5 ref_q_KE5];
    KE2=[KE2 ref_q_KE2];KE4=[KE4 ref_q_KE4]; KE6=[KE6 ref_q_KE6];
    max_heating_time_average=[max_heating_time_average ref_time_average_for_maxm_q];
    all_time_period=[all_time_period ref_period];
end
%%{ 

    num_cases(1)=length(find( max_heating_time_average==1));
num_cases(2)=length(find( max_heating_time_average==5));
num_cases(3)=length(find( max_heating_time_average==10));
num_cases(4)=length(find( max_heating_time_average==15));
num_cases(5)=length(find( max_heating_time_average==20));
num_cases(6)=length(find( max_heating_time_average==30));
figure()
bar([1 5 10 15 20 30],num_cases)
xlabel('Time period of averaged data(sec)')
ylabel('Disturbed Events')
set(gcf,'color','white')
set(gca,'fontsize',30,'fontweight','bold')
%}
edges=(-25:0.25:-10);


%makign histogram for diffrent time period average of B field used for
%calculating heating rate density

%%% for one second histogram of heating density
figure()
histogram(log10(KE1),edges);
xlabel('log_{10}(q_{MHD})')
set(gca,'xtick',[-25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -10])
set(gca,'fontsize',30,'fontweight','bold')
set(gcf,'color','white')
title('One second')
%%% for 5 second histogram of heating density
figure()
histogram(log10(KE2),edges);
xlabel('log_{10}(q_{MHD})')
set(gca,'xtick',[-25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -10])
set(gca,'fontsize',30,'fontweight','bold')
set(gcf,'color','white')
title('Five second')
%%% for 10 second histogram of heating density
figure()
histogram(log10(KE3),edges);
xlabel('log_{10}(q_{MHD})')
set(gca,'xtick',[-25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -10])
set(gca,'fontsize',30,'fontweight','bold')
set(gcf,'color','white')
title('Ten second')
%%% for 15 second histogram of heating density
figure()
histogram(log10(KE4),edges);
xlabel('log_{10}(q_{MHD})')
set(gca,'xtick',[-25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -10])
set(gca,'fontsize',30,'fontweight','bold')
set(gcf,'color','white')
title('Fifteen second')



%%% for 20 second histogram of heating density
figure()
histogram(log10(KE5),edges);
xlabel('log_{10}(q_{MHD})')
set(gca,'xtick',[-25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -10])
set(gca,'fontsize',30,'fontweight','bold')
set(gcf,'color','white')
title('Twenty second')


%%% for 30 second histogram of heating density
figure()
histogram(log10(KE6),edges);
xlabel('log_{10}(q_{MHD})')
set(gca,'xtick',[-25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -10])
set(gca,'fontsize',30,'fontweight','bold')
set(gcf,'color','white')
title('Thirty second')
%}

%This is saving all the heating rate density and start date for 10 minute
%sliding window for R>7Rs
%save all_heating_rate_density_for_mva_with_VK_info.mat all_time_period MHD1 MHD2 MHD3 MHD4 MHD5 MHD6
 save all_KE_heating_rate_density_for_mva_10_second.mat all_time_period KE1 KE2 KE3 KE4 KE5 KE6

%Note heating_rate_density1 le vitaly ko calculation dincha but
%heating_rate_density le brandon ko calculation dincha

