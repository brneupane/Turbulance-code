%This file find out the change in B_t for given disturbed events
clear all
clc
close all
load minimum_variance_value_for_given_events_of_10min_using_one_min_mag_data_with_CS_crossing_id.mat
%loading all the mag data 
fileID = fopen('All_magnetic_field_data_for_magnetosphere.txt','r');
formatSpec = '%d %d %d %d %d %f %f %f %f %f %f %f';
sizeA = [12 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
 mag_Data=A;
mag_yr=mag_Data(1,:);mag_mon=mag_Data(2,:);mag_day=mag_Data(3,:);
mag_hr=mag_Data(4,:); mag_min=mag_Data(5,:); mag_sec=mag_Data(6,:);
% %%Taking all the magnetic field 
Br=mag_Data(7,:); Bt=mag_Data(8,:);  Bp=mag_Data(9,:); 
% % All the mag location data 
mag_R=mag_Data(10,:); mag_LT=mag_Data(11,:); mag_lat=mag_Data(12,:);

% % Converting all the mag date into numbre of minute using (2004,1,1) as the reference date  
mag_date_in_min=24*60*(datenum(mag_yr,mag_mon,mag_day,mag_hr,mag_min,mag_sec)) - 24*60*datenum(2004,1,1);


%checking all the neg_Bt disturb events events
 plasmoid=nan(1,length(mva_date));
for i=1:length(mva_date);
    I=find(mag_date_in_min==mva_date(i));
    ref_Bt=Bt(I:I+10);
    I_neg_Bt=find(ref_Bt<0);
    if ~isempty(I_neg_Bt)
        plasmoid(i)=true;
    else
        plasmoid(i)=false;
    end
end

disturb_events=find(CS_id==true & mva_value>0.25);
good_events=find(CS_id==true &  plasmoid==true & mva_value>0.25);
colormap(jet(12))
scatter(mva_R(good_events),mva_lat(good_events),200,mva_LT(good_events),'filled')
xlabel('R_s');ylabel('Lat')
h1=colorbar;
ylabel(h1,'Local Time')
set(gca,'fontsize',40);set(gcf,'color','w');
%%%% looking at the events for the -ve BT inside 23 Rs for quiet case 
quiet_events=find(((CS_id==true & mva_value<0.25) |CS_id==false) & mva_R<23);
quiet_good_events=find(((CS_id==true & mva_value<0.25) |CS_id==false) & mva_R<23 &  plasmoid==true);
figure()
colormap(hsv)
scatter(mva_R(quiet_good_events),mva_lat(quiet_good_events),200,mva_LT(quiet_good_events),'filled')
xlabel('R_s');ylabel('Lat')
h1=colorbar;
ylabel(h1,'Local Time')
set(gca,'fontsize',40);set(gcf,'color','w');