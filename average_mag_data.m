function [ average_bx,average_by,average_bz ] = average_mag_data(bx,by,bz,time_interval)
%bx by bz are the 15 minute one second data 
n=floor(length(bx)./time_interval);
k=1; 
for i=1:n
    average_bx(i)=mean(bx(k:k+time_interval-1));
    average_by(i)=mean(by(k:k+time_interval-1));
    average_bz(i)=mean(bz(k:k+time_interval-1));
k=time_interval+k;
end

end

