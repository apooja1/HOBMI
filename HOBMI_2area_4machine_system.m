%% HOBMI Code %%

clc 
clearvars -except 
addpath('/Users.../HOBMI_2area_4machine_system')  %add path 
addpath('/Users/..../psat')   % add path of PSAT Matlab toolox


%% Specify the model
%For kundur model
%sampling frequency 20 Hz
mdl=1.3; 
if mdl==1.05  %fault clearing at 1.05s
window=26:135;   
run('d_kundur1_mdl_05.m')
output_data=output_data(window,:);
elseif mdl==1.1  %Case A: fault clearing at 1.1s
run('d_kundur1_mdl_1_1.m')
window=37:146;   
output_data=output_data(window,:);  
elseif mdl==1.3    %Case B: fault clearing at 1.3s 
run('d_kundur1_mdl_1_3_20Hz.m') %20Hz
window=43:152;    
output_data=output_data(window,:);
end
  
% Collect the angular frequency data of four machines
w_01=output_data(:,2); 
w_02=output_data(:,3);
w_03=output_data(:,4);
w_04=output_data(:,5);

fs=1/0.05;  
Mmax=12; %set the maximum number of modes
t=output_data(:,1); %time

%% Obtaining relative angular frequency measurements w_12, w_14, w_24
% Since we do not have access to the relative angular frequency of
% measurements, we obtain those by subtracting its components. For e.g.
% w_12 is obtained by w_1-w_2.
w0012 =w_01-w_02;
w0014 =w_01-w_04;
w0024 =w_02-w_04;

%% 
Y_big=zeros(2*Mmax,52,Mmax);
delay=3;
F_big=zeros(Mmax,2*Mmax);

for mdl=1:3
 for m=2:Mmax
 d=2*m;
 j=1;
 jj=1;

    if mdl==1
X=get_delay_vector(w0012,d,delay);
    elseif mdl==2
X=get_delay_vector(w0014,d,delay);
    elseif mdl==3
X=get_delay_vector(w0024,d,delay);
    end
       
tt=t(1:length(X(:,1)));
[winv, Y] = HOBI(X',d); %uncomment for HOBMI
% [winv, Y] = sobi_2(X',d);   %uncomment for  (comparable)

Y=real(Y);
Y_big(1:d,1:40,m)=Y(:,1:40);

%% This is done for plotting after the mode identification
if ((mdl==1)&&(m==5))
    Y_20=Y;
end

if ((mdl==2)&&(m==6))
    Y_10=Y;
end

if ((mdl==3)&&(m==4))
    Y_40=Y;
end
%%

%% frequency and amplitude
for i=1:length(Y(:,1))
ifq=(instfreq(Y(i,:),fs,'Method','hilbert')); 
if rem(length(ifq),2)==0
    Nn=(2);
    fr=mean(ifq((end-Nn)/2+1:(end+Nn)/2));
else
     Nn=3;
     fr=mean(ifq((end-Nn)/2+1:(end+Nn)/2));
end
if ((0.1<abs(fr))&&(abs(fr)<3))
  freq(i,1)=abs(fr);
else 
  freq(i,1)=0;
end

end
F_big(m,:)=[sort(freq)' zeros(1,Mmax*2-length(freq))];
freq(freq==0) = [];
% freq
%damping computation 
H=[ones(length(tt),1) tt];
for i=1:length(Y(:,1))
    a(i,:)=hilbert(Y(i,:));
amp(i,:)=sqrt((a(i,:)).^2+Y(i,:).^2);
sigg=var(log(amp(i,:))); 
Kk=sigg*eye(length(amp(i,:)));
% beta=(H'*Kk^(-1)*H)^(-1)*H'*Kk^(-1)*log(amp(i,:))';  %least squares to fit the ln(amplitude)
beta=(H'*H)^(-1)*H'*log(amp(i,:))';
sigma(i,1)=real(beta(2,1));
end

[f,idx]=sort(freq);  %sigma=sort(sigma);
f_diff=diff(f);
[~,idx1]=sort(f_diff);
% f=f(idx);
s=sigma(idx1);
while 2*j<length(f)
f_avg(m-1,j)=(f(jj)+f(jj+1))/2;
sigma_avg(m-1,j)=(s(jj)+s(jj+1))/2;
d_avg(m-1,j)=sqrt(((abs(s(jj))-abs(s(jj+1)))^(2)+(abs(f(jj))-abs(f(jj+1)))^(2)));
j=j+1; jj=jj+2;
end
a=[]; amp=[]; sigma=[];
 end

f_avg
sigma_avg
d_avg
d_avg(d_avg==0) = NaN;
 %% mode identification 
[m_val1, m_det1]=min(d_avg,[],1,'omitnan');
for kk=1:length(m_det1)
f_det1(kk,1)=f_avg(m_det1(kk),kk);
sigma_det1(kk,1)=sigma_avg(m_det1(kk),kk);
end
d_avg
f_det1
sigma_det1 
m_val1'
Results=[f_det1 sigma_det1 m_val1']

%% For plotting purposes (ignore)
if mdl==1
Y_2=Y_big(6,:,2);
elseif mdl==2
Y_1=Y_big(4,:,2);
elseif mdl==3
Y_4=Y_big(m_det1(5),:,5);
end
%%
a=[]; amp=[]; sigma=[]; 
clearvars -except mdl X delay Y_big w0012 w0014 w0024 m Y_2 Y_1 Y_4 Mmax t fs tt Y_10 Y_20 Y_40


end

%% plots
tt=t(1:52);
Y_2=Y_20(4,1:52);
Y_1=Y_10(7,1:52);
Y_4=Y_40(1,1:52);
Y_2=Y_2./ 10^(floor(log10(Y_2(1))));
Y_1=Y_1./10^(floor(log10(Y_1(1))));
Y_4=Y_4./10^(floor(log10(Y_4(1))));
Y_4=Y_4/16;
Y_2=smooth(Y_2);
Y_1=smooth(Y_1);
Y_4=smooth(Y_4);


figure()
plot(tt,Y_1)
grid on
hold on 
plot(tt,Y_2)
plot(tt,Y_4)
% plot(tt,Y_24(2,:))
% plot(tt,w_04(1:length(tt)))
hold off 
xlabel('Time (s)')
ylabel('\omega')
legend('y_1','y_2','y_3')
legend Box off
% title('')
set(gcf,'color','w')


%% Toy model 2

clc 
clear all 

addpath('/Users/poojaalgikar/Downloads/stuff/BSS') 
addpath('/Users/poojaalgikar/Downloads/stuff/psat')
addpath('/Users/poojaalgikar/Downloads/stuff/psat/bmid_pkg') 


M=2;
tt=1:0.01:10;tt=tt';
T=length(tt);
for t=1:T
if (tt(t)<=2)||(tt(t)<=10)&&(tt(t)>=6)
    x1(t,1)=exp(-0.01*tt(t))*cos(8*tt(t));
elseif (tt(t)<6)&&(tt(t)>2)
    x1(t,1)=0;

end
end 



for t=1:length(tt)
if (tt(t)<=10)&&(tt(t)>=6)
    x2(t,1)=0.6*exp(-0.03*tt(t))*cos(17*tt(t));
elseif (tt(t)<6)&&(tt(t)>0)
    x2(t,1)=0;

end
end


% alpha=0.9;  %Gaussian
%    Gumebel alpha=2;
alpha=10;

U0 = copularnd('Clayton',alpha,length(tt));        
s1=U0(:,1);s2=U0(:,2);

x1=x1+0.01*s1;
x2=x2+0.01*s2;


% x1=x1(501:end); x2=x2(501:end);   
% tt=tt(501:end);

xx1=5*x1+8*x2;  xx2=4*x1-2*x2;

figure()
plot(tt,x1) 
hold on 
plot(tt,x2) 
hold off 
legend('x_1','x_2')
title('The original signal')

figure()
plot(tt,xx1) 
hold on 
plot(tt,xx2) 
hold off 
legend('xx_1','xx_2')
title('The observed signal')


fs=100;
% Fo1=8/(2*pi); Fo2=17/(2*pi);  


%frequency computation
disp('Original system')
ifq1=instfreq(x1,fs,'Method','hilbert'); fo1=mean(ifq1) 
ifq2=instfreq(x2,fs,'Method','hilbert'); fo2=mean(ifq2)


%Amplitude computation
 a1 = hilbert(x1); a2 = hilbert(x2); a=[a1';a2'];

ifq=[ifq1';ifq2'];x=[x1';x2'];    
H=[ones(length(tt),1) tt(1:end,1)];
for i=1:M
amp(i,:)=sqrt((a(i,:)).^2+x(i,:).^2);
beta(i,:)=(H'*H)^(-1)*H'*log(amp(i,:))';  %least squares to fit the ln(amplitude)
end
% sigo1=-0.01; sigo2=-0.03;
sigmao1=real(beta(1,2));
sigmao2=real(beta(2,2));

X=[xx1 xx2];

%% HOBMI 
disp('HOBI')
[winv, Y] = HOBI(X');
Y=real(Y);
%frequency computation
for i=1:M
ifq1=instfreq(Y(i,:),fs,'Method','hilbert'); 
f1(i,1)=mean(ifq1);
end
f_D1=f1(1) ;
f_D2=f1(2) ;
%Amplitude computation   
H=[ones(length(tt),1) tt(1:end,1)];
for i=1:M
    a(i,:)=hilbert(Y(i,:));
amp(i,:)=sqrt((a(i,:)).^2+Y(i,:).^2);
beta(i,:)=(H'*H)^(-1)*H'*log(amp(i,:))';  %least squares to fit the ln(amplitude)
end
% sigo1=-0.01; sigo2=-0.03;
sigma_D1=real(beta(1,2))
sigma_D2=real(beta(2,2))

figure()
plot(tt,Y(1,:))
hold on 
plot(tt,x1)
hold off 
legend('the estimated signal 1 (SOBI-I)','x_1')
title('The estimated source signal 1 from SOBI-I')
set(gcf,'color','w')

figure()
plot(tt,Y(2,:))
hold on 
plot(tt,x2)
hold off 
title('The estimated source signal 2 from SOBI-I')
legend('the estimated signal 2 (SOBI-I)','x_2')
set(gcf,'color','w')


%%SOBI  
disp('SOBI')
[winv, Y] = sobi_2(X');
Y=real(Y);

%frequency computation
for i=1:M
ifq1=instfreq(Y(i,:),fs,'Method','hilbert'); 
f1(i,1)=mean(ifq1);
end
f_S1=f1(1);
f_S2=f1(2);
%Amplitude computation   
H=[ones(length(tt),1) tt(1:end,1)];
for i=1:M
    a(i,:)=hilbert(Y(i,:));
amp(i,:)=sqrt((a(i,:)).^2+Y(i,:).^2);
beta(i,:)=(H'*H)^(-1)*H'*log(amp(i,:))';  %least squares to fit the ln(amplitude)
end
sigo1=-0.01; sigo2=-0.03;
sigma_S1=real(beta(1,2))
sigma_S2=real(beta(2,2))

figure()
plot(tt,Y(1,:))
hold on 
plot(tt,x1)
hold off 
legend('the estimated signal 1 (SOBI)','x_1')
title('The estimated source signal 1 (SOBI)')
set(gcf,'color','w')

figure()
plot(tt,Y(2,:))
hold on 
plot(tt,x2)
hold off 
title('The estimated source signal 2 (SOBI)')
legend('the estimated signal 2 (SOBI)','x_2')
set(gcf,'color','w')
%% For table construction
AAA=[[1;2] [fo1;fo2] [sigmao1;sigmao2] [f_D1;f_D2] [sigma_D1;sigma_D2] [f_S1; f_S2] [sigma_S1;sigma_S2]];



