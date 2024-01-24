clc;
clear all;
close all;

% Sampling frequency and sampling period
fs=61000
ts=1/fs
% Duration
td=0.1
tau=10*(10^(-3))
snr=3.5 %given snr
a=300000
N=(td/ts)+1; % Number of observed samples, obtained using Ts(N-1)=Td
t=linspace(0,td,N).'; % Generate Nx1 vector of time in seconds
Mobs=20;
stat_test=zeros(Mobs,1);

n0=ceil(tau/ts)%delay
t1=linspace(0,td-tau,N-n0).' ;
chirpsignal=cos(2*pi*(a/2)*(t1.^2)); %transmitted chirp signal
zero_vector=zeros(n0, 1); %zeros matrix of dimentions n0x1
s=[zero_vector;chirpsignal]; %chirp in matrix

%loading given covariance and inverse of data
load('covariance.mat') ;
%recieved data set
dataset_1=load('IPR_Set1.mat') ;
dataset_2=load('IPR_Set2.mat') ;
x_1 = cell2mat(table2array(cell2table(struct2cell(dataset_1))));
x_2 = cell2mat(table2array(cell2table(struct2cell(dataset_2))));
set1 = load("IPR_Set1",'-mat')    % variables: 'set1.Received_Data'
set2 = load("IPR_Set2",'-mat')     % variables: 'set2.Received_Data'
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%----------------1st question(a)---------------------------
%choose a value of pfa and find the threshold using the pfa 
pfa=0.009
qinv_pfa=qfuncinv(pfa);
temp_var=(sqrt((s.')*(Cw_inv)*(s)));
threshold_value=temp_var*qinv_pfa
%---------------1st question(b)----------------------------
%choose a minimum value of pfa and find maximised value of pd
pd=qfunc((qinv_pfa)-temp_var)

%--------------2nd question---------------------------------
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

count1  = 0;
count2  = 0;
for iter= 1:Mobs
    x1 = set1.Received_Data(:,iter);
    x2 = set2.Received_Data(:,iter);

    L1(iter)=x1.'*Cw_inv*s;
    L2(iter)=x2.'*Cw_inv*s;

    if(L1(iter) > threshold_value)
        count1  = count1 + 1;
    end

    if(L2(iter) > threshold_value)
        count2  = count2 + 1;
    end

end
figure
plot(0:Mobs,threshold_value*ones(Mobs+1),'LineWidth',2,'Color',"black")
hold on
legend('Threshold')
scatter(1:Mobs, L1(1:Mobs),'DisplayName','set 1: Lx','MarkerFaceColor', "red")
scatter(1:Mobs, L2(1:Mobs),'DisplayName','set 2: Lx','MarkerFaceColor', "green")
xlabel("Oberservation Number")
ylabel("Test Statistic")
hold off
legend
title("Detection on 2 datasets")

tempvarv1=(sqrt(x1.'*Cw_inv*s))
d1=real(tempvarv1)

pd1=qfunc(qinv_pfa-d1)
temp_varv2=(sqrt(x2.'*Cw_inv*s))
d2=real(temp_varv2)
pd2=qfunc(qinv_pfa-d2)
pdest1=count1/Mobs
pdest2=count2/Mobs

%---------------3rd question--------------------------------
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%plot an Roc curve of pfa Vs pd with a constant SNR
PFA=linspace(0,0.1,100);
for i=1:length(PFA)
    PD(i)=qfunc(qfuncinv(PFA(i))-sqrt((s.') * Cw_inv *s));
end
figure(2)
plot(PFA,PD,LineWidth=2)
title("ROC curve of chirp for SNR=3.5db")
xlabel("probability of false alarm (Pfa)")
ylabel("probabiliy of detection (Pd)")
%-------------4th question(a)---------------------------------
%design a pre-whitening transformation; derive matrix D theoritically
%-------------4th question(b)---------------------------------
%compute D and determine detector performance
%make new x and s
 %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[V,A]=eig(Cw);
D=(sqrt(A^(-1)))*(V.');

%threshold
threshold2=(sqrt((s.')*(D.')*(D)*(s))).*(qinv_pfa);

                count3  = 0;
count4  = 0;
for iter= 1:Mobs
    x3 = set1.Received_Data(:,iter);
    x4 = set2.Received_Data(:,iter);

    L3(iter)=x3.'*(D.')*(D)*s;
    L4(iter)=x4.'*(D.')*(D)*s;

    if(L3(iter) > threshold2)
        count3 = count3 + 1;
    end

    if(L4(iter) > threshold2)
        count4  = count4 + 1;
    end

end
figure(3)
plot(0:Mobs,threshold2*ones(Mobs+1),'LineWidth',2,'Color',"black")
hold on
legend('threshold')
scatter(1:Mobs, L3(1:Mobs),'DisplayName','set 1: Lx','MarkerEdgeColor',"[0 .5 .5]",'MarkerFaceColor',"magenta")
scatter(1:Mobs, L4(1:Mobs),'DisplayName','set 2: Lx','MarkerEdgeColor',"[0 .5 .5]",'MarkerFaceColor',"cyan")
xlabel("Oberservation Number")
ylabel("Test Statistic")
hold off
legend
title("Detection performance on datasets using D matrix")
%prob of detection
prob_detection1=qfunc(qinv_pfa-(sqrt((s.')*(D.')*(D)*(s))))
%prob_detection2=qfunc(qinv_pfa-sqrt(x1.'*(D.')*(D)*s))
%prob_detection3=qfunc(qinv_pfa-sqrt(x2.'*(D.')*(D)*s))
pdest3=count3/Mobs
pdest4=count4/Mobs

%-------------5th question---------------------------------
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

es=s.'*s;
vmin=min(diag(A));
min=find(A==vmin);
n1=ceil(min/6101);
sn=sqrt(es)*V(:,n1);
thresh_new=sqrt(sn.'*Cw_inv*sn)*qinv_pfa
prob_det_new=qfunc(qinv_pfa-(sqrt((sn.')*(Cw_inv)*(sn))))
 count5  = 0;

for iter= 1:Mobs
   
    L5(iter)=(sn.')*(Cw_inv)*sn;

    if(L5(iter) > thresh_new)
        count5= count5 + 1;
    end

end

%figure()
%plot(0:Mobs,thresh_new*ones(Mobs+1),'LineWidth',2,'Color',"red")
%hold on
%legend('threshold')
%scatter(1:Mobs, L5(1:Mobs),'DisplayName','set 1: Lx','MarkerEdgeColor',"[0 .5 .5]",'MarkerFaceColor',"#EDB120")
%xlabel("Oberservation Number")
%ylabel("Test Statistic")
%hold off
%legend

PFA=linspace(0,0.1,100);
for i=1:length(PFA)
    PD(i)=qfunc(qfuncinv(PFA(i))-sqrt((s.') * Cw_inv *s));
end
figure()
plot(PFA,PD,'DisplayName','initial roc','LineWidth',2,Color='red')
hold on
PFA=linspace(0,0.1,100);
for j=1:length(PFA)
    prob_det_new (j)=qfunc(qfuncinv(PFA(j))-sqrt((sn.') * Cw_inv *sn));
end
plot(PFA,prob_det_new,'DisplayName','new signal roc','LineWidth',2,Color='black')
hold off
title("ROC curve for chirp signal Vs new signal")
xlabel("probability of false alarm (Pfa)")
ylabel("probabiliy of detection (Pd)")
legend
%----------------5(c)------------------------------------------------------------------

Lwin=100; h=hamming(Lwin); overlap=Lwin/2;
[Wn1,freq1]=pwelch(x1,h,overlap, [ ], fs,'centered');
[Wn2,freq2]=pwelch(s,h,overlap, [ ], fs,'centered');

figure()
plot(freq1,Wn1)
title("noise signal")
figure()
plot(freq2,Wn2)
title("chirp signal")
figure()
[Wn3,freq3]=pwelch(sn,h,overlap,[],fs,'centered');
figure()
plot(10*log10(pwelch(sn,h,overlap,[],fs,'centered')))
xlabel('hz/sample')
title("psd of new signal")

Nfft=2.^(nextpow2(N));
S=abs(fftshift(fft(sn,Nfft)));
freq3=linspace(-fs/2,fs/2,Nfft);
figure()
plot(S,freq3)






