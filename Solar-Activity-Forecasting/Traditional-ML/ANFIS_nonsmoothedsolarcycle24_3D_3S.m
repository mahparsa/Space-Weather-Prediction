
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     test size 1800 t0 2000
UP1=3101%%%year 2007 MAY
LP1=7%%%year 1866
UP2=3189%%%yaer 2014 September
LP2=3101
monthssn=month_new;



TS = monthssn(LP1:UP1,:);
TSS= monthssn(LP2:UP2,:);
L=UP2-LP1+1

in=9;
for i=1+in:UP1-LP1+1
    M(i-in,:) = [TS(i-9,4) TS(i-6,4) TS(i-3,4) TS(i,4)];
end

for j=1+in:UP2-LP2+1
    MM(j-in,:) = [ TSS(j-9,4) TSS(j-6,4) TSS(j-3,4) TSS(j,4)];
end



DATA=[M;MM];
N=size(DATA,1);
p=size(DATA,2)-1;
U = zeros(N,1:p)
aa1=DATA(:,1:p)
aan=aa1;
% U=[aa1]
% minU=min(U)
% maxU=max(U)
% for i=1:size(U,2)
%     aan(:,i)=(U(:,i)- minU(i))./(maxU(i)-minU(i))
% end
N_TRAIN=size(M,1);
N_TEST=size(MM,1);
S=aan(1:N_TRAIN,1:p);
ST=aan(N_TRAIN+1:N,1:p);

%%%for target

t=DATA(1:N_TRAIN,p+1)
% tn=(t- min(t))./( max(t)-min(t))
Out=t;
OUT=t;
OUTT=DATA(N_TRAIN+1:N,p+1);



trnN=S;
CekN1=ST;


%%%%%%%%%%%%%%%%%%%%


Ath=max(S')'

Atht1=max(ST')'
   

SS=[Ath,S]

P=SS
PP=S
%%%%%%%%%%%%%%%%%%%%
SST=[Atht1,ST]
PT1=[SST]
PPT1=[ST]


%%%%%%%%%%%%%%%%%%%%
tic;
ll=[PP,Out];  

fismat=genfis1(ll,2);
fismatA = anfis(ll,fismat,[4000,0.00001,0.0001,0.9,1.01],[])
Ea=evalfis(PPT1,fismatA);




EET1=Ea;
ERROR1=EET1-OUTT;
PER1=mse(ERROR1);
ROOTPER1=norm(EET1-OUTT)/sqrt(length(EET1));
NMSE1=(N_TEST*PER1)/(mse(EET1-mean(OUTT))*N_TEST);
% AVE1=(1/NT)*sum((abs(ERROR1))./abs(chk_data1(:, p+1)))*100
Corr1=corrcoef(EET1,OUTT);
NRMSE1=ROOTPER1/std(OUTT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on 
plot(EET1,'g');
hold on;
plot(OUTT,'b');
