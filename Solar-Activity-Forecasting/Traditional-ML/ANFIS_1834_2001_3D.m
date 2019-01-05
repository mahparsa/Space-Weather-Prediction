%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Monthly Sunspot Number
%%GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOODDDDDDDDDDDDDDDDDDDDDDD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     test size 1800 t0 2000
UP1=2030%%%year 1949
LP1=1030%%%year 1866
UP2=3030%%%yaer 1965
LP2=2030
monthssn=month_new;
TS = monthssn(LP1:UP1,:);
TSS= monthssn(LP2:UP2,:);
L=UP2-LP1+1

in=3;
for i=in+1:UP1-LP1+1
    M(i-in,:) = [TS(i-3,4)  TS(i-2,4) TS(i-1,4) TS(i,4)];
end

for j=1+in:UP2-LP2+1
    MM(j-in,:) = [ TSS(j-3,4) TSS(j-2,4) TSS(j-1,4) TSS(j,4)];
end


DATA=[M;MM];
N=size(DATA,1);
p=size(DATA,2)-1;
U = zeros(N,1:p)
aa1=DATA(:,1:p)
aan=aa1;

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



PP=S
%%%%%%%%%%%%%%%%%%%%

PPT1=[ST]


%%%%%%%%%%%%%%%%%%%%
tic;
ll=[PP,Out];  

fismat=genfis1(ll,2);
fismatA = anfis(ll,fismat,[1000,0.01,0.0001,0.9,1.01],[])
EET1=evalfis(PPT1,fismatA);

ERROR1=EET1-OUTT;
PER1=mse(ERROR1);
ROOTPER1=norm(EET1-OUTT)/sqrt(length(EET1));


NMSE1=((N_TEST-1)/N_TEST)*(var(ERROR1)/var(OUTT));
% AVE1=(1/NT)*sum((abs(ERROR1))./abs(chk_data1(:, p+1)))*100
Corr1=corrcoef(EET1,OUTT);
NRMSE1=ROOTPER1/std(OUTT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on 
plot(EET1,'g');
hold on;
plot(OUTT,'b');
