% function [NMSE]=ANFIS_AE_1DAY(h,NM,IT)
% load('ONE_MINUTE_1992');
TRAIN_ONE=[];
for i=145:168
    TRAIN_ONE=[TRAIN_ONE,TRAIN1(i,:) ];
end
TEST_ONE=[];
for ii=1:24
    TEST_ONE=[TEST_ONE,TEST(ii,:) ];
end

UP1= size(TRAIN_ONE,2);
LP1=1
LP2=1
UP2=size(TEST_ONE,2)
TS=[]
TSS=[]
TS = TRAIN_ONE';
TSS= TEST_ONE';
M=[]
MM=[]
for h=1:1
in=2;
A=[];
AA=[];
for i=in*h+1:UP1-LP1+1
    A(i-in*h,:) = [  TS(i-2*h,1) TS(i-h,1) TS(i,1)];
end
for j=1+in*h:UP2-LP2+1
    AA(j-in*h,:) = [ TSS(j-2*h,1) TSS(j-h,1) TSS(j,1)];
end
p=2
M=A;
MM=AA;
p=2;
N=size(M,1);
NT=size(MM,1);
Nt=N+NT;
A=[M];
  %%p show the number of input data 
U = zeros(N, p+1);
aa1=A;
U=[aa1];
minU=min(U);
maxU=max(U);
aan=[];
for i=1:size(U,2)
    aan(:,i)=(U(:,i)- minU(i))./(maxU(i)-minU(i));
end  
AT=[MM];
  %%p show the number of input data 
aa=AT;
UT=[aa];
minUT=min(UT);
maxUT=max(UT);
aanT=[];
for i=1:size(UT,2)
    aanT(:,i)=(UT(:,i)- minUT(i))./(maxUT(i)-minUT(i));
end  
Pn=aan;
PnT=aanT;
chk_data=AT(:,p+1);
trnN=Pn(1:N,1:p);
CekN=PnT(1:NT,1:p);
S1=trnN;
ST=CekN;
OUTT=PnT(:,p+1);
Ea=zeros(1,NT);
Eo=zeros(1,NT);
Out=Pn(1:N,p+1);
tic
for i=1:N
    Ath1(i)=max(S1(i,:));
end
    
fismat=genfis1([trnN,Out],2);
fismat2 = anfis([trnN,Out],fismat,[1000,0.000001,0.01,0.9,1.01],[]);
OO= evalfis(CekN, fismat2);
toc
PERT=mse(OUTT-OO);
ROOTPERT=norm(OO-OUTT)/sqrt(length(OO));
ET=OO;
EET =ET.*(maxU(p+1)-minU(p+1))+minU(p+1); 
ERROR=EET-chk_data;
PER=mse(ERROR);
ROOTPER=norm(EET-chk_data)/sqrt(length(EET));
NMSE(h)=(NT*PER)/(mse(chk_data-mean(chk_data))*NT);
end
AVE=(1/NT)*sum( abs(ERROR)./abs(chk_data) )*100;
Corr=corrcoef(EET,chk_data);
plot(chk_data)
hold on 
plot(EET,'--r')
% 
% ROOTPERT=norm(OO-TESTDAT(:,5))/sqrt(length(OO))
