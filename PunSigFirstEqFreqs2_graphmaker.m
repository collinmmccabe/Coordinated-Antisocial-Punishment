%This program computes expected fitness of Punishers and NonPunishers in a simple discrete generation model which incorporates some
% of the structural assumptions of Herb's agent based model

clear
clc

c = .01;                                                                   %cost of contributing to PG
b = 2*c;                                                                   %per capita benefit of public good
p = 1.5*c;                                                                   %cost of being punished
k = p;                                                                   %cost of punishing one individual by a single punisher
a = 2;                                                                     %scale parameter for cost of punishmen
q = p;                                                                     %signal cost
T = 25;                                                                    %number of interactions
n = 18;                                                                    %group sizetau = 9;                                                                   %threshold
e = .1;                                                                    %error rates                                                                %punishment threshold
r = 0.07;                                                                   %relatedness

nFreqPts = 101;
taulist = [1,5,10,15];
ntauval = length(taulist);

Freqs = ((1:nFreqPts)-1)/(nFreqPts - 1);
dq = zeros(ntauval,nFreqPts);
df = zeros(n,nFreqPts);

                         
for ip = 1:ntauval
    
    tau = taulist(ip);
    
    for iq = 1:nFreqPts

        x = Freqs(iq);
        dq(ip,iq) = WpMinusWnPunSigFirst(Freqs(iq),b,c,p,k,a,tau,q,T,n,e,r);

    end   %for iq
end       %for ip

for jt = 1:n
    
    tau = jt-1;
    
    for iq = 1:nFreqPts

        x = Freqs(iq);
        df(jt,iq) = WpMinusWnPunSigFirst(Freqs(iq),b,c,p,k,a,tau,q,T,n,e,r);

    end   %for iq
end       %for jt


[max_df,max_x] = max(df');

max_x = (max_x -1)/nFreqPts;

ntauval = n;
xequ = zeros(1,ntauval);
xeq = zeros(1,ntauval);
MutAdvan = zeros(1,ntauval);
Wav = zeros(1,ntauval);
FreqCoopGroups = zeros(1,ntauval);

for jt = 1:n                                                                
    
    tau = jt-1;
    
    dw0 = WpMinusWnPunSigFirst(0.0,b,c,p,k,a,tau,q,T,n,e,r);
    
    dw1 = WpMinusWnPunSigFirst(1.0,b,c,p,k,a,tau,q,T,n,e,r);
    
          
    if dw1 < 0 && dw0 > 0
        
        xeq(jt) = fzero(@(x) WpMinusWnPunSigFirst(x,b,c,p,k,a,tau,q,T,n,e,r),[0 1]);
        xequ(jt) = 0;
            
    elseif dw1 > 0 && dw0 < 0 
        
        xeq(jt) = 1;
        xequ(jt) = fzero(@(x) WpMinusWnPunSigFirst(x,b,c,p,k,a,tau,q,T,n,e,r),[0 1]);
        
    elseif dw1 > 0 && dw0 > 0 
        
        xeq(jt) = 1;
        xequ(jt) = 0;
        
    elseif dw1 < 0 && dw0 < 0 && max_df(jt) <= 0
        
        xeq(jt) = 0;
        xequ(jt) = 0;
        
    elseif dw1 < 0 && dw0 < 0 && max_df(jt) > 0
          
        xmax = max_x(jt);
        wmax= max_df(jt);
        
        dw1_test = WpMinusWnPunSigFirst(1,b,c,p,k,a,tau,q,T,n,e,r);
        dwInt_test = WpMinusWnPunSigFirst(xmax,b,c,p,k,a,tau,q,T,n,e,r);
        
        xeq(jt) = fzero(@(x) WpMinusWnPunSigFirst(x,b,c,p,k,a,tau,q,T,n,e,r),[xmax 1]);
        xequ(jt) = fzero(@(x) WpMinusWnPunSigFirst(x,b,c,p,k,a,tau,q,T,n,e,r),[0 xmax]);
        
    end
    
       
    
    Wav(jt) = WbarPunSigFirst(xeq(jt),b,c,p,k,a,tau,q,T,n,e,r)-WbarPunSigFirst(0,b,c,p,k,a,tau,q,T,n,e,r);    
    
    LiarCoopInvade(jt) = q - binopdf(jt,n,xeq(jt))*(b*jt/n-p) ;
    LiarDefectInvade(jt) = q - binopdf(jt,n,xeq(jt))*(b*(jt+1)/n-c)- binocdf(jt+1,n,xeq(jt))*(p - c + b/n);
    
end    



ptau = (1:n) -1;
xzero = zeros(1,n);

plot(ptau,xzero,'ko','MarkerSize',5,'MarkerFaceColor','k')
hold on
plot(ptau,xequ,'ko','MarkerSize',5,'MarkerFaceColor','w')
plot(ptau,xeq,'ko','MarkerSize',5,'MarkerFaceColor','k')
hold off
axis square
title(['b/c = ' num2str(b/c) ', k/p = ' num2str(k/p) ', a = ' num2str(a) ', q = ' num2str(q) ', n = ' num2str(n) ', e = ' num2str(e)  ', T = ' num2str(T) ', r = ' num2str(r) ] )
xlim([0 n-1])
xlabel('Threshold number of punishers (\tau)')
ylabel('equilibrium frequency')

