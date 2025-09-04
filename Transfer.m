function [TXX, TXY, TYX, TYY, TZX, TZY, fc] = Transfer(x, y, X, Y, Z)

fs=1;
N=size(x,1);
t=(0:N-1)/fs;%Sampling time

xr = highpass(x,1/940,fs);
yr = highpass(y,1/940,fs);

xs = highpass(X,1/940,fs);
ys = highpass(Y,1/940,fs);
zs = highpass(Z,1/940,fs);

%% ISTF calculation
[wxr,fc]=cwt(xr,'amor');
wxrj = conj (wxr);
wyr=cwt(yr,'amor');% 
wyrj = conj (wyr);
wxs=cwt(xs,'amor');% 
wys=cwt(ys,'amor');% 
wzs=cwt(zs,'amor');% 
wxsj = conj (wxs);
wysj = conj (wys);
wzsj = conj (wzs);

period=zeros(size(wxr,1),7);
period(:,1)=1./fc(:,1);
X = 0:0.1:9.9;
CC=4*exp(-X/2)+1;

for i=1:size(period,1)
    period(i,2)=period(i,1)*CC(i);
end
for i=1:size(period,1)
    period(i,3)=ceil(period(i,2));
    period(i,4)=ceil((period(i,3))/2);
    period(i,5)=7200-period(i,4);
end
wxsxr=wxs.*wxrj;
wyryr=wyr.*wyrj;
wxsyr=wxs.*wyrj;
wxrxr=wxr.*wxrj;
wyrxr=wyr.*wxrj;
wxryr=wxr.*wyrj;
wysxr=wys.*wxrj;
wysyr=wys.*wyrj;
wzsyr=wzs.*wyrj;
wzsxr=wzs.*wxrj;
wysys=wys.*wysj;
wxsxs=wxs.*wxsj;
wzszs=wzs.*wzsj;

for i=1:size(wxr,1)
    for j=1:period(i,4)-1
        sxsxr(i,j)=sum(wxsxr(i,1:period(i,3))); 
        syryr(i,j)=sum(wyryr(i,1:period(i,3))); 
        sxryr(i,j)=sum(wxryr(i,1:period(i,3))); 
        sxsyr(i,j)=sum(wxsyr(i,1:period(i,3))); 
        sxrxr(i,j)=sum(wxrxr(i,1:period(i,3))); 
        syrxr(i,j)=sum(wyrxr(i,1:period(i,3))); 
        sysxr(i,j)=sum(wysxr(i,1:period(i,3))); 
        sysyr(i,j)=sum(wysyr(i,1:period(i,3))); 
        szsxr(i,j)=sum(wzsxr(i,1:period(i,3))); 
        szsyr(i,j)=sum(wzsyr(i,1:period(i,3))); 
        sxsxs(i,j)=sum(wxsxs(i,1:period(i,3))); 
        sysys(i,j)=sum(wysys(i,1:period(i,3))); 
        szszs(i,j)=sum(wzszs(i,1:period(i,3))); 
    end 
end
for i=1:size(wxr,1)
    for j=period(i,4):period(i,5)
        sxsxr(i,j)=sum(wxsxr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        syryr(i,j)=sum(wyryr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sxryr(i,j)=sum(wxryr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sxsyr(i,j)=sum(wxsyr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sxrxr(i,j)=sum(wxrxr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        syrxr(i,j)=sum(wyrxr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sysxr(i,j)=sum(wysxr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3)));  
        sysyr(i,j)=sum(wysyr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        szsxr(i,j)=sum(wzsxr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3)));  
        szsyr(i,j)=sum(wzsyr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sxsxs(i,j)=sum(wxsxs(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sysys(i,j)=sum(wysys(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        szszs(i,j)=sum(wzszs(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
    end 
end

for i=1:size(wxr,1)
    for j=period(i,5)+1:7200
        sxsxr(i,j)=sum(wxsxr(i,(7200-period(i,3)+1:7200)));
        syryr(i,j)=sum(wyryr(i,(7200-period(i,3)+1:7200)));
        sxryr(i,j)=sum(wxryr(i,(7200-period(i,3)+1:7200)));
        sxsyr(i,j)=sum(wxsyr(i,(7200-period(i,3)+1:7200)));
        sxrxr(i,j)=sum(wxrxr(i,(7200-period(i,3)+1:7200)));
        syrxr(i,j)=sum(wyrxr(i,(7200-period(i,3)+1:7200))); 
        sysxr(i,j)=sum(wysxr(i,(7200-period(i,3)+1:7200)));
        sysyr(i,j)=sum(wysyr(i,(7200-period(i,3)+1:7200)));
        szsxr(i,j)=sum(wzsxr(i,(7200-period(i,3)+1:7200)));
        szsyr(i,j)=sum(wzsyr(i,(7200-period(i,3)+1:7200)));
        sxsxs(i,j)=sum(wxsxs(i,(7200-period(i,3)+1:7200)));
        sysys(i,j)=sum(wysys(i,(7200-period(i,3)+1:7200)));
        szszs(i,j)=sum(wzszs(i,(7200-period(i,3)+1:7200)));
    end
end


Cxsxr=((sxsxr))./(sqrt(sxsxs).*sqrt(sxrxr));
Cxsyr=((sxsyr))./(sqrt(sxsxs).*sqrt(syryr));
Cysxr=((sysxr))./(sqrt(sysys).*sqrt(sxrxr));
Cysyr=((sysyr))./(sqrt(sysys).*sqrt(syryr));
Czsxr=((szsxr))./(sqrt(szszs).*sqrt(sxrxr));
Czsyr=((szsyr))./(sqrt(szszs).*sqrt(syryr));
Cyrxr=((syrxr))./(sqrt(syryr).*sqrt(sxrxr));

C2xsxr=Cxsxr.*conj(Cxsxr);
C2xsyr=Cxsyr.*conj(Cxsyr);
C2ysxr=Cysxr.*conj(Cysxr);
C2ysyr=Cysyr.*conj(Cysyr);
C2zsxr=Czsxr.*conj(Czsxr);
C2zsyr=Czsyr.*conj(Czsyr);
C2yrxr=Cyrxr.*conj(Cyrxr);

% Multivariate coherence
Cxsxryr=((C2xsxr+C2xsyr-2.*real(Cxsxr.*conj(Cxsyr).*conj(Cyrxr)))./(1-(C2yrxr)));
Cysxryr=((C2ysxr+C2ysyr-2.*real(Cysxr.*conj(Cysyr).*conj(Cyrxr)))./(1-(C2yrxr)));
Czsxryr=((C2zsxr+C2zsyr-2.*real(Czsxr.*conj(Czsyr).*conj(Cyrxr)))./(1-(C2yrxr)));
%

clear wxsxr wyryr wxsyr wxrxr wyrxr wxryr wysxr wysyr wzsyr wzsxr wysys wxsxs wzszs;
wxrxs=wxrj.*wxs;
wyryr=wyrj.*wyr;
wxryr=wxrj.*wyr;
wyrxs=wyrj.*wxs;
wxrxr=wxrj.*wxr;
wyrxr=wyrj.*wxr;
wxrys=wxrj.*wys;
wyrys=wyrj.*wys;
wyrzs=wyrj.*wzs;
wxrzs=wxrj.*wzs;
wxsxs=wxsj.*wxs;
wysys=wysj.*wys;
wzszs=wzsj.*wzs;

Qxx=zeros(size(fc,1),period(1,7));
Qxy=zeros(size(fc,1),period(1,7));
Qyx=zeros(size(fc,1),period(1,7));
Qyy=zeros(size(fc,1),period(1,7));
Qzx=zeros(size(fc,1),period(1,7));
Qzy=zeros(size(fc,1),period(1,7));
for i=1:size(wxr,1)
    for j=1:period(i,4)-1
        sxrxs(i,j)=sum(wxrxs(i,1:period(i,3))); 
        syryr(i,j)=sum(wyryr(i,1:period(i,3))); 
        syrxr(i,j)=sum(wyrxr(i,1:period(i,3))); 
        syrxs(i,j)=sum(wyrxs(i,1:period(i,3))); 
        sxrxr(i,j)=sum(wxrxr(i,1:period(i,3))); 
        sxryr(i,j)=sum(wxryr(i,1:period(i,3))); 
        sxrys(i,j)=sum(wxrys(i,1:period(i,3))); 
        syrys(i,j)=sum(wyrys(i,1:period(i,3))); 
        sxrzs(i,j)=sum(wxrzs(i,1:period(i,3))); 
        syrzs(i,j)=sum(wyrzs(i,1:period(i,3))); 
        sxsxs(i,j)=sum(wxsxs(i,1:period(i,3))); 
        sysys(i,j)=sum(wysys(i,1:period(i,3))); 
        szszs(i,j)=sum(wzszs(i,1:period(i,3))); 
    end 
end
for i=1:size(wxr,1)
    for j=period(i,4):period(i,5)
        sxrxs(i,j)=sum(wxrxs(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        syryr(i,j)=sum(wyryr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        syrxr(i,j)=sum(wyrxr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        syrxs(i,j)=sum(wyrxs(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sxrxr(i,j)=sum(wxrxr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sxryr(i,j)=sum(wxryr(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sxrys(i,j)=sum(wxrys(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3)));  
        syrys(i,j)=sum(wyrys(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sxrzs(i,j)=sum(wxrzs(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3)));  
        syrzs(i,j)=sum(wyrzs(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sxsxs(i,j)=sum(wxsxs(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        sysys(i,j)=sum(wysys(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
        szszs(i,j)=sum(wzszs(i,j-floor((period(i,3))/2)+1:j-floor((period(i,3))/2)+period(i,3))); 
    end 
end

for i=1:size(wxr,1)
    for j=period(i,5)+1:7200
        sxrxs(i,j)=sum(wxrxs(i,(7200-period(i,3)+1:7200)));
        syryr(i,j)=sum(wyryr(i,(7200-period(i,3)+1:7200)));
        syrxr(i,j)=sum(wyrxr(i,(7200-period(i,3)+1:7200)));
        syrxs(i,j)=sum(wyrxs(i,(7200-period(i,3)+1:7200)));
        sxrxr(i,j)=sum(wxrxr(i,(7200-period(i,3)+1:7200)));
        sxryr(i,j)=sum(wxryr(i,(7200-period(i,3)+1:7200))); 
        sxrys(i,j)=sum(wxrys(i,(7200-period(i,3)+1:7200)));
        syrys(i,j)=sum(wyrys(i,(7200-period(i,3)+1:7200)));
        sxrzs(i,j)=sum(wxrzs(i,(7200-period(i,3)+1:7200)));
        syrzs(i,j)=sum(wyrzs(i,(7200-period(i,3)+1:7200)));
        sxsxs(i,j)=sum(wxsxs(i,(7200-period(i,3)+1:7200)));
        sysys(i,j)=sum(wysys(i,(7200-period(i,3)+1:7200)));
        szszs(i,j)=sum(wzszs(i,(7200-period(i,3)+1:7200)));
    end
end
Qxx=(sxrxs.*syryr-sxryr.*syrxs)./(sxrxr.*syryr-sxryr.*syrxr);
Qxy=(syrxs.*sxrxr-syrxr.*sxrxs)./(sxrxr.*syryr-sxryr.*syrxr);
Qyx=(sxrys.*syryr-sxryr.*syrys)./(sxrxr.*syryr-sxryr.*syrxr);
Qyy=(syrys.*sxrxr-syrxr.*sxrys)./(sxrxr.*syryr-sxryr.*syrxr);
Qzx=(sxrzs.*syryr-sxryr.*syrzs)./(sxrxr.*syryr-sxryr.*syrxr);
Qzy=(syrzs.*sxrxr-syrxr.*sxrzs)./(sxrxr.*syryr-sxryr.*syrxr);
Txxnew=zeros(size(Cxsxryr,1),size(Cxsxryr,2));
Txynew=zeros(size(Cxsxryr,1),size(Cxsxryr,2));
for i=1:size(Cxsxryr,1)
    k=1;
    for j=1:size(Cxsxryr,2)
        if Cxsxryr(i,j)>=0.95
            Txxnew(i,k)=Qxx(i,j);
            k=k+1;
        end
    end
end
for i=1:size(Txxnew,1)
    TXX(i,1)=sum(Txxnew(i,:))/numel(find(Txxnew(i,:)~=0));
end

for i=1:size(Cxsxryr,1)
    k=1;
    for j=1:size(Cxsxryr,2)
        if Cxsxryr(i,j)>=0.95
            Txynew(i,k)=Qxy(i,j);
            k=k+1;
        end
    end
end
for i=1:size(Txynew,1)
    TXY(i,1)=sum(Txynew(i,:))/numel(find(Txynew(i,:)~=0));
end

Tyxnew=zeros(size(Cysxryr,1),size(Cysxryr,2));
Tyynew=zeros(size(Cysxryr,1),size(Cysxryr,2));
for i=1:size(Cysxryr,1)
    k=1;
    for j=1:size(Cysxryr,2)
        if Cysxryr(i,j)>=0.95
            Tyxnew(i,k)=Qyx(i,j);
            k=k+1;
        end
    end
end
for i=1:size(Tyxnew,1)
    TYX(i,1)=sum(Tyxnew(i,:))/numel(find(Tyxnew(i,:)~=0));
end
for i=1:size(Cysxryr,1)
    k=1;
    for j=1:size(Cysxryr,2)
        if Cysxryr(i,j)>=0.95
            Tyynew(i,k)=Qyy(i,j);
            k=k+1;
        end
    end
end
for i=1:size(Tyynew,1)
    TYY(i,1)=sum(Tyynew(i,:))/numel(find(Tyynew(i,:)~=0));
end

Tzxnew=zeros(size(Czsxryr,1),size(Czsxryr,2));
Tzynew=zeros(size(Czsxryr,1),size(Czsxryr,2));
for i=1:size(Czsxryr,1)
    k=1;
    for j=1:size(Czsxryr,2)
        if Czsxryr(i,j)>=0.95
            Tzxnew(i,k)=Qzx(i,j);
            k=k+1;
        end
    end
end
for i=1:size(Tzxnew,1)
    TZX(i,1)=sum(Tzxnew(i,:))/numel(find(Tzxnew(i,:)~=0));
end
for i=1:size(Czsxryr,1)
    k=1;
    for j=1:size(Czsxryr,2)
        if Czsxryr(i,j)>=0.95
            Tzynew(i,k)=Qzy(i,j);
            k=k+1;
        end
    end
end
for i=1:size(Tzynew,1)
    TZY(i,1)=sum(Tzynew(i,:))/numel(find(Tzynew(i,:)~=0));
end
end

