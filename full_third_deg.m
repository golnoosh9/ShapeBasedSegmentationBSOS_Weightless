varst37=[];



for xx=1:mon_num
    varst2=zeros(1,mon_num+cnum+2+vercount+cnum+3);
    varst2(xx)=3;
    varst37=[varst37;varst2];
    
end

for ii=1:mon_num
    for jj=1:mon_num
        if ii==jj
            continue;
        end
        varst2=zeros(1,mon_num+cnum+2+vercount+cnum+3);
        varst2(ii)=2;
        varst2(jj)=1;
       % varst2(mm)=1;
        varst37=[varst37;varst2];
     end

end


for ii=1:mon_num-2
    for jj=ii+1:mon_num-1
        for mm=jj+1:mon_num
        varst2=zeros(1,mon_num+cnum+2+vercount+cnum+3);
        varst2(ii)=1;
        varst2(jj)=1;
        varst2(mm)=1;
        varst37=[varst37;varst2];
        end
    end
end
