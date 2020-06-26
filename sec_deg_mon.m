varst7=[];



for xx=1:mon_num
    varst2=zeros(1,mon_num+cnum+2);
    varst2(xx)=2;
    varst7=[varst7;varst2];    
    
end

for ii=1:mon_num
    for jj=ii+1:mon_num
        varst2=zeros(1,mon_num+cnum+2);
        varst2(ii)=1;
        varst2(jj)=1;
        varst7=[varst7;varst2];
    end
end
