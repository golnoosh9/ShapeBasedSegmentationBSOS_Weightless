varst7=[];
for xx=1:mon_num
    varst2=zeros(1,mon_num+cnum+2);
    varst2(xx)=1;
    varst2(mon_num+j)=1; %for weight
    varst7=[varst7;varst2];    
    
end

