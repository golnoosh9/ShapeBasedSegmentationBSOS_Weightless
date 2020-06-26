function [fvar]=generate_deg1_mean(fvar,vnum,mon_num,evald,multiplier)


for xx=1:mon_num
     gi=zeros(1,vnum+1);%%%vars plus multiplier
     gi(xx)=1;
     gi(vnum+1)=evald(xx)*multiplier;
     fvar=[fvar;gi];  
    
end



