function [fvar]=generate_deg2_monomials(fvar,vnum,mon_num,evald,maximum,multiplier,const)


for xx=1:mon_num
     gi=zeros(1,vnum+1);%%%vars plus multiplier
     gi(xx)=2;
     gi(vnum+1)=evald(xx)^2*multiplier;
     fvar=[fvar;gi];  
    
end



for ii=1:mon_num
    for jj=ii+1:mon_num
         gi=zeros(1,vnum+1);%%%vars plus multiplier
     gi(ii)=1;
     gi(jj)=1;
     gi(vnum+1)=2*evald(ii)*evald(jj)*multiplier;
     fvar=[fvar;gi];  
    end
end
if(const==1)
gi=zeros(1,vnum+1);
gi(vnum+1)=-(2/maximum)^2;
fvar=[fvar;gi];
end
