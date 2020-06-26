function [fvar]=generate_deg2_monomials(fvar,vnum,mon_num,evald1,evald2,multiplier)


for i=1:mon_num
    for j=1:mon_num
     gi=zeros(1,vnum+1);%%%vars plus multiplier
     gi(i)=1;
     gi(j)=1;
     gi(vnum+1)=evald1(i)*evald2(j)*multiplier;
     fvar=[fvar;gi];  
    
end
end