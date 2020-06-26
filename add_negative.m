function []=add_negative(begin,pos,POP,mon_num)
%neg_add=coeffs;
add=1;
load('coeffs.mat');
%begin=0;
if begin==1
    sequence=' ';
    pos_new=[];
    neg_new=[];
    pos_count=1;
    ex_count=1;
    
else
    %load('rand_four_example.mat');
    load('rand_air_example.mat');
end

if add==1
    if pos==1
        pos_new(:,pos_count)=coeffs;
        pos_count=pos_count+1;
        sequence=strcat(sequence,'p');
    else
        sequence=strcat(sequence,'n');
        neg_new(:,ex_count)=coeffs;
        ex_count=ex_count+1;
    end
end

% save('rand_four_example.mat','pos_new','neg_new','ex_count','pos_count');
save('rand_air_example.mat','pos_new','neg_new','ex_count','pos_count','sequence');
pnum=size(pos_new);
nnum=size(neg_new);
% neg_new(:,60:70)=[];
if(pos==0 | pnum+nnum==1 )
    
    inv_coeffs=learn_inv(pnum(2),nnum(2),pos_new,neg_new);
    %save('four_random_inv','inv_coeffs');
    save('air_random_inv','inv_coeffs');
end
if(pos==1)
;
end
