function[final_state]= Multiply_Statements_BSOS(state1,state2,vnum)
final_state=[];
size1=size(state1);
size2=size(state2);
size1=size1(1);
size2=size2(1);
for i=1:size1
    for j=1:size2
       exp1=state1(i,:);
       exp2=state2(j,:);
       newexp(1:vnum)=exp1(1:vnum)+exp2(1:vnum);
       newexp(vnum+1)=exp1(vnum+1)*exp2(vnum+1);
       final_state=[final_state;newexp];
    end
end
