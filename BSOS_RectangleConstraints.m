%%%0.01<x<5 => 0<x-0.01<4.99 => 0<(x-0.01)/4.99<1
for i=1:maxy-miny%%%%2-eval is bounded 
    eval=(subs(VEC,[x,y],[(minx-meanx)/(maximum),sminx(i)]));
    evald=double(eval);
    
    gtotal=[];
        for mon=1:mon_num
            gi=zeros(1,vnum+1);%%%vars plus multiplier
      
           gi(mon)=1; 
           gi(vnum+1)=evald(mon)/4.99;
           gtotal=[gtotal;gi];
        end
         gi=zeros(1,vnum+1);
         gi(vnum+1)=-0.01/4.99;
pop.G{ineq_count}=[gtotal;gi];

        ineq_count=ineq_count+1;
    
    eval=(subs(VEC,[x,y],[(maxx-meanx)/(maximum),smaxx(i)]));
    
    evald=double(eval);
    
     gtotal=[];
        for mon=1:mon_num
            gi=zeros(1,vnum+1);%%%vars plus multiplier

           gi(mon)=1; 
           gi(vnum+1)=evald(mon)/4.99;
           gtotal=[gtotal;gi];
        end
         gi=zeros(1,vnum+1);
         gi(vnum+1)=-0.01/4.99;
pop.G{ineq_count}=[gtotal;gi];

        ineq_count=ineq_count+1;
end
return

for i=1:maxx-minx%%%%2-eval is bounded
    eval=(subs(VEC,[x,y],[sminy(i),(miny-meany)/(maximum)]));
    evald=double(eval);
    
    gtotal=[];
        for mon=1:mon_num
            gi=zeros(1,vnum+1);%%%vars plus multiplier
     
           gi(mon)=1; 
           gi(vnum+1)=evald(mon)/4.99;
           gtotal=[gtotal;gi];
        end
         gi=zeros(1,vnum+1);
         gi(vnum+1)=-0.01/4.99;
        pop.G{ineq_count}=[gtotal];
        
        ineq_count=ineq_count+1;
    
    eval=(subs(VEC,[x,y],[smaxy(i),(maxy-meany)/(maximum)]));
    
    evald=double(eval);
    
     gtotal=[];
        for mon=1:mon_num
            gi=zeros(1,vnum+1);%%%vars plus multiplier

           gi(mon)=1; 
           gi(vnum+1)=evald(mon)/4.99;
           gtotal=[gtotal;gi];
        end
         gi=zeros(1,vnum+1);
         gi(vnum+1)=-0.01/4.99;
        pop.G{ineq_count}=[gtotal];
        
        ineq_count=ineq_count+1;
end
