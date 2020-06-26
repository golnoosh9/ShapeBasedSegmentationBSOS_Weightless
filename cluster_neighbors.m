function [xp,yp]= cluster_neighbors(bs,cnum,csize)
    sbs=size(bs);
    count=1;
    for i=1:sbs(1)%%%% for each neighborhood cluster into csize untill finished
        i
        points=bs{i};
        spoints=size(points);
        pcount=1;
        for j=1:ceil(spoints(1)/csize)
            pcount
            if pcount+csize<=spoints(1)
                 xp{count}=points(pcount:pcount+csize-1,1);
                 yp{count}=points(pcount:pcount+csize-1,2);
            else
                pcount
                spoints(1)
                
               xp{count}=points(pcount:spoints(1),1);
                yp{count}=points(pcount:spoints(1),2);
            end
            count=count+1;
             pcount=pcount+csize;
        end
        
        
        
        
        
    end
        
        
        
        
        






end