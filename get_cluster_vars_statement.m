init_w=1;
for i=1:size_bs(1)
    tsize=size(cell2mat(bs(i)));
    cluster_num=tsize(1)/cluster_size;
    resid_cluster=mod(tsize(1),cluster_size);
    for j=1:cluster_num
      end_w=end_w+cluster_size;
     cluster_len=end_w-init_w+1;
%fvar=[];
        fmean=[];
        for xx=init_w:end_w
            xs=xp{xx};
             ys=yp{xx};
    
             xs=image_multiplier*(xs-meanx)/maximum;
              ys=image_multiplier*(ys-meany)/maximum;
  
                eval=(subs(VEC,[x,y],[xs,ys]));
              evald=double(eval);
            fmean=Add_statements_BSOS(fmean, calculate_taylor_expansion(taylor_negative_order,vnum,mon_num,evald/cnum,1/maximum));

          count=count+1;
         tcount=tcount+1;
        end
        fvariance_cluster=Add_statements_BSOS(fvariance_cluster,Multiply_Statements_BSOS(fmean,fmean));
         for xx=init_w:end_w
            xs=xp{xx};
             ys=yp{xx};
    
             xs=image_multiplier*(xs-meanx)/maximum;
              ys=image_multiplier*(ys-meany)/maximum;
  
                eval=(subs(VEC,[x,y],[xs,ys]));
              evald=double(eval);
              fstate=calculate_taylor_expansion(taylor_negative_order,vnum,mon_num,evald,1/maximum);
              
            fvariance_cluster=Add_Statements_BSOS(fvariance_cluster,Multiply_Statements_BSOS(fstate,fstate));
              fstate=calculate_taylor_expansion(taylor_negative_order,vnum,mon_num,-2*evald,1/maximum);
            fvariance_cluster=Add_Statements_BSOS(fvariance_cluster,Multiply_Statements_BSOS(fstate,fmean));

          count=count+1;
         tcount=tcount+1;
        end
      init_w=init_w+cluster_size;
     
    end
    %%%%%%residue
    if resid_cluster>0
    end_w=end_w+resid_cluster;
      for xx=init_w:end_w
            xs=xp{xx};
             ys=yp{xx};
    
             xs=image_multiplier*(xs-meanx)/maximum;
              ys=image_multiplier*(ys-meany)/maximum;
  
                eval=(subs(VEC,[x,y],[xs,ys]));
              evald=double(eval);
            fmean=Add_statements_BSOS(fmean, calculate_taylor_expansion(taylor_negative_order,vnum,mon_num,evald/cnum,1/maximum));

          count=count+1;
         tcount=tcount+1;
      end
        fvariance_cluster=Add_statements_BSOS(fvariance_cluster,Multiply_Statements_BSOS(fmean,fmean));  
      for xx=init_w:end_w
            xs=xp{xx};
             ys=yp{xx};
    
             xs=image_multiplier*(xs-meanx)/maximum;
              ys=image_multiplier*(ys-meany)/maximum;
  
                eval=(subs(VEC,[x,y],[xs,ys]));
              evald=double(eval);
               fstate=calculate_taylor_expansion(taylor_negative_order,vnum,mon_num,evald,1/maximum);
              
            fvariance_cluster=Add_Statements_BSOS(fvariance_cluster,Multiply_Statements_BSOS(fstate,fstate));
              fstate=calculate_taylor_expansion(taylor_negative_order,vnum,mon_num,-2*evald,1/maximum);
            fvariance_cluster=Add_Statements_BSOS(fvariance_cluster,Multiply_Statements_BSOS(fstate,fmean));

          count=count+1;
         tcount=tcount+1;
     end
    init_w=init_w+resid_cluster;
    end
end
