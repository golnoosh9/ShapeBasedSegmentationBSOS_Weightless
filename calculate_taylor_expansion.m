function [tdegrees]=calculate_taylor_expansion(order,vnum,mon_num,coeff_constant_alpha,e)
% order=4;
% coeff_constant_alpha=neighbor_evald;
% e=1/maximum;

poly_coeffs = sym('c',[mon_num 1]);

raised_degree=1:order;
VEC=monomials_gen(poly_coeffs,raised_degree);

taylor_exp=zeros(vnum+1,1);
sVEC=size(VEC);
sVEC=sVEC(1);

%original_distance_function=(e^2)/((coeff_constant_alpha'*poly_coeffs)^2);
original_distance_function=1/exp(((coeff_constant_alpha'*poly_coeffs)/e)^2);

taylor_exp=0;
vars=symvar(original_distance_function);
taylor_exp = feval(symengine,'mtaylor',original_distance_function,'[c1=0,c2=0,c3=0,c4=0,c5=0,c6=0]',order+1);

 

[total_coeffs, terms]=coeffs(taylor_exp);
tdegrees=generate_monomial_basis_degree(order,vnum);
sterms=size(terms);
sterms=sterms(2);

for i=1:sVEC
    for j=1:sterms
   if(VEC(i)==(terms(j)))

       tdegrees(i,vnum+1)=total_coeffs(j)*(1);
       break;
   end
end
end
end
