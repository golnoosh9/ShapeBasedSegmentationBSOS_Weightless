coef_mat=evald*evald';
coef_mat=coef_mat*2;
row1=coef_mat(1,2:5);
row2=coef_mat(2,3:5);
row3=coef_mat(3,4:5);
row4=coef_mat(4,5);
coefs_sc=[row1,row2,row3,row4];