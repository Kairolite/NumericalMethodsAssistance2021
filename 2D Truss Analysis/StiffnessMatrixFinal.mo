function StiffnessMatrixFinal

input Real [:,:,size(kge_mat,2)] kge_mat;
output Real [size(kge_mat,2),size(kge_mat,2)] kgfinal_mat;

algorithm

for row in 1:size(kgfinal_mat,1) loop
  for col in 1:size(kgfinal_mat,2) loop
  kgfinal_mat[row,col] := sum({kge_mat[i,row,col] for i in 1:size(kge_mat,1)});
  end for;
end for;

end StiffnessMatrixFinal;
