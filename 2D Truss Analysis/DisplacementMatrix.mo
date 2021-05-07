function DisplacementMatrix

input Real node_mat [:,5];
input Real kgfinal_mat [2*size(node_mat,1),2*size(node_mat,1)];
output Real displacement_mat [2*size(node_mat,1)];

protected
Real F [2*size(node_mat,1)];
Real aug_mat [2*size(node_mat,1),2*size(node_mat,1)+1];

// [K]{U} = {F}
algorithm
F := {node_mat[integer((n+rem(n,2))/2),integer(3-rem(n,2))] for n in 1:size(F,1)};

// Construct Augmented Matrix of [kgfinal_mat|F]
for row in 1:size(aug_mat,1) loop
  for col_kg in 1:size(kgfinal_mat,1) loop
    aug_mat[row,col_kg] := kgfinal_mat[row,col_kg];
  end for;
  
  aug_mat[row,end] := F[row];
  
end for;

// Changes rows according to boundary conditions
for row in 1:size(aug_mat,1) loop
  if node_mat[integer((row+rem(row,2))/2),integer(5-rem(row,2))] == 1 then
    for col in 1:size(aug_mat,2)-1 loop
      aug_mat[row,col] := 0;
    end for;
    aug_mat[row,row] := 1;
  end if;
end for;

// Solve into Reduced Row Echelon Form by Gauss Jordan function
aug_mat := GaussJordan(aug_mat);

displacement_mat := {aug_mat[i,end] for i in 1:size(aug_mat,1)};

end DisplacementMatrix;
