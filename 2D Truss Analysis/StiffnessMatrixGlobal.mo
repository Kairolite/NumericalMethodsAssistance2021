function StiffnessMatrixGlobal

input Real member_mat [:,7];
input Real node_mat [:,5];
output Real kge_mat [size(member_mat,1),2*size(node_mat,1),2*size(node_mat,1)];

protected
Real [size(member_mat,1),4,4] Ke_mat;
Integer [2] nodes;
Integer ke_row = 1;
Integer ke_col = 1;
Integer g_row;
Integer g_col;

algorithm

// initialize element stiffness matrix and its global format
Ke_mat := StiffnessMatrixElement(member_mat);
kge_mat := zeros(size(member_mat,1),2*size(node_mat,1),2*size(node_mat,1));

// for each member element stiffness matrix, find the corresponding global matrix coordinates
for member in 1:size(member_mat,1) loop
  nodes:= integer({member_mat[member,node] for node in 2:3}); // nodes {i,j}
  
  ke_row := 1;
  while ke_row <= size(Ke_mat,2) loop
    g_row := (2 * nodes[div(ke_row + rem(ke_row,2),2)]) - rem(ke_row,2);
    
    ke_col := 1;
    while ke_col <= size(Ke_mat,3) loop
      g_col := (2 * nodes[div(ke_col + rem(ke_col,2),2)]) - rem(ke_col,2);
      
      // Transfers values from the element format into its position in the global format
      kge_mat[member, g_row, g_col] := Ke_mat[member, ke_row, ke_col];
      
      ke_col := ke_col + 1;
    end while;
    
    ke_row := ke_row + 1;
  end while;
  
end for;

end StiffnessMatrixGlobal;
