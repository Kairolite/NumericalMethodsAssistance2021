function TrussPostProcess

input Real member_mat [:,7];
input Real node_mat [:,5];
input Real kgfinal_mat [2*size(node_mat,1),2*size(node_mat,1)];
input Real displacement_mat [2*size(node_mat,1)];

// Reaction Forces Vector
output Real reaction_mat [2*size(node_mat,1)];

// Internal Forces and Normal Stresses Matrix
output Real uf_mat [size(member_mat,1),7];

protected
Real F [2*size(node_mat,1)];
Real k_vec [size(member_mat,1)];
Real theta;
Real trig [2];
Real T [4,4];
Real mem_nodes [2];
Real U_mat [4];
Real local_u [4];

algorithm
// Load Matrix
F := {node_mat[integer((n+rem(n,2))/2),integer(3-rem(n,2))] for n in 1:size(F,1)};

// Stiffness Constant Vector
k_vec := {(member_mat[i,5] * member_mat[i,6] / member_mat[i,7]) for i in 1:size(member_mat,1)};

// Reaction Forces
reaction_mat := kgfinal_mat * displacement_mat - F;

// Cleaning float errors
for force in 1:size(reaction_mat,1) loop
  if abs(reaction_mat[force]) <= 10e-6 then
    reaction_mat[force] := 0;
  end if;
end for;

// Internal Forces and Normal Stresses
for row in 1:size(uf_mat,1) loop
  uf_mat[row,1] := member_mat[row,1]; // Transfer member #
  
  // Construct transformation matrix
  theta := Modelica.SIunits.Conversions.from_deg(member_mat[row,4]);
  trig := {Modelica.Math.sin(theta), Modelica.Math.cos(theta)};
  T := [trig[2], trig[1], 0, 0; -trig[1], trig[2], 0, 0;
        0, 0, trig[2], trig[1]; 0, 0, -trig[1], trig[2]];
  
  // Retrieve U displacement values
  mem_nodes := {member_mat[row,n] for n in 2:3}; // Nodes {i,j}
  U_mat := {displacement_mat[integer(2*mem_nodes[integer((k+rem(k,2))/2)] - rem(k,2))] for k in 1:4};
  
  // Calculate internal forces u
  local_u := T * U_mat;
  
  // Transfer values to uf matrix
  for u_col in 2:5 loop
    uf_mat[row,u_col] := local_u[u_col-1];
  end for;
  
  // Calculate internal forces and normal stresses
  uf_mat[row,6] := k_vec[row] * (uf_mat[row,2] - uf_mat[row,4]);
  uf_mat[row,7] := uf_mat[row,6] / member_mat[row,5];
  
end for;

// Cleaning float errors
for value_col in 1:size(uf_mat,1) loop
  for value_row in 2:size(uf_mat,1) loop
    if abs(uf_mat[value_row,value_col]) <= 10e-6 then
      uf_mat[value_row,value_col] := 0;
    end if;
  end for;
end for;

end TrussPostProcess;
