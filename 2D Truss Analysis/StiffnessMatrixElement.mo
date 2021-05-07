function StiffnessMatrixElement

input Real [:,7] member_mat;
output Real [size(member_mat,1),4,4] Ke_mat;

protected
Real theta;
Real [3] StiffTrig;
Real [4,4] StiffTrans;
Real [size(member_mat,1)] k_vec;
Real float_error = 10e-10;

algorithm

k_vec := {(member_mat[i,5] * member_mat[i,6] / member_mat[i,7]) for i in 1:size(member_mat,1)};

// Finding stiffness matrix of each element member
for i in 1:size(member_mat,1) loop

  // Clearing the matrices
  StiffTrig := zeros(3);
  StiffTrans := zeros(4,4);
  
  // Converting degrees to radians
  theta := Modelica.SIunits.Conversions.from_deg(member_mat[i,4]);

  // {cos^2, sin^2, sincos}
  StiffTrig := {(Modelica.Math.cos(theta))^2,
                (Modelica.Math.sin(theta))^2,
                (Modelica.Math.sin(theta)*Modelica.Math.cos(theta))};
  
  // Handle float error elements in StiffTrig
  for t in 1:size(StiffTrig,1) loop
    if abs(StiffTrig[t]) <= float_error then
      StiffTrig[t] := 0;
    end if;
end for;
  
  // Construct stiffness transformation matrix
  StiffTrans := [  StiffTrig[1],    StiffTrig[3], -1*StiffTrig[1], -1*StiffTrig[3];
                   StiffTrig[3],    StiffTrig[2], -1*StiffTrig[3], -1*StiffTrig[2];
                -1*StiffTrig[1], -1*StiffTrig[3],    StiffTrig[1],    StiffTrig[3];
                -1*StiffTrig[3], -1*StiffTrig[2],    StiffTrig[3],    StiffTrig[2]];
  
  // Multiply in stiffness constant of element, add final stiffness matrix to Ke_mat
  for m in 1:4 loop
    for n in 1:4 loop
      Ke_mat[i,m,n] := k_vec[i] * StiffTrans[m,n];
    end for;
  end for;

end for;

end StiffnessMatrixElement;
