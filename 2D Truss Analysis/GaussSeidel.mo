function GaussSeidel

// Gauss-Seidel Algorithm
// Christopher S.E. May 2021

input Real [:,:] A;             // An augmented matrix of m*n
output Real [size(A,1)] X;      // Output solution vector
// output Real [size(A,1)] error;  // Approximate error

protected
Integer m = size(A,1);  // Number of rows in augmented matrix
Integer n = size(A,2);  // Number of columns in augmented matrix
Integer iter = 1000;    // Number of iterations

// For aX = b
Real a [size(A,1),size(A,1)]; // Square matrix (m x m)
Real b [size(A,1)];           // Column vector (m x 1)


Real prev_x;  // Value of x from previous iteration
Real row_sol; // Part of equation numerator when solving for X[i]

algorithm

// Transfer from input augmented matrix form into aX = b form
for i in 1:m loop
  for j in 1:m loop
    a[i,j] := A[i,j];
  end for;
  b[i] := A[i,end];
end for;

// Set initial solution vector and error vector to zero
X := {0 for x in 1:m};
// error := {0 for e in 1:m};

// Gauss Seidel Iterative Algorithm
// X[i] = ( b[i] - [sum(A[i,j]*X[j])] - A[i,i] ) / A[i,i]
// for calculation purposes: row_sol = [sum(A[i,j]*X[j])] - A[i,i]
// therefore: X[i] = ( b[i] - row_sol ) / A[i,i]

while iter > 0 loop
  
  // For each X[i] value
  for i in 1:m loop
    
    // Reset row_sol and store previous x value
    row_sol := 0;
    prev_x := X[i];
    
    // Solving for unknown value X[i]
    for j in 1:m loop
      row_sol := row_sol + (a[i,j] * X[j]);
    end for;
    row_sol := row_sol - (a[i,i] * X[i]);
    
    X[i] := (b[i] - row_sol)/ a[i,i];
    
    // Calculate the error for X[i]
/*    
    if X[i] <> 0 then
      error[i] := abs((X[i] - prev_x)/X[i]) * 100;
    else
      error[i] := 0; // avoid div by zero
    end if;
*/    
  end for;
  
  iter := iter - 1; // Move on to the next iteration
  
end while;

end GaussSeidel;
