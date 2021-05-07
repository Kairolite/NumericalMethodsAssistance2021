function GaussJordan

// Gauss-Jordan Algorithm
// Transforms input matrix A into reduced row echelon form matrix B
// Christopher S.E. November 2020

input Real [:,:] A;                     // An augmented matrix of m*n
output Real [size(A,1),size(A,2)] B;    // Output matrix in reduced row echelon form

// Local variables
protected
Integer h = 1;          // Initialize pivot row
Integer k = 1;          // Initialize pivot column
Integer m = size(A,1);  // Number of rows in matrix
Integer n = size(A,2);  // Number of columns in matrix
Integer c = 0;          // Index counter
Integer max_row;        // Row index of max number in pivot column

Real [:] pivot_column;  // Vector containing pivot column data
Real [:] pivot_row;     // Vector containing backwards pivot row data rows
Real r;                 // Ratio value used for row operations

// Limit to handle floating point errors
Real float_error = 10e-6;

algorithm

// Transfer input matrix A into variable B
B := A;

while h <= m and k <= n loop
  
  // Dealing with floating point errors
  for i in 1:m loop
    for j in 1:n loop
      if abs(B[i,j]) <= float_error then
        B[i,j] := 0;
      end if;
    end for;
  end for;

  // Finding the pivot
    pivot_column := {B[i,h] for i in h:m};
  
    // Get position index of lowest row with greatest pivot number
    c:= h-1;
    for element in pivot_column loop
      c := c+1;
      if abs(element) == max(abs(pivot_column)) then
        max_row := c;
      end if;
    end for;
  
  // No pivot in this column, move on to next column
  if B[max_row, k] == 0 then
    k := k+1;
  
  else
    // Swap Rows h <-> max_row
    if sum({(B[h,i] - B[max_row,i]) for i in 1:n}) <> 0 then
      for i in 1:n loop
        B[h,i] := B[h,i] - B[max_row,i];
        B[max_row,i] := B[h,i] + B[max_row,i];
        B[h,i] := B[max_row,i] - B[h,i];
      end for;
    end if;
    
    // Divide pivot row by pivot number
    r := B[h,k];
    for i in 1:n loop
      B[h,i] := B[h,i] / r;
    end for;
    
    // For all rows below the pivot
    for i in (h+1):m loop
      
      // Store the ratio of the row to the pivot
      r := B[i,k] / B[h,k];
      
      // Set lower part of pivot column to zero
      B[i,k] := 0;
      
      // Operations on the remaining row elements
      for j in (k+1):n loop
          B[i,j] := B[i,j] - B[h,j] * r;
      end for;
      
    end for;
    
    // Move on to next pivot row and column
    h := h+1;
    k := k+1;
    
  end if;
  
end while;

// The matrix is now in row echelon form

// Set values of (h,k) to (m,n)
h := m;
k := n;

while h >= 1 and k >=1 loop

  // Dealing with floating point errors
  for i in 1:m loop
    for j in 1:n loop
      if abs(B[i,j]) <= float_error then
        B[i,j] := 0;
      end if;
    end for;
  end for;

  // Finding the pivot
    pivot_row := {B[h,i] for i in 1:k};
    
    // Get position index k of pivot
    c := 0;
    for element in pivot_row loop
      c := c+1;
      if element <> 0 then
        break;
      end if;
    end for;
    k := c;

  // No pivot in this row, move on to next row
  if B[h, k] == 0 then
    h := h-1;
    
  else
    
    // Perform row operations
    for i in 1:(h-1) loop
      r := B[i,k];
      for g in 1:n loop
      B[i,g] := B[i,g] - B[h,g] * r;
      end for;
    end for;
    
    // Move on to next pivot row and column
    h := h-1;
    k := k-1;
    
  end if;

end while;

// The matrix is now in reduced row echelon form

end GaussJordan;
