class TrussExample

// PROCEDURE
/*
Discretize Nodes and Elements
Table = [element, node i, node j, theta, area, modulus, length]

Assume Element Stiffness Constants
k = AE/L; for each group of elements

Find Element Stiffness Matrix
[K](e) = k[Tstiffness]; for each element

Assemble Global Stiffness Matrix
[K](e) -> [K](eG); Convert into global format
[K](G) = sum(K(eG))

Boundary Conditions and Loads
[K](G) * {U} = {F}; insert external loads to F
optional: eliminate any empty rows where U = 0

Solve Global Displacements
[K](G) * {U} = {F}; solve for {U} using solve()

Reaction Forces
{R} = [K](G) * {U} - {F}; simply calculate {R}

Internal Forces
{U} = [T]{u}; find {u}
{fix, fjx} = k * {(uix-ujx),(ujx-uix)}; for each member

Normal Stresses
sigma = fix/A OR sigma = (E/L)*(uix-ujx); for each member
*/

// DECLARATIONS
// Data for each member: [element #, node i, node j, theta, area, modulus, length]
// SI Units: [integer, integer, integer, degrees, meters^2, pascals, meters]
// Imperial Units: [integer, integer, integer, degrees, inches^2, lb/inches, inches]
parameter Real [:,7] member = [1, 1, 2,   0, 8, 1.9e6, 36;
                               2, 2, 3, 135, 8, 1.9e6, 50.9;
                               3, 3, 4,   0, 8, 1.9e6, 36;
                               4, 2, 4,  90, 8, 1.9e6, 36;
                               5, 2, 5,  45, 8, 1.9e6, 50.9;
                               6, 4, 5,   0, 8, 1.9e6, 36];

// External loads and boundary conditions for each node: [node #, FX, FY, UX(1:Fixed), UY(1:Fixed)]
parameter Real [:,5] node_load = [1, 0, 0   , 1, 1;
                                  2, 0, 0   , 0, 0;
                                  3, 0, 0   , 1, 1;
                                  4, 0, -500, 0, 0;
                                  5, 0, -500, 0, 0];

// Vector for equivalent stiffness constant of each member
Real k [size(member,1)];

// Arrays for Stiffness Matrices
Real Ke [size(member,1), 4, 4]; // Element
Real KGe [size(member,1), 2*size(node_load,1), 2*size(node_load,1)]; // Element in Global Form
Real KGt [2*size(node_load,1), 2*size(node_load,1)]; // Total Global

// Displacement Matrix
Real U [2*size(node_load,1)];

// Postprocessing Phase
Real R [2*size(node_load,1)]; // Reaction forces for each node
Real InForce_NormStress [size(member,1),7]; // Internal forces and normal stress for each member
// Format: [member #, uix, uiy, ujx, ujy, internal force, normal stress]

// EQUATIONS
equation

// Calculate stiffness constants for each member element
k = {(member[i,5] * member[i,6] / member[i,7]) for i in 1:size(member,1)};

// Find stiffness matrix for each element
Ke = StiffnessMatrixElement(member);

// Assemble the global stiffness matrix
KGe = StiffnessMatrixGlobal(member, node_load);
KGt = StiffnessMatrixFinal(KGe);

// Solve for the displacement matrix
U = DisplacementMatrix(node_load, KGt);

// Postprocessing Phase
(R, InForce_NormStress) = TrussPostProcess(member, node_load, KGt, U);

end TrussExample;
