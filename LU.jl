# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%% Stable LU Decomposition %%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#   % For numerical stability you want to swap rows with
#   % largest pivot. Note there are other ways to make this 
#   % algorithm more efficient

# % Example [eps 1; 1 1]x = [1; 2] 

# format rat;
s = [1 2 3]
# % Gaussian Elimination (LU decomposition)
A = [1  2  3; 
     2  1 -1; 
     4  0  1];

# % Preform Gaussian elimination step on Identity matrix
# % I_3 (swap R3 <-> R1)
s1 = [3 2 1]
P1 = [0 0 1;
      0 1 0;
      1 0 0];
  
P1*A

# %      4     0     1
# %      2     1    -1
# %      1     2     3

# % I_3 (-2/4 R1 + R2)
E1 = [1    0  0;
     -2/4  1  0;
      0    0  1];
  
E1*P1*A
 
# %        4              0              1       
# %        0              1             -3/2     
# %        1              2              3       

# % I_3 (-1/4 R1 + R3)
E2 = [  1   0  0;
        0   1  0;
      -1/4  0  1];
  
A1 = E2*E1*P1*A

L1 = inv(E2*E1)

# %        4              0              1       
# %        0              1             -3/2     
# %        0              2             11/4   

# % Swap for largest pivot
# % I_3 (R3 <-> R2)
s2 = [3 1 2]
P2 =  [1 0 0;
       0 0 1;
       0 1 0];
  
L1 = P2*L1
P2*E2*E1*P1*A

# %        4              0              1       
# %        0              2             11/4     
# %        0              1             -3/2    

# % I_3 (-1/2 R2 + R3)
E3 = [  1   0   0;
        0   1   0;
        0  -1/2 1];
  
U = E3*P2*E2*E1*P1*A

# %        4              0              1       
# %        0              2             11/4     
# %        0              0            -23/8     

# L2 = 

L = inv(E3*P2*E2*E1*P1) #% not lower triangular?

# %        1/4            1              0       
# %        1/2            1/2            1       
# %        1              0              0      

# % The issue is that the swaps E1 and E4 are not lower triangular, however
# % its still easy to solve Ax = b

b = [23;0;0];
x = A\b # % solution

# %       -1       
# %        6       
# %        4     

U = E3*P2*E2*E1*P1*A;

b_hat = E3*P2*E2*E1*P1*b;

# % Augemnted system [U | b_hat] can be solved with backsubstitution
# horzcat(U, b_hat)

# %        4              0              1              0       
# %        0              2             11/4           23       
# %        0              0            -23/8          -23/2    

x3 = (-8/23 * -23/2) # % = 4
# % => x2 = 6
# % => x1 = -1

x = U\b_hat

# %       -1       
# %        6       
# %        4     

# %%%%% Alternatively

# % Permutation Matrix (The swaps we performed
P = P2*P1

# % P =
# % 
# %        0              0              1       
# %        1              0              0       
# %        0              1              0     

# % Upper triangular
U = E3*P2*E2*E1*P1*A

# %        4              0              1       
# %        0              2             11/4     
# %        0              0            -23/8     

# % But L is NOT lower triangular
L1 = inv(E3*P2*E2*E1*P1) 

# % L =
# % 
# %        1/4            1              0       
# %        1/2            1/2            1       
# %        1              0              0    

# % IS lower triangular
L = P*L1
# uniqueidx(v) = unique(i -> v[i], eachindex(v))


# % L =
# % 
# %        1              0              0       
# %        1/4            1              0       
# %        1/2            1/2            1     


# % Yeilds permuted decomposition
# % L*U = 
# % 
# %        4              0              1       
# %        1              2              3       
# %        2              1             -1       


# % A =
# % 
# %        1              2              3       
# %        2              1             -1       
# %        4              0              1    

# % Since P'P = I, we can undo permutation.
# % P'*L*U = P'P*L1*U = L1*U = A =
# % 
# %        1              2              3       
# %        2              1             -1       
# %        4              0              1    


# % Check this is what matlab produces
# [l,u,p] = lu(A)

# % l =
# % 
# %        1              0              0       
# %        1/4            1              0       
# %        1/2            1/2            1       
# % 
# % 
# % u =
# % 
# %        4              0              1       
# %        0              2             11/4     
# %        0              0            -23/8     
# % 
# % 
# % p =
# % 
# %        0              0              1       
# %        1              0              0       
# %        0              1              0       

s = [1;2;3;4]
A = [2 1 1 0;
     4 3 3 1;
     8 7 9 5;
     6 7 9 8]

# swap rows 3 and 1
s1 = [3;2;1;4]
A = A[s1,:]

E1 = [1 0 0 0; 
      -4/8 1 0 0;
      -2/8 0 1 0;
      -6/8 0 0 1]

A1 = E1*A
L1 = inv(E1)

# swap rows 2 and 4
s2 = [3;4;1;2]
# A1[:,2:end] = A1[s2,2:end]

E2 = [1 0 0 0; 
      0 1 0 0;
      0 0 1 0;
      0 0 0 1]

A2 = E2*A1