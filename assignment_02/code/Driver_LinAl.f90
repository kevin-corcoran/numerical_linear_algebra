Program Driver_LinAl

  ! use LinAl, only: mat, msize, nsize, readMat
  use LinAl, only: mat, msize, nsize, readMat, traceMat, twoNorm, printMat, &
    printColumnNorm, partGaussElim, print_mat, LUdecomp, LU_partial, print_vec, &
    backSubstitutionLU

  implicit none
  
  character(len=100) :: myFileName
  integer :: i,j
  real :: traceA
  real :: p_i, e_

  real :: norm
  integer, dimension(2) :: sizeA, sizeB
  real, dimension(:,:), allocatable :: A, B, L, U, X, P_, b_, x_
  integer, dimension(:), allocatable :: s
  logical :: SINGLR
  SINGLR = .FALSE.


  p_i = 3.141592653589793
  e_ = 2.718281828459045

  myFileName = 'Amat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  allocate(A(msize,nsize))
  allocate(L(msize,nsize))
  allocate(U(msize,nsize))
  ! Always initialize with zeros
  mat = 0.0
  L = 0.0
  U = 0.0
  call readMat(myFileName)
  A = mat
  sizeA = [msize, nsize]

  write(*,*) "A"
  ! call printMat(A, [msize, nsize])
  call print_mat(A)
  call traceMat(A, msize, traceA)
  print *, ""
  write(*,*) "Trace of A"
  write(*,*) (traceA)
  call printColumnNorm(A)
  deallocate(mat)


  myFileName = 'Bmat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  allocate(B(msize,nsize))
  ! Always initialize with zeros
  mat = 0.0
  call readMat(myFileName)
  B = mat
  sizeB = [msize, nsize]

  print *, ""
  write(*,*) "B"
  call print_mat(B)
  !call printMat(B, sizeB)

  deallocate(mat)

  ! Gaussian elimination
  call partGaussElim(A, B, sizeA, sizeB, SINGLR)







  ! LU decomposition 
  deallocate(A)
  myFileName = 'Amat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  allocate(A(msize,nsize))
  mat = 0.0
  call readMat(myFileName)
  A = mat
  deallocate(mat)

  allocate(s(sizeA(1)))
  SINGLR = .FALSE.

  print *, ""
  print *, ""
  print *, ""
  print *, "A before LU"
  call print_mat(A)
  call LU_partial(A, sizeA(1), s, SINGLR)
  print *, "A after LU"
  call print_mat(A)
  print *, "Permutation vector"
  call print_vec(s)








  deallocate(B)
  myFileName = 'Bmat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  allocate(B(msize,nsize))
  allocate(X(msize,nsize))
  ! Always initialize with zeros
  mat = 0.0
  X = 0.0
  call readMat(myFileName)
  B = mat
  sizeB = [msize, nsize]

  deallocate(mat)

  call backSubstitutionLU(A, msize, B, s, X)



  !P_ = [1, 2, 3, 1; -3, 2, 5, 1]
  myFileName = 'P.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize))
  allocate(P_(msize,nsize))
  allocate(b_(msize,1))
  allocate(x_(msize,1))
  mat = 0.0
  P_ = 0.0
  b_ = 0.0
  x_ = 0.0
  call readMat(myFileName)
  P_ = mat
  !P_ = RESHAPE( SOURCE = (/ 1., 2., 3., 1., &
  !                         -3., 2., 5., 1., &
  !                         p_i, e_, -sqrt(2.0), 1., &
  !                         2., 4., 6., 2. /), SHAPE = (/ 4,4 /) )

  b_ = RESHAPE([0.,0.,0.,0.], SHAPE = [4,1])
  ! Gaussian elimination
  call print_mat(P_)
  call print_mat(b_)
  call partGaussElim(P_, b_, [4,4], [4,1], SINGLR)


End Program Driver_LinAl
