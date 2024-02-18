module LinAl
  implicit none
  
  integer, save :: msize, nsize
  real, dimension(:,:), allocatable, save :: mat
  

contains

  !********************************************************
  subroutine read_mat(A, myFileName)

    real, dimension(:,:), allocatable :: A
    character(len=100) :: myFileName
    ! Read matrix A
    ! myFileName = 'Amat.dat'

    open(10,file=myFileName)
    read(10,*) msize,nsize
    close(10)


    allocate(mat(msize,nsize))
    allocate(A(msize,nsize))
    ! Always initialize with zeros
    mat = 0.0
    
    call readMat(myFileName)
    A = mat
    ! sizeA(1) = msize
    ! sizeA(2) = nsize
    deallocate(mat)


  end subroutine read_mat

  subroutine print_mat(A)
    implicit none
    real, dimension(:,:), intent(IN) :: A

    integer :: i, j, n, m

    m = size(A, dim=1)
    n = size(A, dim=2)

    print *
    print '(i5,a5,i5)', m, "x",n
    do i = 1, m
      do j = 1,n
        write(6,"(F11.4)", ADVANCE="NO") A(i,j)
      end do
      print *
    end do
    print *, ""

  end subroutine print_mat

  subroutine print_vec(v)
    ! for integers pass in vec + 0.0
    implicit none
    real, dimension(:), intent(IN) :: v

    integer :: i, j, n

    n = size(v)

    print *
    print '(i5)', n
    do j = 1,n
      write(6,"(F11.4)", ADVANCE="NO") v(j)
    end do
    print *
    print *, ""

  end subroutine print_vec


  !********************************************************

  subroutine readMat(filename)

    implicit none
    character(len=*) :: filename

    integer :: i,j

    ! Reads a file containing the matrix A 
    ! Sample file:
    !
    ! 4 4 
    ! 2.0 1.0 1.0 0.0
    ! 4.0 3.0 3.0 1.0
    ! 8.0 7.0 9.0 5.0
    ! 6.0 7.0 9.0 8.0
    !
    ! Note that the first 2 numbers in the first line are the matrix dimensions, i.e., 4x4,
    ! then the next msize lines are the matrix entries. This matrix is found in Eq. 2.18 of the lecture note.
    ! Note that entries must be separated by a tab.


    open(10,file=filename)

    ! Read the matrix dimensions
    read(10,*) i,j

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)
    
  end subroutine readMat

  ! read image with msize, and nsize already initialized
  subroutine readMat2(myFileName)

    ! real, dimension(:,:), allocatable :: A
    character(len=100) :: myFileName
    integer :: i, j


    allocate(mat(msize,nsize))
    ! allocate(A(nsize,msize))

    ! Always initialize with zeros
    mat = 0.0

    open(10,file=myFileName)

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)
    ! deallocate(mat)


  end subroutine

  subroutine printMat(A, dims)

    ! Write a function (and/or subroutine) that takes two input arguments of
    ! (i) an m × n matrix A, and (ii) both of its dimensions. This routine then
    ! prints the matrix and its dimensions (i.e., m and n) to the screen in a
    ! human readable form (i.e., screen output).

    implicit none
    integer, dimension(:), intent(IN) :: dims
    real, dimension(:,:), intent(IN) :: A
    integer :: i, j

    ! write(*,*) "Dimensions"
    write(*,*) dims(1), "x", dims(2)
    ! write(*,*) "Matrix"

    do i = 1, dims(1)
      write(*,*) (A(i,j) , j = 1, dims(2) )
    end do

  end subroutine printMat


  subroutine traceMat(A, m, trace)

    implicit none
    integer, intent(IN) :: m
    real, dimension(:,:), intent(IN) :: A

    real, intent(OUT) :: trace
    integer :: i

    ! A: m × m square matrix
    ! m: the first dimension of A
    ! trace: return value, the trace of A

    do i=1,m
       trace = trace + A(i,i)
    enddo
    
  end subroutine traceMat

  subroutine printColumnNorm(A)

    implicit none
    real, dimension(:,:), intent(IN) :: A

    real :: norm
    integer :: j, nsize, msize

    ! Calculate two norm of columns of A_(m x n)

    nsize = size(A, dim=2)
    msize = size(A, dim=1)
    norm = 0.0
    print *, ""
    do j = 1, nsize
      norm = 0.0
      call twoNorm(A(:,j), nsize, norm)
      write(*,*) "Norm column: ", j
      write(*,*) norm
    end do

  end subroutine printColumnNorm

  subroutine two_norm(vec, n, norm)

    implicit none
    real, dimension(:,:), intent(IN) :: vec
    integer, intent(IN) :: n

    real, intent(OUT) :: norm

    integer :: i

    norm = 0.0

    ! Write a function (and/or subroutine) that takes three arguments. Two
    ! input arguments are a vector and its dimension, and one output argument
    ! is its Euclidean norm (i.e., 2-norm).

    ! TODO: allow this to work for size vec(i) and vec(i,1)
    do i=1,n
       norm = norm + vec(i,1)**2
    enddo
    norm = sqrt(norm)

  end subroutine two_norm

  subroutine twoNorm(vec, n, norm)

    implicit none
    real, dimension(:), intent(IN) :: vec
    integer, intent(IN) :: n

    real, intent(OUT) :: norm

    integer :: i

    norm = 0.0

    ! Write a function (and/or subroutine) that takes three arguments. Two
    ! input arguments are a vector and its dimension, and one output argument
    ! is its Euclidean norm (i.e., 2-norm).

    ! TODO: allow this to work for size vec(i) and vec(i,1)
    do i=1,n
       norm = norm + vec(i)**2
    enddo
    norm = sqrt(norm)

  end subroutine twoNorm

  subroutine partGaussElim(A, B, dimsA, dimsB, SINGLR)

    implicit none
    integer, dimension(:) :: dimsA
    integer, dimension(:) :: dimsB
    real, dimension(:,:) :: A
    real, dimension(:,:) :: B
    real, dimension(:,:), allocatable :: As ! copy of A
    real, dimension(:,:), allocatable:: Bs ! copy of B
    real, dimension(:,:), allocatable:: X ! solution 
    real, dimension(:,:), allocatable :: D ! augemented matrix
    real, dimension(:), allocatable :: temp ! for swapping rows
    integer :: i, j, m, kmax 
    integer, dimension(2) :: dimsD
    integer, dimension(2) :: dimsX
    real :: maximum
    logical :: SINGLR ! flag for whether matrix A is invertible 
  

    ! Should probably assert dimsA(1) = dimsA(2) = dimsB(1)
    m = dimsA(1)
    ! X = 0.0 ! set = to 0?

    ! Concatanate matrices (augmented matrix)
    dimsD(1) = dimsA(1)
    dimsD(2) = dimsA(2) + dimsB(2)

    allocate(As(dimsA(1),dimsA(2)))
    allocate(Bs(dimsB(1),dimsB(2)))
    allocate(X(dimsA(2), dimsB(2)))
    X = 0.0
    dimsX(1) = dimsA(2)
    dimsX(2) = dimsB(2)
    As = A ! copy A
    Bs = B ! copy B
    allocate(D(dimsD(1),dimsD(2)))
    allocate(temp(dimsD(2)))
    do i=1,dimsD(2)
      if (i <= dimsA(2)) then
        D(:,i)=A(:,i)
      else
        D(:,i)=B(:,i-dimsA(2))
      end if
    enddo

    ! Gaussian Elimination
    do i = 1,m
      ! need to add +(i-1) to account for the array "shrinking"
      kmax = maxloc(abs(D(i:m,i)), dim=1) + (i-1)
      maximum = D(kmax, i)
    !   ! if (maximum < 1e-14*norma) then !singular

      ! Swap rows so pivot is maximum
      if (i /= kmax) then 
        temp = D(kmax,:)
        D(kmax,:) = D(i,:)
        D(i,:) = temp
      end if

      ! Matrix is singular (stop)
      if (D(i, i) == 0) then
        SINGLR = .TRUE.
        stop
      end if

      ! Eliminiation step
      do j = i+1,m
        D(j,:) = D(j,:) - D(j,i)*(D(i,:)/D(i,i))
      enddo

    enddo

    A = D(:, 1:m)
    B = D(:, m+1:dimsB(2))
    deallocate(D)


    ! I still need to teach myself how these print statements work..
    ! print A
    print *, ""
    print *, ""
    print *, "matrix A after gauss"
    ! call printMat(A, dimsA)
    call print_mat(A)
    ! print B
    print *, ""
    print *, ""
    print * , "matrix B after gauss"
    !call printMat(B, dimsB)
    call print_mat(B)
    print *, SINGLR
    call backSubstitutionU(A, B, X)
    print *, ""
    print *, ""
    print *, "solution X to AX = B"
    !call printMat(X, dimsX)
    call print_mat(X)


    print *, ""
    print *, ""
    print *, "Error Matrix"
    ! call printMat2(matmul(As,X)-Bs)
    call print_mat(matmul(As,X)-Bs)

    print *, ""
    print *, ""
    print *, "Two norm of columns of error matrix"
    call printColumnNorm(matmul(As,X)-Bs)

  end subroutine partGaussElim

  subroutine LUdecomp(A, L, U)
    real, dimension(:,:) :: A
    real, dimension(:,:) :: L
    real, dimension(:,:) :: U
    real :: sum_
    integer :: i, j, k, n
    L = 0.0
    U = 0.0
    n = size(A, dim=1)

    do i = 1, n

      ! upper triangular
      do k = i, n
        ! summation of L(i,j) * U(j,k)
        sum_ = 0.0
        do j = 1, i
          sum_ = sum_ + L(i,j) * U(j,k)
        enddo

        ! evaluating U(i,k)
        U(i,k) = A(i,k) - sum_
      enddo

      ! lower triangular
      do k = i, n
        if (i == k) then
          L(i,i) = 1
        else
          ! summation of L(k,j) * U(j,i)
          sum_ = 0.0
          do j = 1, i
            sum_ = sum_ + L(k,j) * U(j,i)
          enddo
          L(k,i) = (A(k,i) - sum_)/U(i,i)
        end if
      enddo

    enddo

  end subroutine LUdecomp

  subroutine LU_partial(A, m, s, SINGLR)
    real, dimension(:,:) :: A
    integer, intent(IN) :: m
    integer, dimension(:) :: s
    logical :: SINGLR

    real, dimension(m) :: temp ! for swapping rows
    integer :: i, j, k, kmax, temp_s
    real :: maximum, l_
    real, dimension(m,m) :: L, U


    L = 0.0
    ! U = A
    
    do i = 1,m
      s(i) = i
    enddo
    ! Gaussian Elimination
    do i = 1,m
      ! need to add +(i-1) to account for the array "shrinking"
      kmax = maxloc(abs(A(i:m,i)), dim=1) + (i-1)
      maximum = A(kmax, i)
    !   ! if (maximum < 1e-14*norma) then !singular

      ! Swap rows so pivot is maximum
      if (i /= kmax) then 
        temp = A(kmax,:)
        A(kmax,:) = A(i,:)
        A(i,:) = temp

        temp_s = s(i)
        s(i) = s(kmax)
        s(kmax) = temp_s

        ! swap L
        if (i > 1) then
          temp = L(kmax,:)
          L(kmax,:) = L(i,:)
          L(i,:) = temp
        end if
      end if


      ! Matrix is singular (stop)
      if (A(i, i) == 0) then
        SINGLR = .TRUE.
        stop
      end if

      ! Eliminiation step
      L(i,i) = 1.0
      do j = i+1,m
        L(j,i) = A(j,i)/A(i,i)
        A(j,:) = A(j,:) - A(j,i)*(A(i,:)/A(i,i))
      enddo
    enddo

    ! Add L to A
    do i = 1,m
      do j = i+1,m
        A(j,i) = L(j,i)
      enddo
    enddo

  end subroutine LU_partial

  subroutine backSubstitutionLU(A, m, B, s, X)
    real, dimension(:, :) :: A, B, X
    real, dimension(:, :), allocatable :: Y, Bs
    integer, dimension(:) :: s
    real, dimension(m, m) :: L, U
    integer :: m, i, j

    allocate(Y(m, size(B, dim=2)))
    allocate(Bs(m, size(B, dim=2)))
    Y = 0.0
    Bs = 0.0

    B = B(s,:)
    Bs = B ! copy of B

    do i = 1,m
      do j = 1,i
        U(j,i) = A(j,i)
      enddo
      L(i,i) = 1.0
      do j = i+1,m
        L(j,i) = A(j,i)
      enddo
    enddo

    call backSubstitutionL(L,B,Y)
    call backSubstitutionU(U,Y,X)

    ! print *, ""
    ! print *, ""
    ! print *, "Solution AX = B using LU decomposition: "
    ! call print_mat(X)

    ! print *, ""
    ! print *, ""
    ! print *, "Two norm of error matrix using LU decomposition"
    ! ! call printMat2(matmul(As,X)-Bs)

    ! call printColumnNorm(matmul(matmul(L,U),X)-Bs)

  end subroutine backSubstitutionLU


  ! Back substitution when L is upper triangular
  subroutine backSubstitutionL(L,B,X)
    implicit none
    real, dimension(:, :) :: L
    real, dimension(:, :) :: B
    real, dimension(:, :) :: X
    integer :: n, i, j

    n = size(B,dim=1)

    do i=1,n ! start, stop [,step]

      X(i,:) = B(i,:)/L(i,i) 

      do j = i+1,n
        B(j,:) = B(j,:) - X(i,:)*L(j,i)
      enddo
      ! B(1:i-1,:) = B(1:i-1,:) - X(i,:)*U(1:i-1,i)
    enddo


  end subroutine backSubstitutionL


  ! Back substitution when U is upper triangular
  subroutine backSubstitutionU(U,B,X)
    ! Solve UX = B for X
    implicit none
    real, dimension(:, :) :: U
    real, dimension(:, :) :: B
    real, dimension(:, :) :: X
    integer :: n, i, j

    n = size(B,dim=1)

    do i=n,1,-1 ! start, stop [,step]

      X(i,:) = B(i,:)/U(i,i) 

      do j = 1,i-1
        B(j,:) = B(j,:) - X(i,:)*U(j,i)
      enddo
      ! B(1:i-1,:) = B(1:i-1,:) - X(i,:)*U(1:i-1,i)
    enddo


  end subroutine backSubstitutionU

  subroutine choleskyDecom(A, flag)
    ! Replace A with an upper triangular matrix L^T
    implicit none
    real, dimension(:,:) :: A
    real, dimension(:,:), allocatable :: R
    logical :: flag
    real, dimension(:,:), allocatable :: v
    real :: h
    integer :: n,j

    ! should probably assert n = m (square)
    n = size(A, dim=1)
    allocate(R(n,n))
    R = 0.0 ! always initialize to 0!
    allocate(v(1,n))

    do j =1,n
      v = A(j:j,j:n)
      if (j>1) then
        ! v = A(j,j:n) - matmul(R(j-1:1:-1,j),R(1:j-1,j:n)) ! Why can't we just use transpose()?
        v = A(j:j,j:n) - matmul(transpose(R(1:j-1,j:j)), R(1:j-1,j:n))
      end if
      if (v(1,1)<=0) then
        write(*,*) "matrix is not positive definite"
        flag = .TRUE.
        stop
      else
        h = 1/sqrt(v(1,1))
      end if
      R(j:j,j:n) = v*h
    enddo

    A = R
    deallocate(R)
    deallocate(v)


  end subroutine choleskyDecom

  subroutine backSubChol(A, B, X, dim)
    ! backsubstitution routine that solves LL T x = b. The routine
    ! takes as the following information as arguments:
    ! • the matrix A already Cholesky-decomposed (input)
    ! • its first dimension (input)
    ! • a vector b (Note : you can generalize this and input a matrix B containing
    ! n rhs vectors if you prefer, i.e., B = [b 1 t b 2 t · · · t b n ], b i ∈ R m ) (input)

    real, dimension(:,:) :: A ! L^t => upper triangular
    real, dimension(:,:) :: B ! AX=B
    real, dimension(:,:) :: X ! AX=B
    integer :: dim
    ! real, dimension(dim,size(B,dim=2)) :: X ! Solution
    real, dimension(dim,size(B,dim=2)) :: Y ! Solution

    
    ! First solve L(Y) = B
    call backSubstitutionL(transpose(A), B, Y)
    ! Then solve L^t(X) = Y
    call backSubstitutionU(A, Y, X)


    ! Solution
    ! print *, "X"
    ! call printMat2(X)

  end subroutine backSubChol


  subroutine householderQR(A, d, dimsA)
    ! calculates full Q_mxm, R_mxn

    ! Input A, d, dimensions of A
    ! Computes the implicit QR. Matrix A is replaced by the upper triangular matrix R and the normalized 
    ! householder vectors are
    ! stored in the lower half of A. The diagonal R is stored in d (column vector).

    implicit none
    real, dimension(:,:) :: A
    real, dimension(:,:) :: d ! size n?
    integer, dimension(:) :: dimsA
    integer :: m, n, j
    real :: alpha, nrm
    m = dimsA(1) 
    n = dimsA(2)

    ! Iterate over columns of A
    do j = 1,n
      call twoNorm(A(j:m,j), m-j+1, alpha) ! stores the norm in alpha

      if (A(j,j) >= 0) then
        d(j,1) = -alpha
      else
        d(j,1) = alpha
      end if

      nrm = sqrt(alpha*(alpha+abs(A(j,j)))) ! nrm = ||v||
      A(j,j) = A(j,j) - d(j,1);
      A(j:m,j) = A(j:m,j)/nrm ! store u = v/||v||
      ! A(j:m,j) = A(j:m,j) ! store v

      ! transform the rest of the matrix A := A-u*(u'*A)
      if (j < n) then
        A(j:m,j+1:n) = A(j:m, j+1:n) - matmul(A(j:m,j:j), matmul(transpose(A(j:m,j:j)), A(j:m, j+1:n)))
      end if 
    enddo
    
  end subroutine householderQR

  subroutine HouseholderQTy(A_qr, y)
    ! replaces y <- Q^Ty
    implicit none
    real, dimension(:,:), intent(IN) :: A_qr ! size mxn
    real, dimension(:,:) :: y

    ! local variables
    integer :: m, n, j
    m = size(A_qr, dim=1); n = size(A_qr, dim=2)

    do j = 1,n
      ! z(j:m)=z(j:m)-A(j:m,j)*(A(j:m,j)’*z(j:m));
      y(j:m,1) = y(j:m,1) - matmul(matmul(A_qr(j:m,j:j),transpose(A_qr(j:m,j:j))), y(j:m,1))
    enddo
  end subroutine HouseholderQTy

  subroutine householderQy(A, y, dimsA)
    ! replaces y <- Qy
    ! computes Qy using the Householder
    ! reflections Q stored as vectors in the matrix A by
    ! housholderQR(A). Replace input y with Qy
    implicit none
    real, dimension(:,:) :: A ! size mxn
    real, dimension(:,:) :: y ! size nx1
    integer, dimension(:) :: dimsA
    integer :: m, n, j
    m = dimsA(1) 
    n = dimsA(2)

    do j = n,1,-1 ! start, stop [,step]
      y(j:m,1) = y(j:m,1) - matmul(matmul(A(j:m,j:j),transpose(A(j:m,j:j))),y(j:m,1))
    enddo

  end subroutine householderQy

  subroutine eQR(A, Q, R)
    ! explicitely calculate QR decomposition. Takes input A mxn, Q mxm, R mxn and outputs QR factorization
    ! store in Q and R. Leaves A unchanged.

    real, dimension(:,:) :: A ! mXn
    real, dimension(:,:) :: Q ! mXm
    real, dimension(:,:) :: R ! mXn
    real, dimension(:,:), allocatable :: dR, Id
    integer, dimension(2) :: dimsA
    integer :: msize, nsize, i

    dimsA(1) = size(A, dim=1)
    dimsA(2) = size(A, dim=2)
    msize = dimsA(1)
    nsize = dimsA(2)

    allocate(dR(nsize, 1))
    allocate(Id(msize, msize)) ! identity matrix
    Id(1:msize,1:msize) = 0.0
    forall (i = 1:msize) Id(i,i) = 1.0

    call householderQR(A, dR, dimsA)

    ! create Q and R explicitely (Make this a subroutine?)
    do i = 1,msize
      ! Apply Q onto the columns of the Identity matrix, using the orthonormal vectors, u,
      ! from the Householder reflector, H = I-uu^T (which are stored in the lower part of E)
      call householderQy(A, Id(:,i:i), dimsA)
      Q(:,i:i) = Id(:,i:i)
      if (i <= nsize) then
        R(:, i) = A(1:i,i) ! upper triangular part of E stores R (excluding the diagonal part)
        R(i,i) = dR(i,1) ! diagonal part of R
      end if
    enddo
    deallocate(Id)
    deallocate(dR)

  end subroutine

  subroutine expliciteQR(A_qr, dR, Q, R)
    ! create Q and R explicitely
    ! Apply Q onto the columns of the Identity matrix, using the orthonormal vectors, u,
    ! from the Householder reflector, H = I-uu^T (which are stored in the lower part of A_qr)
    implicit none
    real, dimension(:,:), intent(IN) :: A_qr, dR ! size mxn
    real, dimension(:,:), intent(OUT) :: Q, R ! Q_mxmH

    ! local variables
    real, dimension(:,:), allocatable :: id
    integer :: m, n, i
    m = size(A_qr, dim=1); n = size(A_qr, dim=2)
    
    allocate(Id(m,m)) ! identity matrix
    Id(1:m,1:m) = 0.0
    forall (i = 1:m) Id(i,i) = 1.0

    ! create Q and R explicitely
    do i = 1,m
      call householderQy(A_qr, Id(:,i:i), [m,n])
      Q(:,i:i) = Id(:,i:i)
      if (i <= n) then
        R(:, i) = A_qr(1:i,i) ! upper triangular part of E stores R (excluding the diagonal part)
        R(i,i) = dR(i,1) ! diagonal part of R
      end if
    enddo
  end subroutine expliciteQR

  subroutine hessenberg(A)

    ! uses householder transformations to compute hessenberg decomposition of A
    ! factors a symmetric matrix A into a tridiagonal matrix

    real, dimension(:,:) :: A ! size mxn
    real, dimension(:,:), allocatable :: v ! size mxn
    integer :: m, n, j
    real :: alpha, nrm
    m = size(A, dim=1) 
    n = size(A, dim=2) 

    allocate(v(m,1))

    ! Iterate over columns of A
    do j = 1,n-1
      v = 0.0

      call twoNorm(A(j+1:m,j), m-j, alpha) ! stores the norm in alpha
      if (A(j+1,j) >= 0) then
        alpha = -alpha
      end if

      ! compute householder vector
      v(j+1:m,1) = A(j+1:m,j)
      v(j+1,1) = A(j+1,j) - alpha
      ! normalize vector
      nrm = 0.0
      call twoNorm(v(1:m,1), m, nrm) ! nrm = ||v||
      v(1:m,1) = v(1:m,1)/nrm

      ! update A from the left
      A = A - 2 * matmul(matmul(v, transpose(v)),A)
      ! update A from the right
      A = A - 2 * matmul(matmul(A, v),transpose(v))
    enddo
    
    deallocate(v)
  end subroutine

  subroutine eigQR(A)

    ! Without shifts computes the eigenvalues of the (tridiagonal) matrix A

    real, dimension(:,:) :: A ! size nxn
    real, dimension(:,:), allocatable :: Q, R, diag
    ! ! vector to store eigenvalues
    ! real, dimension(:,:), allocatable :: v ! size nx1
    integer :: n, j, i, num_iter
    real :: norm, tol
    n = size(A, dim=1) 

!   allocate(Q(msize, msize)) ! Q matrix
    allocate(Q(n, n)) ! Q matrix
    Q = 0.0
  ! allocate(R(msize, nsize)) ! R matrix
    allocate(R(n, n)) ! R matrix
    R = 0.0
    allocate(diag(n-1, 1)) ! identity matrix
    diag = 0.0
    norm  = 1.0
    tol = 10.0**(-12.0)

    num_iter = 0
    ! do i = 1,100 ! loop until error is small
    do while (norm >= tol) ! norm of subdiagonal
    num_iter = num_iter + 1
      ! calculate A = QR
      call eQR(A, Q, R)
      ! update A = RQ
      A = matmul(R, Q)

      diag = diag(1:n-1,1:1)
      forall (i = 1:n-1) diag(i,1) = A(i+1,i)
      call two_norm(diag, n-1, norm)
    enddo
    print *, "num iterations"
    print *, num_iter
    ! enddo

    ! call printMat2(A)

  end subroutine


  subroutine eigQRshift(A)

    ! with shift computes the eigenvalues of A

    real, dimension(:,:) :: A ! size nxn
    real, dimension(:,:), allocatable :: Q, R, Id, B
    ! ! vector to store eigenvalues
    ! real, dimension(:,:), allocatable :: v ! size nx1
    integer :: n, j, i, m, num_iter
    real :: mu, tol, last_subdiag
    m = size(A, dim=1) 
    tol = 10.0**(-12.0)

!   allocate(Q(msize, msize)) ! Q matrix
    allocate(Q(m, m)) ! Q matrix
    Q = 0.0
  ! allocate(R(msize, nsize)) ! R matrix
    allocate(R(m, m)) ! R matrix
    R = 0.0
    allocate(B(m, m)) ! identity matrix
    B = 0.0
    allocate(Id(m, m)) ! identity matrix
    Id(1:m,1:m) = 0.0
    forall (i = 1:m) Id(i,i) = 1.0


    num_iter = 0
    do j = m,2,-1 ! loop until error is small ( norm of super(sub) diagonal?)
    ! Loop while error is small
      last_subdiag = 1.0
      do while (last_subdiag >= tol)

        B = B(1:j,1:j)
        R = R(1:j,1:j)
        Q = Q(1:j,1:j)
        Id = Id(1:j,1:j)
        num_iter = num_iter + 1

        ! calculate A = QR
        B = A(1:j,1:j)
        mu = B(j, j)
        call eQR(B - mu*Id, Q, R)
        ! update A = RQ
        B = matmul(R, Q) + mu*Id
        A(1:j,1:j) = B
        last_subdiag = abs(B(j, j-1))
      enddo
    enddo
    print *, "num iterations"
    print *, num_iter

end subroutine

subroutine inverseIter(A, mu, x)

  ! Calculates eigenvectors of a matrix A given a "guess" mu approximately equal to an eigenvalue
  ! replaces B = A-mu*Id with LU decomposition and stores pivots in s
  ! solves B(xn) = x for xn and calculates two norm error x-xn for the stopping criteria. 

    ! input matrix
    real, dimension(:,:) :: A ! size nxn
    ! x: return value (eigenvector of A)
    real, dimension(:,:) :: x ! initialize default [1 1 ... 1]/sqrt(n)

    ! B = A-mu*Id => same eigenvectors as A
    real, dimension(:,:), allocatable :: Id, B
    real, dimension(:,:), allocatable :: xn, x_copy ! eigenvector
    integer, dimension(:), allocatable :: s ! pivots in LU decomposition
    real, dimension(:,:), allocatable :: err 
    real :: nrm, nrm_err
    logical :: SINGLR


    integer :: n, j, i, num_iter
    real :: mu, tol
    n = size(A, dim=1) 

  ! iterates to eigenvector
    ! allocate(x(n,1))
    allocate(xn(n,1))
    ! allocate(y(n,1))
    allocate(s(n))
    allocate(err(n,1))

    allocate(x_copy(n,1))
    x_copy = 0.0

    ! x = 1.0/sqrt(n + 0.0)
    ! initialize error above tolerance
    nrm_err = 1.0
    tol = 10.0**(-12.0)

    allocate(Id(n, n)) ! identity matrix
    Id(1:n,1:n) = 0.0
    forall (i = 1:n) Id(i,i) = 1.0

    allocate(B(n,n))
    B = A - mu*Id

    SINGLR = .false.
    call LU_partial(B, n, s, SINGLR)
    num_iter = 0
    do while (nrm_err >= tol)
      num_iter = num_iter + 1
    ! solve Bxn = LL^t xn = x
      ! First solve L(y) = x then solve L^t(xn) = y
      x_copy = x !~ routine changes x so we need a copy
      call backSubstitutionLU(B, n, x, s, xn)
      
      call two_norm(xn, n, nrm)
      xn = xn/nrm

      err = abs(xn) - abs(x_copy) ! -x is an eigenvector too
      call two_norm(err, n, nrm_err)
      x = xn
    enddo

    print *, "Number of iterations: "
    print *, num_iter

    deallocate(B)
    deallocate(s)
    deallocate(xn)
    deallocate(Id)
    deallocate(x_copy)
    deallocate(err)


  end subroutine

  subroutine printMatToFile(A, filename)

    implicit none
    real, dimension(:,:) :: A
    character(len=*) :: filename
    integer :: i, j, n, m

    open(10,file=filename)

    m = size(A, dim=1)
    n = size(A, dim=2)

    ! Print in matrix form
    ! write(10,*) ""
    ! ! print '(i5,a5,i5)', m, "x",n
    ! do i = 1, m
    !   write(10,*) (A(i,j) , j = 1, n )
    ! enddo


    ! Print dimensions 
    ! print '(i5,i5)', m, n
    ! write(10,*) m
    ! write(10,*) n

    ! Print as column and then reshape later
    ! write(10,*) (A(j,i), i=1,m, j=1,n) ! switching i and j works better with reshape
    do i = 1, m
        ! write(10,*) (A(j,i), j=1,n) ! switching i and j works better with reshape
      do j = 1,n
        write(10,*) A(i,j)  ! switching i and j works better with reshape
        ! write(10, "(ES28.16)", ADVANCE="NO") A(j,i)
      end do
      ! write(10,*)
      ! write(10,*) ""
    end do
    ! write(10,*) ""


    close(10)


  end subroutine 
  
  subroutine Jacobi(A, b, x_0, x, acc) 

    ! Uses Jacobi iterative algorithm to solve the system Ax = b
    !
    ! Inputs:
    !   A (dimension(m,m)): To be be decomposed as A = D + L + U
    !   b (dimension(m,1)): 
    !   x_0 (dimension(m,1)): Inital guess
    !   acc (double precision): Desired accuracy of algorithm
    ! Input/Output:
    !   x (dimension(m,1)): Solution to system

    double precision, dimension(:,:) :: A, b, x_0, x
    double precision, dimension(:,:), allocatable :: temp
    ! double precision, dimension(:,:), allocatable :: x ! guess
    double precision :: acc, r, s
    integer i, j, m, iter, d_
    character(len=100) :: filename

    d_ = A(1,1) ! used for filenaming

    m = size(A, dim = 1)

    ! allocate(x(m,1)) ! 
    allocate(temp(m,1)) ! size of x (and b)
    temp = 0.0
     
    x = x_0

    call two_norm( matmul(A,x) - b, m,  r)

    ! write (filename, "(A11,I2,A4)") "jac_error_d", d_, ".dat"
    open(10, file ='jac_error.dat') 


    iter = 0
    ! do while ( iter < 1000 ) 
    do while (r > acc) 

      ! write error to file as long as its not huge (can exit long before this)
      if (r < huge(r)) then
        write(10, *) r ! write error to file
      else
        print*, "Error is too large. This algorithm does not converge"
        exit ! no point in continuing 
      end if

      iter = iter + 1
      if (iter > 1000) then
        print *, "This is taking too long. This does not converge"
        exit
      end if

      ! do one iteration over rows
      do i = 1,m
        s = 0.0
        ! Lower part
        do j = 1, i-1
           s = s + A(i, j)*x(j, 1)
        enddo
        ! Upper part
        do j = i+1, m
           s = s + A(i, j)*x(j, 1) 
        enddo
        temp(i,1) = (b(i,1) - s)/(A(i,i))
      enddo
      x = temp

      call two_norm( matmul(A,x) - b, m,  r)
    end do

    close(10) 

    print *, "number of iterations: "
    print *, iter

  end subroutine Jacobi

  subroutine GaussSeidel(A, b, x_0, x, acc) 
    ! Uses Gauss-Seidel iterative algorithm to solve the system Ax = b. Uses
    ! updated values as soon as they become available.
    !
    ! Inputs:
    !   A (dimension(m,m)): To be be decomposed as A = D + L + U
    !   b (dimension(m,1)): Vector b in the system Ax = b
    !   x_0 (dimension(m,1)): Inital guess
    !   acc (double precision): Desired accuracy of algorithm
    ! Input/Output:
    !   x (dimension(m,1)): Solution to system

    double precision, dimension(:,:) :: A, b, x_0, x
    double precision :: acc
    
    ! Local variables 
    !   r: normed residual (exit when this is within desired accuracy)
    !   s: used to calculate products
    double precision :: r, s
    integer i, j, m, iter, maxIter


    m = size(A, dim = 1)

    x = x_0 ! set our solution our initial guess


    call two_norm( matmul(A,x) - b, m,  r)

    open(10, file='gseid_error.dat')
    iter = 0
    ! do while ( iter < 1000 ) 
    do while (r > acc)

      ! Exit criteria
      ! write error to file as long as its not huge (can exit long before this)
      if (r < huge(r)) then
        write(10, *) r ! write error to file
      else
        print *, "Error is too large. This does not converge"
        exit ! no point in continuing 
      end if

      ! stop after 1000 iterations
      iter = iter + 1
      if (iter > 1000) then
        print *, "This is taking too long. Algorithm does not converge"
        exit
      end if

      ! do one iteration over rows
      do i = 1,m
        s = 0.0
        ! Lower part
        do j = 1, i-1
           s = s + A(i, j)*x(j, 1)
        enddo
        ! Upper part
        do j = i+1, m
           s = s + A(i, j)*x(j, 1) 
        enddo
        x(i,1) = (b(i,1) - s)/(A(i,i))
      enddo

      call two_norm( matmul(A,x) - b, m,  r)
    end do

    close(10)

    print *, "number of iterations: "
    print *, iter

    ! call printMat2(x)


  end subroutine GaussSeidel

  subroutine ConjugateGradient(A, b, x_0, x, acc)
    double precision, dimension(:,:) :: A, b, x_0, x
    double precision :: acc
    ! Local variables
    double precision, dimension(:,:), allocatable :: A_copy
    double precision, dimension(:,:), allocatable :: p ! direction
    double precision, dimension(:,:), allocatable :: r ! residual
    double precision, dimension(:,:), allocatable :: y ! temporary vector
    double precision, dimension(1,1) :: pty
    ! Local variables
    double precision :: E, E_new, alpha, beta
    integer i, j, m, iter

    


    m = size(A, dim=1)

    ! Check that A is symmetric positive definite
    allocate(A_copy(m,m))
    A_copy = A

    ! Create a subroutine?
    ! Check if symmetric and
    ! Check positive definite (calculate eigenvalues)
    ! First put A_copy in hessenberg form and run the QR algorithm with shift
    ! to reveal eigenvalues
    call hessenberg(A_copy)
    call eigQRshift(A_copy) ! I need better stopping criteria for this algorithm
    do i = 1,m
      do j = 1,m
        if (.not.(A(i,j) == A(j,i))) then
          print *, "matrix is not symetric"
          exit ! set flag?
        endif
        if (A(j,j) <= 0.0) then
          print *, "matrix is not positive definite"
        end if
      enddo
    enddo

    deallocate(A_copy)

    open(10, file='conjgr_error.dat')
    ! Begin algorithm
    allocate(p(m,1))
    ! p = b
    allocate(r(m,1))
    x = x_0
    r = b - matmul(A,x)

    p = r
    call two_norm(r, m, E)
    E = E**(2)
    iter = 0


    do while (sqrt(E) > acc)
      iter = iter + 1

      ! Exit criteria
      ! write error to file as long as its not huge (can exit long before this)
      if (sqrt(E) < huge(sqrt(E))) then
        write(10, *) sqrt(E) ! write error to file
      else
        print *, "Error is too large. This does not converge"
        exit ! no point in continuing 
      end if

      if (iter > 1000) then
        print *, "This does not converge"
        ! stop
        exit
      end if

      ! "Smart"
      y = matmul(A, p)
      pty = matmul(transpose(p), y)
      alpha = E/pty(1,1)
      x = x + alpha* p
      r = r - alpha*y
      call two_norm(r, m, E_new)
      E_new = E_new**(2)
      beta = E_new/E
      p = r + beta*p
      E = E_new
    end do

    close(10)

    print *, "number of iterations"
    print *, iter

  end subroutine ConjugateGradient

end module LinAl
