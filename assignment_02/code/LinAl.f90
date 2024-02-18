module LinAl
  implicit none
  
  integer, save :: msize, nsize
  real, dimension(:,:), allocatable, save :: mat
  

contains

  !********************************************************

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
    implicit none
    integer, dimension(:), intent(IN) :: v

    integer :: i, j, n

    n = size(v)

    print *
    print '(i5)', n
    do j = 1,n
      write(6,"(i4)", ADVANCE="NO") v(j)
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

    print *, ""
    print *, ""
    print *, "Solution AX = B using LU decomposition: "
    call print_mat(X)

    print *, ""
    print *, ""
    print *, "Two norm of error matrix using LU decomposition"
    ! call printMat2(matmul(As,X)-Bs)

    call printColumnNorm(matmul(matmul(L,U),X)-Bs)

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

end module LinAl
