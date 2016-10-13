subroutine lglnodes(x,w,Pout,Pxout,n,np) 
  !
  ! 
  ! lglnodes.m
  !
  ! Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde 
  ! matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
  ! integration and spectral methods. 
  !
  ! Reference on LGL nodes and weights:  
  !   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
  !   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
  !
  ! Written by Greg von Winckel - 04/17/2004
  ! Contact: gregvw@chtm.unm.edu
  !
  !
  use type_defs
  implicit none
  integer :: n,n1,np
  real(kind = dp) :: w(0:n), x(0:n), xold(0:n)
  integer :: i,k
  real(kind = dp) :: P(1:n+1,1:n+1),eps
  real(kind = dp) :: Pout(0:n,0:np),Pxout(0:n,0:np)
  ! Truncation + 1
  N1 = N+1
  eps = 2.2204d-16
  
  ! Use the Chebyshev-Gauss-Lobatto nodes as the first guess
  do i = 0,n
     x(i) = -cos(pi*real(i,dp)/real(N,dp))
  end do
  
  ! The Legendre Vandermonde Matrix
  !  P=zeros(N1,N1);
  
  ! Compute P_(N) using the recursion relation
  ! Compute its first and second derivatives and 
  ! update x using the Newton-Raphson method.
  
  xold = 2.0_dp
  
  ! do while (maxval(abs(x-xold)) .gt. eps) 
  do i = 1,100   
     xold = x
     P(:,1) = 1.0_dp
     P(:,2) = x
     do  k=2,n
        P(:,k+1)=( real(2*k-1,dp)*x*P(:,k)-real(k-1,dp)*P(:,k-1) )/real(k,dp)
     end do
     x = xold-( x*P(:,N1)-P(:,N))/(real(N1,dp)*P(:,N1))
     if (maxval(abs(x-xold)) .lt. eps ) exit
  end do
  w=2.0_dp/(real(N*N1,dp)*P(:,N1)**2)
  
  Pout(:,0) = 1.0_dp
  Pout(:,1) = x
  ! k*P_k(x)=(2*k-1)*x*P_{k-1}(x)-(k-1)*P_{k-2}(x)
  do k=2,np
     Pout(:,k)=(real(2*k-1,dp)*x*Pout(:,k-1)-real(k-1,dp)*Pout(:,k-2) )/real(k,dp)
  end do

  Pxout(:,0) = 0.0_dp
  Pxout(:,1) = 1.0_dp
  ! P_k'(x)=P_{k-2}'(x)+(2k-1)P_{k-1}(x)
  do k=2,np
     Pxout(:,k) = Pxout(:,k-2)+real(2*k-1,dp)*Pout(:,k-1)
  end do

end subroutine lglnodes
