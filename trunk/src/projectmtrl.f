      subroutine PROJECTMTRL( ib,  ie,  jb,  je,  kb,  ke, 
     *                        iba, iea, jba, jea, kba, kea, 
     *                        rho, mu, lambda, dt,
     *                        h, cfl, vsmin, rhoscale, muscale,
     *                        lascale, info )
***********************************************************************
***
*** Projection of the material arrays to satisfy the following constraints
***
***    1. CFL condition max(vp) < h*cflmax/dt
***    2. Smallest resolvable speed, min(vs) > vsmin
***    3. Positive density rho > eps
***    4. Positive mu,  mu > eps
***    5. Non-negative lambda, lambda >= 0
***      
*** Input: ib, ie, jb, je, kb, ke - Declared dimensions of arrays
***        iba, iea, jba, jea, kba, kea - Apply projection to this subarray
***        rho - Density
***        mu  - Lam\'e parameter
***        lambda - Lam\'e parameter
***        h - Grid spacing
***        cfl - Maximum stable cfl number
***        vsmin - Smallest resolvable wave speed, condition is
***                not enforced if vsmin < 0
***        rhoscale, muscale, lascale - Scaling parameters for rho, mu, lambda
***
*** Output: rho, mu, lambda - Material is projected to satisfy the constraints.
***         info - 0   --> No correction was needed
***               -1   --> Projection was unsuccessful
***                1   --> Negative density or mu was corrected.
***                10  --> CFL limit was enforced
***                100 --> smallest resolvable wave speed was enforced.
***       ( can be a sum of 1, 10, 100 if more than one condition was enforced )
***
***********************************************************************

      implicit none
      integer ib, ie, jb, je, kb, ke, iba, iea, jba, jea, kba, kea
      integer i, j, k, info, zcorr, cflcorr, lowvscorr
      real*8 rho(ib:ie,jb:je,kb:ke), mu(ib:ie,jb:je,kb:ke)
      real*8 lambda(ib:ie,jb:je,kb:ke), acof, bcof, c1, c2
      real*8 dt, h, cfl, rhoscale, muscale, lascale
      real*8 vsmin, k1, k2, idet, a11, a12, a22

      acof = (cfl*h/dt)**2
      if( vsmin .gt. 0 )then
         bcof = vsmin**2
      else
         bcof = -1
      endif
      zcorr     = 0
      cflcorr   = 0
      lowvscorr = 0
      do k=kba,kea
         do j=jba,jea
            do i=iba,iea
               if( mu(i,j,k).le.0 )then
                  mu(i,j,k) = 1e-6*muscale
                  zcorr = 1
               endif
               if( rho(i,j,k).le.0 )then
                  rho(i,j,k) = 1e-6*rhoscale
                  zcorr = 1
               endif
               if( lambda(i,j,k).lt.0 )then
                  lambda(i,j,k) = 0
                  zcorr = 1
               endif
               c1 = acof*rho(i,j,k)-4*mu(i,j,k)-lambda(i,j,k)
               c2 = mu(i,j,k)-bcof*rho(i,j,k)
               if( c1.lt.0 .and. c2.ge.0 )then
*** correct first condition
                  k1 = c1/( (acof*rhoscale)**2+16*muscale**2+lascale**2)
                  mu(i,j,k) = mu(i,j,k) +muscale**2*4*k1
                  lambda(i,j,k) = lambda(i,j,k)+lascale**2*k1
                  rho(i,j,k) = rho(i,j,k) -rhoscale**2*acof*k1
                  c1 = acof*rho(i,j,k)-4*mu(i,j,k)-lambda(i,j,k)
                  c2 = mu(i,j,k)-bcof*rho(i,j,k)
                  cflcorr = 1
               endif

               if( c1.ge.0 .and. c2.lt.0 )then
*** correct second condition
                  k2 = c2/(muscale**2+(bcof*rhoscale)**2)
                  mu(i,j,k)  = mu(i,j,k)  - muscale**2*k2
                  rho(i,j,k) = rho(i,j,k) + rhoscale**2*bcof*k2
                  c1 = acof*rho(i,j,k)-4*mu(i,j,k)-lambda(i,j,k)
                  c2 = mu(i,j,k)-bcof*rho(i,j,k)
                  lowvscorr = 1
               endif

               if( c1.lt.0 .or. c2.lt.0 )then
*** correct both conditions
                  a11 = -( (acof*rhoscale)**2+16*muscale**2+lascale**2 )
                  a12 =  acof*bcof*rhoscale**2 + 4*muscale**2
                  a22 = -( muscale**2 + (bcof*rhoscale)**2 )
                  idet = 1/(a11*a22-a12*a12)
                  k1 = idet*(-a22*c1+a12*c2)
                  k2 = idet*( a12*c1-a11*c2)
                  cflcorr = 1
                  lowvscorr = 1
               endif

*** Was it successful ?
               if( mu(i,j,k).le.0 .or. lambda(i,j,k).lt.0 .or.
     *                                    rho(i,j,k).le.0 )then
                  info = -1
                  return
               endif
            enddo
         enddo
      enddo

      info = 0
      if( zcorr.eq.1 )then
         info = 1
      endif
      if( cflcorr.eq.1 )then
         info = info + 10
      endif
      if( lowvscorr.eq.1 )then
         info = info + 100
      endif
      
      end         

c-----------------------------------------------------------------------
      subroutine PROJECTMTRLC( ib,  ie,  jb,  je,  kb,  ke, 
     *                         iba, iea, jba, jea, kba, kea, 
     *                         rho, mu, lambda, dt,
     *                         met, jac, cfl, vsmin, rhoscale, muscale,
     *                         lascale, info )

***********************************************************************
***
*** Same as PROJECTMTRL, but for the curvilinear grid.
***
***********************************************************************

      implicit none
      integer ib, ie, jb, je, kb, ke, iba, iea, jba, jea, kba, kea
      integer i, j, k, info, zcorr, cflcorr, lowvscorr
      real*8 rho(ib:ie,jb:je,kb:ke), mu(ib:ie,jb:je,kb:ke)
      real*8 lambda(ib:ie,jb:je,kb:ke), jac(ib:ie,jb:je,kb:ke)
      real*8 met(4,ib:ie,jb:je,kb:ke)
      real*8 acof, bcof, c1, c2, acofbase, acof1, acof2, as2, d
      real*8 dt, cfl, rhoscale, muscale, lascale
      real*8 vsmin, k1, k2, idet, a11, a12, a22

      acofbase = (cfl/dt)**2
      if( vsmin .gt. 0 )then
         bcof = vsmin**2
      else
         bcof = -1
      endif
      zcorr     = 0
      cflcorr   = 0
      lowvscorr = 0
      do k=kba,kea
         do j=jba,jea
            do i=iba,iea
               if( mu(i,j,k).le.0 )then
                  mu(i,j,k) = 1e-6*muscale
                  zcorr = 1
               endif
               if( rho(i,j,k).le.0 )then
                  rho(i,j,k) = 1e-6*rhoscale
                  zcorr = 1
               endif
               if( lambda(i,j,k).lt.0 )then
                  lambda(i,j,k) = 0
                  zcorr = 1
               endif
               as2 = met(2,i,j,k)**2+met(3,i,j,k)**2+met(4,i,j,k)**2
               d = SQRT( (met(1,i,j,k)**2+as2)**2 -
     *                     4*met(1,i,j,k)**2*met(4,i,j,k)**2 )
               if( as2+d .gt. met(1,i,j,k)**2 )then
                  acof1 = 0.5d0*(5*met(1,i,j,k)**2+3*as2+d)
                  acof2 = 0.5d0*(  met(1,i,j,k)**2+  as2+d)
               else
                  acof1 = 3*met(1,i,j,k)**2+as2
                  acof2 =   met(1,i,j,k)**2
               endif
               acof = acofbase*jac(i,j,k)
               c1 = acof*rho(i,j,k)-acof1*mu(i,j,k)-acof2*lambda(i,j,k)
               c2 = mu(i,j,k)-bcof*rho(i,j,k)
               if( c1.lt.0 .and. c2.ge.0 )then
*** correct first condition
                  k1 = c1/( (acof*rhoscale)**2+(acof1*muscale)**2+
     *                                        (acof2*lascale)**2 )
                  mu(i,j,k)     = mu(i,j,k)    +muscale**2*acof1*k1
                  lambda(i,j,k) = lambda(i,j,k)+lascale**2*acof2*k1
                  rho(i,j,k) = rho(i,j,k) -rhoscale**2*acof*k1
                  c1=acof*rho(i,j,k)-acof1*mu(i,j,k)-acof2*lambda(i,j,k)
                  c2 = mu(i,j,k)-bcof*rho(i,j,k)
                  cflcorr = 1
               endif

               if( c1.ge.0 .and. c2.lt.0 )then
*** correct second condition
                  k2 = c2/(muscale**2+(bcof*rhoscale)**2)
                  mu(i,j,k)  = mu(i,j,k)  - muscale**2*k2
                  rho(i,j,k) = rho(i,j,k) + rhoscale**2*bcof*k2
                  c1=acof*rho(i,j,k)-acof1*mu(i,j,k)-acof2*lambda(i,j,k)
                  c2 = mu(i,j,k)-bcof*rho(i,j,k)
                  lowvscorr = 1
               endif

               if( c1.lt.0 .or. c2.lt.0 )then
*** correct both conditions
                  a11 = -( (acof*rhoscale)**2+(acof1*muscale)**2+
     *                                            (acof2*lascale)**2 )
                  a12 =  acof*bcof*rhoscale**2 + acof1*muscale**2
                  a22 = -( muscale**2 + (bcof*rhoscale)**2 )
                  idet = 1/(a11*a22-a12*a12)
                  k1 = idet*(-a22*c1+a12*c2)
                  k2 = idet*( a12*c1-a11*c2)
                  cflcorr   = 1
                  lowvscorr = 1
               endif

*** Was it successful ?
               if( mu(i,j,k).le.0 .or. lambda(i,j,k).lt.0 .or.
     *                                    rho(i,j,k).le.0 )then
                  info = -1
                  return
               endif
            enddo
         enddo
      enddo

      info = 0
      if( zcorr.eq.1 )then
         info = 1
      endif
      if( cflcorr.eq.1 )then
         info = info + 10
      endif
      if( lowvscorr.eq.1 )then
         info = info + 100
      endif
      
      end         

