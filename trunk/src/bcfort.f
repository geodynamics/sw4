c-----------------------------------------------------------------------
      subroutine bcfort( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     wind, nx, ny, nz,
     +     u, h, bccnd, sbop, mu, la, t,
     *     bforce1, bforce2, bforce3, bforce4, bforce5, bforce6, 
     +     om, ph, cv )
      implicit none
      real*8 d4a, d4b
      parameter( d4a=2d0/3, d4b=-1d0/12 )
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz
      integer s, wind(6,6), i, j, k, bccnd(6), w, kl, qq
      real*8 x, y, z, h, sbop(0:4)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8 bforce(3,ifirst:ilast,jfirst:jlast)
c note that the numbering of bforce adds one from C (side goes from 1 in Fortran)
      real*8 bforce1(3,*),  bforce2(3,*)
      real*8 bforce3(3,*),  bforce4(3,*)
      real*8 bforce5(3,*),  bforce6(3,*)
      real*8 ux, vy, wx, wy, uz, vz, wz, t, om, ph, cv

c the boundary window 'wind' is now an input argument

c loop over all sides of the 3-D domain
      do s=1,6
*** dirichlet condition, bccnd=1
*** supergrid condition, bccnd=2
c now assigning the forcing arrays outside of this routine!
        if( bccnd(s).eq.1 .or. bccnd(s).eq.2)then

            qq=1
            if (s.eq.1) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce1 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce1(1,qq)
                    u(2,i,j,k) = bforce1(2,qq)
                    u(3,i,j,k) = bforce1(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.2) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce2 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce2(1,qq)
                    u(2,i,j,k) = bforce2(2,qq)
                    u(3,i,j,k) = bforce2(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.3) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce3 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce3(1,qq)
                    u(2,i,j,k) = bforce3(2,qq)
                    u(3,i,j,k) = bforce3(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.4) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce4 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce4(1,qq)
                    u(2,i,j,k) = bforce4(2,qq)
                    u(3,i,j,k) = bforce4(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.5) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce5 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce5(1,qq)
                    u(2,i,j,k) = bforce5(2,qq)
                    u(3,i,j,k) = bforce5(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.6) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce6 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce6(1,qq)
                    u(2,i,j,k) = bforce6(2,qq)
                    u(3,i,j,k) = bforce6(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            endif

          elseif( bccnd(s).eq.3 )then
*** Periodic condition, bccnd=3
            if (s.eq.1) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i+nx,j,k)
                    u(2,i,j,k) = u(2,i+nx,j,k)
                    u(3,i,j,k) = u(3,i+nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.2) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i-nx,j,k)
                    u(2,i,j,k) = u(2,i-nx,j,k)
                    u(3,i,j,k) = u(3,i-nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.3) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j+ny,k)
                    u(2,i,j,k) = u(2,i,j+ny,k)
                    u(3,i,j,k) = u(3,i,j+ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.4) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j-ny,k)
                    u(2,i,j,k) = u(2,i,j-ny,k)
                    u(3,i,j,k) = u(3,i,j-ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.5) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k+nz)
                    u(2,i,j,k) = u(2,i,j,k+nz)
                    u(3,i,j,k) = u(3,i,j,k+nz)
                  enddo
                enddo
              enddo
            elseif (s.eq.6) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k-nz)
                    u(2,i,j,k) = u(2,i,j,k-nz)
                    u(3,i,j,k) = u(3,i,j,k-nz)
                  enddo
                enddo
              enddo
            endif
          elseif( bccnd(s).eq.0 )then
*** Free surface condition, bccnd=0
            if( s.ne.5 .and. s.ne.6 )then
              write(*,*) 
     *             'ERROR: Free surface condition ',
     +             'not implemented for side ', s
              stop
            endif

c moved the assignment of bforce5/6 into its own routine

**** Do the free surface condition (kl is the direction)
          if( s.eq.5 )then
            k = 1
            kl= 1

            do j=jfirst+2,jlast-2
              do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
*** Compute wx, wy, ux, vy by centered differences along boundary
                wx = 
     *               d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *               d4b*(u(3,i+2,j,k)-u(3,i-2,j,k)) 
                ux = 
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))+
     *               d4b*(u(1,i+2,j,k)-u(1,i-2,j,k)) 

                wy = 
     *               d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *               d4b*(u(3,i,j+2,k)-u(3,i,j-2,k)) 
                vy = 
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))+
     *               d4b*(u(2,i,j+2,k)-u(2,i,j-2,k)) 

                uz = 0
                vz = 0
                wz = 0
c interior contribution to uz, vz, wz (kl is the direction)
                do w=1,4
                  uz = uz + sbop(w)*u(1,i,j,k+kl*(w-1))
                  vz = vz + sbop(w)*u(2,i,j,k+kl*(w-1))
                  wz = wz + sbop(w)*u(3,i,j,k+kl*(w-1))
                enddo
                u(1,i,j,k-kl)=(-uz-kl*wx+kl*h*bforce5(1,qq)/mu(i,j,k))
     *               /sbop(0)
                u(2,i,j,k-kl)=(-vz-kl*wy+kl*h*bforce5(2,qq)/mu(i,j,k))
     *               /sbop(0)
                u(3,i,j,k-kl)=(-wz+(-kl*la(i,j,k)*(ux+vy)+
     *               kl*h*bforce5(3,qq))/(2*mu(i,j,k)+la(i,j,k)))
     *               /sbop(0)
              enddo
            enddo
          else
c s=6
            k = nz
            kl= -1

            do j=jfirst+2,jlast-2
              do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
*** Compute wx, wy, ux, vy by centered differences along boundary
                wx = 
     *               d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *               d4b*(u(3,i+2,j,k)-u(3,i-2,j,k)) 
                ux = 
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))+
     *               d4b*(u(1,i+2,j,k)-u(1,i-2,j,k)) 

                wy = 
     *               d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *               d4b*(u(3,i,j+2,k)-u(3,i,j-2,k)) 
                vy =
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))+
     *               d4b*(u(2,i,j+2,k)-u(2,i,j-2,k)) 

                uz = 0
                vz = 0
                wz = 0
c interior contribution to uz, vz, wz (kl is the direction)
                do w=1,4
                  uz = uz + sbop(w)*u(1,i,j,k+kl*(w-1))
                  vz = vz + sbop(w)*u(2,i,j,k+kl*(w-1))
                  wz = wz + sbop(w)*u(3,i,j,k+kl*(w-1))
                enddo
                u(1,i,j,k-kl)=(-uz-kl*wx+kl*h*bforce6(1,qq)/mu(i,j,k))
     *               /sbop(0)
                u(2,i,j,k-kl)=(-vz-kl*wy+kl*h*bforce6(2,qq)/mu(i,j,k))
     *               /sbop(0)
                u(3,i,j,k-kl)=(-wz+(-kl*la(i,j,k)*(ux+vy)+
     *               kl*h*bforce6(3,qq))/(2*mu(i,j,k)+la(i,j,k)))
     *               /sbop(0)
              enddo
            enddo
          endif

        endif
      enddo
      end


c-----------------------------------------------------------------------
      subroutine bcfortsg( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     wind, nx, ny, nz,
     +     u, h, bccnd, sbop, mu, la, t,
     *     bforce1, bforce2, bforce3, bforce4, bforce5, bforce6, 
     +     om, ph, cv, strx, stry )
      implicit none
      real*8 d4a, d4b
      parameter( d4a=2d0/3, d4b=-1d0/12 )
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz
      integer s, wind(6,6), i, j, k, bccnd(6), w, kl, qq
      real*8 x, y, z, h, sbop(0:4)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8 bforce(3,ifirst:ilast,jfirst:jlast)
c note that the numbering of bforce adds one from C (side goes from 1 in Fortran)
      real*8 bforce1(3,*),  bforce2(3,*)
      real*8 bforce3(3,*),  bforce4(3,*)
      real*8 bforce5(3,*),  bforce6(3,*)
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      real*8 ux, vy, wx, wy, uz, vz, wz, t, om, ph, cv

c the boundary window 'wind' is now an input argument

c loop over all sides of the 3-D domain
      do s=1,6
*** dirichlet condition, bccnd=1
*** supergrid condition, bccnd=2
c now assigning the forcing arrays outside of this routine!
        if( bccnd(s).eq.1 .or. bccnd(s).eq.2)then

            qq=1
            if (s.eq.1) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce1 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce1(1,qq)
                    u(2,i,j,k) = bforce1(2,qq)
                    u(3,i,j,k) = bforce1(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.2) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce2 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce2(1,qq)
                    u(2,i,j,k) = bforce2(2,qq)
                    u(3,i,j,k) = bforce2(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.3) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce3 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce3(1,qq)
                    u(2,i,j,k) = bforce3(2,qq)
                    u(3,i,j,k) = bforce3(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.4) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce4 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce4(1,qq)
                    u(2,i,j,k) = bforce4(2,qq)
                    u(3,i,j,k) = bforce4(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.5) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce5 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce5(1,qq)
                    u(2,i,j,k) = bforce5(2,qq)
                    u(3,i,j,k) = bforce5(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.6) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce6 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce6(1,qq)
                    u(2,i,j,k) = bforce6(2,qq)
                    u(3,i,j,k) = bforce6(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            endif

          elseif( bccnd(s).eq.3 )then
*** Periodic condition, bccnd=3
            if (s.eq.1) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i+nx,j,k)
                    u(2,i,j,k) = u(2,i+nx,j,k)
                    u(3,i,j,k) = u(3,i+nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.2) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i-nx,j,k)
                    u(2,i,j,k) = u(2,i-nx,j,k)
                    u(3,i,j,k) = u(3,i-nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.3) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j+ny,k)
                    u(2,i,j,k) = u(2,i,j+ny,k)
                    u(3,i,j,k) = u(3,i,j+ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.4) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j-ny,k)
                    u(2,i,j,k) = u(2,i,j-ny,k)
                    u(3,i,j,k) = u(3,i,j-ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.5) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k+nz)
                    u(2,i,j,k) = u(2,i,j,k+nz)
                    u(3,i,j,k) = u(3,i,j,k+nz)
                  enddo
                enddo
              enddo
            elseif (s.eq.6) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k-nz)
                    u(2,i,j,k) = u(2,i,j,k-nz)
                    u(3,i,j,k) = u(3,i,j,k-nz)
                  enddo
                enddo
              enddo
            endif
          elseif( bccnd(s).eq.0 )then
*** Free surface condition, bccnd=0
            if( s.ne.5 .and. s.ne.6 )then
              write(*,*) 
     *             'ERROR: Free surface condition ',
     +             'not implemented for side ', s
              stop
            endif

c moved the assignment of bforce5/6 into its own routine

**** Do the free surface condition (kl is the direction)
          if( s.eq.5 )then
            k = 1
            kl= 1

            do j=jfirst+2,jlast-2
              do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
*** Compute wx, wy, ux, vy by centered differences along boundary
                wx = strx(i)*(
     *               d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *               d4b*(u(3,i+2,j,k)-u(3,i-2,j,k)) )
                ux = strx(i)*(
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))+
     *               d4b*(u(1,i+2,j,k)-u(1,i-2,j,k)) )

                wy = stry(j)*(
     *               d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *               d4b*(u(3,i,j+2,k)-u(3,i,j-2,k)) )
                vy = stry(j)*(
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))+
     *               d4b*(u(2,i,j+2,k)-u(2,i,j-2,k)) )

                uz = 0
                vz = 0
                wz = 0
c interior contribution to uz, vz, wz (kl is the direction)
                do w=1,4
                  uz = uz + sbop(w)*u(1,i,j,k+kl*(w-1))
                  vz = vz + sbop(w)*u(2,i,j,k+kl*(w-1))
                  wz = wz + sbop(w)*u(3,i,j,k+kl*(w-1))
                enddo
                u(1,i,j,k-kl)=(-uz-kl*wx+kl*h*bforce5(1,qq)/mu(i,j,k))
     *               /sbop(0)
                u(2,i,j,k-kl)=(-vz-kl*wy+kl*h*bforce5(2,qq)/mu(i,j,k))
     *               /sbop(0)
                u(3,i,j,k-kl)=(-wz+(-kl*la(i,j,k)*(ux+vy)+
     *               kl*h*bforce5(3,qq))/(2*mu(i,j,k)+la(i,j,k)))
     *               /sbop(0)
              enddo
            enddo
          else
c s=6
            k = nz
            kl= -1

            do j=jfirst+2,jlast-2
              do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
*** Compute wx, wy, ux, vy by centered differences along boundary
                wx = strx(i)*(
     *               d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *               d4b*(u(3,i+2,j,k)-u(3,i-2,j,k)) )
                ux = strx(i)*(
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))+
     *               d4b*(u(1,i+2,j,k)-u(1,i-2,j,k)) )

                wy = stry(j)*(
     *               d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *               d4b*(u(3,i,j+2,k)-u(3,i,j-2,k)) )
                vy = stry(j)*(
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))+
     *               d4b*(u(2,i,j+2,k)-u(2,i,j-2,k)) )

                uz = 0
                vz = 0
                wz = 0
c interior contribution to uz, vz, wz (kl is the direction)
                do w=1,4
                  uz = uz + sbop(w)*u(1,i,j,k+kl*(w-1))
                  vz = vz + sbop(w)*u(2,i,j,k+kl*(w-1))
                  wz = wz + sbop(w)*u(3,i,j,k+kl*(w-1))
                enddo
                u(1,i,j,k-kl)=(-uz-kl*wx+kl*h*bforce6(1,qq)/mu(i,j,k))
     *               /sbop(0)
                u(2,i,j,k-kl)=(-vz-kl*wy+kl*h*bforce6(2,qq)/mu(i,j,k))
     *               /sbop(0)
                u(3,i,j,k-kl)=(-wz+(-kl*la(i,j,k)*(ux+vy)+
     *               kl*h*bforce6(3,qq))/(2*mu(i,j,k)+la(i,j,k)))
     *               /sbop(0)
              enddo
            enddo
          endif

        endif
      enddo
      end


c----------------------------------------------------------------------
      subroutine TWDIRBDRY( wind, h, t, om, cv, ph, bforce )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, wind(6)
      real*8 bforce(3,*), h, t, om, cv, ph, x, y, z
      integer i, j, k, qq
c
c NOTE: pass in the window for one side, i.e., wind(1,side) in the calling routine
c
*** Twilight forced dirichlet condition
      qq = 1
      do k=wind(5),wind(6)
c need to add zmin to work in a composite grid setting
        z = (k-1)*h
        do j=wind(3),wind(4)
          y = (j-1)*h
          do i=wind(1),wind(2)
            x = (i-1)*h
            bforce(1,qq) = sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph)
            bforce(2,qq) = sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph)
            bforce(3,qq) = sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t))
            qq = qq+1
          enddo
        enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine TWFRSURFZ( ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     ni, nj, nk, h, kz, t, omega, c, phase, bforce,
     *                      mu, lambda )
c     *        momega,   mphase, ampmu, amplambda, fo )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 bforce(3,ifirst:ilast,jfirst:jlast), h
      integer ni, nj, nk, i, j, kz
      doubleprecision mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision x
      doubleprecision y
      doubleprecision z
      doubleprecision t
      doubleprecision omega
      doubleprecision c
      doubleprecision phase
c      doubleprecision momega
c      doubleprecision mphase
c      doubleprecision ampmu
c      doubleprecision amplambda

      doubleprecision forces(3)
c      doubleprecision t10
      doubleprecision t13
      doubleprecision t15
      doubleprecision t16
      doubleprecision t19
c      doubleprecision t2
      doubleprecision t20
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t28
      doubleprecision t29
c      doubleprecision t3
      doubleprecision t32
      doubleprecision t33
      doubleprecision t34
      doubleprecision t37
      doubleprecision t38
      doubleprecision t43
      doubleprecision t44
      doubleprecision t49
c      doubleprecision t54
c      doubleprecision t56
c      doubleprecision t6
      doubleprecision t60
      doubleprecision t62
      doubleprecision t65
c      doubleprecision t9

      z = (kz-1)*h
c the do loops should span jfirst,jlast and ifirst,ilast
c      do j=1,nj
      do j=jfirst,jlast
         y = (j-1)*h
c         do i=1,ni
         do i=ifirst,ilast
            x=(i-1)*h
c        t2 = momega*x+mphase
c        t3 = cos(t2)
c        t6 = sin(momega*y+mphase)
c        t9 = momega*z+mphase
c        t10 = sin(t9)
        t13 = mu(i,j,kz)
        t15 = omega*x+phase
        t16 = cos(t15)
        t19 = omega*y+phase
        t20 = sin(t19)
        t21 = c*t
        t23 = omega*(z-t21)
        t24 = sin(t23)
        t28 = omega*(x-t21)
        t29 = sin(t28)
        t32 = omega*z+phase
        t33 = cos(t32)
        t34 = t33*omega
        forces(1) = t13*(t16*omega*t20*t24+t29*t20*t34)
        t37 = sin(t15)
        t38 = cos(t19)
        t43 = omega*(y-t21)
        t44 = sin(t43)
        forces(2) = t13*(t37*t38*omega*t24+t37*t44*t34)
        t49 = cos(t23)
c        t54 = sin(t2)
c        t56 = cos(t9)
        t60 = cos(t28)
        t62 = sin(t32)
        t65 = cos(t43)
        forces(3) = 2*t13*t37*t20*t49*omega+lambda(i,j,kz)*(t6
     #0*omega*t20*t62+t37*t65*omega*t62+t37*t20*t49*omega)
        bforce(1,i,j) = forces(1)
        bforce(2,i,j) = forces(2)
        bforce(3,i,j) = forces(3)
      enddo
      enddo
      end


c-----------------------------------------------------------------------
      subroutine TWFRSURFZSGSTR( ifirst, ilast, jfirst, jlast, kfirst, 
     *                  klast, h, kz, t, om, c, ph, omstrx, omstry,
     *                  bforce, mu, lambda )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 bforce(3,ifirst:ilast,jfirst:jlast), h
      integer i, j, kz
      doubleprecision mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision x
      doubleprecision y
      doubleprecision z
      doubleprecision t
      doubleprecision om
      doubleprecision c
      doubleprecision ph
      doubleprecision omstrx
      doubleprecision omstry

      doubleprecision forces(3)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t11
      doubleprecision t12
      doubleprecision t15
      doubleprecision t17
      doubleprecision t19
      doubleprecision t20
      doubleprecision t22
      doubleprecision t24
      doubleprecision t25
      doubleprecision t29
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t36
      doubleprecision t39
      doubleprecision t4
      doubleprecision t40
      doubleprecision t46
      doubleprecision t51
      doubleprecision t53
      doubleprecision t56
      doubleprecision t6
      doubleprecision t7

      z = (kz-1)*h
      do j=jfirst,jlast
         y = (j-1)*h
         do i=ifirst,ilast
            x=(i-1)*h
            t1 = c*t
            t3 = om*(x-t1)
            t4 = sin(t3)
            t6 = om*y+ph
            t7 = sin(t6)
            t10 = om*z+ph
            t11 = cos(t10)
            t12 = t11*om
            t15 = sin(omstrx*x)
            t17 = 1+t15/2
            t19 = om*x+ph
            t20 = cos(t19)
            t22 = om*t7
            t24 = om*(z-t1)
            t25 = sin(t24)
            forces(1) = mu(i,j,kz)*(t4*t7*t12+t17*t20*t22*t25)
            t29 = sin(t19)
            t31 = om*(y-t1)
            t32 = sin(t31)
            t36 = sin(omstry*y)
            t39 = (1+t36/2)*t29
            t40 = cos(t6)
            forces(2) = mu(i,j,kz)*(t29*t32*t12+t39*t40*om*t25)
            t46 = cos(t24)
            t51 = cos(t3)
            t53 = sin(t10)
            t56 = cos(t31)
            forces(3) = 2*mu(i,j,kz)*t29*t7*t46*om+lambda(i,j,kz)*
     #(t17*t51*t22*t53+t39*t56*om*t53+t29*t7*t46*om)
            bforce(1,i,j) = forces(1)
            bforce(2,i,j) = forces(2)
            bforce(3,i,j) = forces(3)
         enddo
      enddo
      return
      end

