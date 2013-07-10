subroutine timestep(b, lwork)
! calculate the pseudo time step
! ==================================================================
  use typelib4
  use typepara
  use ctrlib
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), intent(inout) :: b
  logical, intent(in) :: lwork

  integer  :: i, j, k 
  integer  :: il,jl,kl  !nx,ny,nz,
  integer  :: ie,je,ke,ib,jb,kb
  real(dp) :: cc, lx, ly, lz, uu, uv, uw
  real(dp) :: sx, sy, sz, qs, ds, cs
  real(dp) :: fict, dtime, cflrat, ex, reval
  real(dp) :: qq, eqq
  real(dp) :: plim, dpi, dpj, dpk, r
  real(dp) :: cfac, radvi, radvj, radvk
  real(dp) :: fmue, f1, f2, fac, dtv
! ******************************************************************

  if (.not. lwork) return

  il = b%il
  jl = b%jl
  kl = b%kl
  ie = b%ie
  je = b%je
  ke = b%ke
  ib = b%ib
  jb = b%jb
  kb = b%kb

  b%dtim = 10.0d0
  plim   = 0.001d0

  if (eqntype == 'e' .or. eqntype == 'E') then

    do k=2,kl
      do j=2,jl
        do i=2,il
          cc = gamma*b%p(i,j,k)/b%w(i,j,k,1)
          if(cc < 0.0d0 .or. cc .ne. cc) then
             write(*,'(3I4,3F10.5," Pressure ",F10.5," Density ",5F10.5)') &
                  i,j,k,b%x(i,j,k,:),b%p(i,j,k),b%w(i,j,k,:)
             stop
          end if
          if(abs(b%w(i,j,k,1)) < 1.0e-8) then
            print *, "Zero density in step : ",i,j,k,b%w(i,j,k,1)
          end if
          cc = sqrt(max(cc,0.0d0))
          uu = b%w(i,j,k,2)/b%w(i,j,k,1)
          uv = b%w(i,j,k,3)/b%w(i,j,k,1)
          uw = b%w(i,j,k,4)/b%w(i,j,k,1)

          ! calculate spectral radii
          sx = 0.5d0*(b%si(i,j,k,1) + b%si(i+1,j,k,1))
          sy = 0.5d0*(b%si(i,j,k,2) + b%si(i+1,j,k,2))
          sz = 0.5d0*(b%si(i,j,k,3) + b%si(i+1,j,k,3))          
          qs = abs(uu*sx +uv*sy +uw*sz)
          ds = sx*sx +sy*sy +sz*sz
          cs = cc*sqrt(ds)
          b%radi(i,j,k) = qs +cs
          
          sx = 0.5d0*(b%sj(i,j,k,1) +b%sj(i,j+1,k,1))
          sy = 0.5d0*(b%sj(i,j,k,2) +b%sj(i,j+1,k,2))
          sz = 0.5d0*(b%sj(i,j,k,3) +b%sj(i,j+1,k,3))
          qs = abs(uu*sx +uv*sy +uw*sz)
          ds = sx*sx +sy*sy +sz*sz
          cs = cc*sqrt(ds)
          b%radj(i,j,k) = qs +cs

          sx = 0.5d0*(b%sk(i,j,k,1) +b%sk(i,j,k+1,1))
          sy = 0.5d0*(b%sk(i,j,k,2) +b%sk(i,j,k+1,2))
          sz = 0.5d0*(b%sk(i,j,k,3) +b%sk(i,j,k+1,3))
          qs = abs(uu*sx +uv*sy +uw*sz)
          ds = sx*sx +sy*sy +sz*sz
          cs = cc*sqrt(ds)
          b%radk(i,j,k) = qs +cs          

          ! --- CALC TIME STEP

          fict  = b%cellvolume/(b%radi(i,j,k) + b%radj(i,j,k) + b%radk(i,j,k))
          dtime = 2.0d0*dtstep/(3.0d0*cfl)

          if (execmode==0) then
            !
            ! TIME STEP FOR STEADY-STATE FLOW
            !
            b%vdtim(i,j,k) = fict

          else

            !if(dtime < fict) print *, fict, dtime,dtstep
            b%vdtim(i,j,k) = min(fict,dtime)

          end if
        
          ! adaptive time step
          dpi    = abs((b%p(i+1,j,k)  -2.0d0*b%p(i,j,k)  +b%p(i-1,j,k))/            &
                      (b%p(i+1,j,k)  +2.0d0*b%p(i,j,k)  +b%p(i-1,j,k)  +plim))
          dpj    = abs((b%p(i,j+1,k)  -2.0d0*b%p(i,j,k)  +b%p(i,j-1,k))/            &
                      (b%p(i,j+1,k)  +2.0d0*b%p(i,j,k)  +b%p(i,j-1,k)  +plim))
          dpk    = abs((b%p(i,j,k+1)  -2.0d0*b%p(i,j,k)  +b%p(i,j,k-1))/            &
                      (b%p(i,j,k+1)  +2.0d0*b%p(i,j,k)  +b%p(i,j,k-1)  +plim))
          r      = 1.0d0/(1.0d0+min(dim(cfl,1.0d0),2.0d0*(dpi +dpj +dpk*(1.0d0-i2d))))
          b%vdtim(i,j,k) = b%vdtim(i,j,k)*r

          ! global time step
          if(b%vdtim(i,j,k)<b%dtim) b%dtim = b%vdtim(i,j,k)
        end do
      end do
    end do

  else if(eqntype == 'n' .or. eqntype == 'N') then

    if (distype.eq.'c' .or. distype.eq.'C') then
      cfac = 4.0d0
    else if(distype.eq.'r2' .or. distype.eq.'R2') then 
      cfac = 1.0d0
    else
      cfac = 2.0d0
    end if

    do k=2,kl
      do j=2,jl
        do i=2,il
          cc = gamma*b%p(i,j,k)/b%w(i,j,k,1)
          if(cc < 0.0d0 .or. cc .ne. cc) then
             write(*,'(3I4,3F10.5," Pressure ",F10.5," Density ",5F10.5)') &
                  i,j,k,b%x(i,j,k,:),b%p(i,j,k),b%w(i,j,k,:)
             stop
          end if
          cc = sqrt(max(cc,0.0d0))
          uu = b%w(i,j,k,2)/b%w(i,j,k,1)
          uv = b%w(i,j,k,3)/b%w(i,j,k,1)
          uw = b%w(i,j,k,4)/b%w(i,j,k,1)

          fmue = viscl/prlmean                !+ b%vt(i,j,k)/prtmean    ! ?
          f1   = 4.0d0/b%w(i,j,k,1)/3.0d0
          f2   = gamma/b%w(i,j,k,1)
          fac  = max(f1,f2)
          dtv  = fac*fmue/b%cellvolume

          ! calculate spectral radii
          sx = 0.5d0*(b%si(i,j,k,1) + b%si(i+1,j,k,1))
          sy = 0.5d0*(b%si(i,j,k,2) + b%si(i+1,j,k,2))
          sz = 0.5d0*(b%si(i,j,k,3) + b%si(i+1,j,k,3))          
          qs = abs(uu*sx +uv*sy +uw*sz)
          ds = sx*sx +sy*sy +sz*sz
          cs = cc*sqrt(ds)
          b%radi(i,j,k) = qs +cs
          radvi         = dtv*ds
          
          sx = 0.5d0*(b%sj(i,j,k,1) +b%sj(i,j+1,k,1))
          sy = 0.5d0*(b%sj(i,j,k,2) +b%sj(i,j+1,k,2))
          sz = 0.5d0*(b%sj(i,j,k,3) +b%sj(i,j+1,k,3))
          qs = abs(uu*sx +uv*sy +uw*sz)
          ds = sx*sx +sy*sy +sz*sz
          cs = cc*sqrt(ds)
          b%radj(i,j,k) = qs +cs
          radvj         = dtv*ds

          sx = 0.5d0*(b%sk(i,j,k,1) +b%sk(i,j,k+1,1))
          sy = 0.5d0*(b%sk(i,j,k,2) +b%sk(i,j,k+1,2))
          sz = 0.5d0*(b%sk(i,j,k,3) +b%sk(i,j,k+1,3))
          qs = abs(uu*sx +uv*sy +uw*sz)
          ds = sx*sx +sy*sy +sz*sz
          cs = cc*sqrt(ds)
          b%radk(i,j,k) = qs +cs          
          radvk         = dtv*ds
        
          ! --- CALC TIME STEP

          fict  = b%cellvolume/(b%radi(i,j,k) + b%radj(i,j,k) + b%radk(i,j,k) &
                  + cfac*(radvi + radvj + radvk))

          dtime = 2.0d0*dtstep/(3.0d0*cfl)

          if (execmode==0) then
            !
            !
            !
            b%vdtim(i,j,k) = fict

          else

            b%vdtim(i,j,k) = min(fict,dtime)

          end if

          ! adaptive time step
          dpi    = abs((b%p(i+1,j,k)  -2.0d0*b%p(i,j,k)  +b%p(i-1,j,k))/            &
                      (b%p(i+1,j,k)  +2.0d0*b%p(i,j,k)  +b%p(i-1,j,k)  +plim))
          dpj    = abs((b%p(i,j+1,k)  -2.0d0*b%p(i,j,k)  +b%p(i,j-1,k))/            &
                      (b%p(i,j+1,k)  +2.0d0*b%p(i,j,k)  +b%p(i,j-1,k)  +plim))
          dpk    = abs((b%p(i,j,k+1)  -2.0d0*b%p(i,j,k)  +b%p(i,j,k-1))/            &
                      (b%p(i,j,k+1)  +2.0d0*b%p(i,j,k)  +b%p(i,j,k-1)  +plim))
          r      = 1.0d0/(1.0d0+min(dim(cfl,1.0d0),2.0d0*(dpi +dpj +dpk*(1.0d0-i2d))))
          b%vdtim(i,j,k) = b%vdtim(i,j,k)*r

          ! global time step
          if(b%vdtim(i,j,k)<b%dtim) b%dtim = b%vdtim(i,j,k)

        end do
      end do
    end do
  end if

  do k=0, kb
    do j=0, jb
      do i=0, ib
        if (b%blank(i,j,k) == 0) then
          b%vdtim(i,j,k) = 0.0d0
          b%radi(i,j,k)  = 0.0d0
          b%radj(i,j,k)  = 0.0d0
          b%radk(i,j,k)  = 0.0d0
        end if
      end do
    end do
  end do


! --- set values of spectral radii in dummy nodes (needed by "dflux") --------
  do i=0,ib
     do j=0,jb
        b%radi(i,j, 0) = b%radi(i,j, 2)
        b%radi(i,j, 1) = b%radi(i,j, 2)
        b%radi(i,j,ke) = b%radi(i,j,kl)
        b%radi(i,j,kb) = b%radi(i,j,kl)

        b%radj(i,j, 0) = b%radj(i,j, 2)
        b%radj(i,j, 1) = b%radj(i,j, 2)
        b%radj(i,j,ke) = b%radj(i,j,kl)
        b%radj(i,j,kb) = b%radj(i,j,kl)

        b%radk(i,j, 0) = b%radk(i,j, 2)
        b%radk(i,j, 1) = b%radk(i,j, 2)
        b%radk(i,j,ke) = b%radk(i,j,kl)
        b%radk(i,j,kb) = b%radk(i,j,kl)
     end do
  end do

  do j=2,jl
     do k=0,kb
        b%radi(0 ,j,k) = b%radi(2 ,j,k)
        b%radi(1 ,j,k) = b%radi(2 ,j,k)
        b%radi(ie,j,k) = b%radi(il,j,k)
        b%radi(ib,j,k) = b%radi(il,j,k)

        b%radj(0 ,j,k) = b%radj(2 ,j,k)
        b%radj(1 ,j,k) = b%radj(2 ,j,k)
        b%radj(ie,j,k) = b%radj(il,j,k)
        b%radj(ib,j,k) = b%radj(il,j,k)

        b%radk(0 ,j,k) = b%radk(2 ,j,k)
        b%radk(1 ,j,k) = b%radk(2 ,j,k)
        b%radk(ie,j,k) = b%radk(il,j,k)
        b%radk(ib,j,k) = b%radk(il,j,k)
     end do
  end do

  do i=2,il
     do k=2,kl
        b%radi(i,0 ,k) = b%radi(i,2, k)
        b%radi(i,1 ,k) = b%radi(i,2, k)
        b%radi(i,je,k) = b%radi(i,jl,k)
        b%radi(i,jb,k) = b%radi(i,jl,k)

        b%radj(i,0 ,k) = b%radj(i,2, k)
        b%radj(i,1 ,k) = b%radj(i,2, k)
        b%radj(i,je,k) = b%radj(i,jl,k)
        b%radj(i,jb,k) = b%radj(i,jl,k)

        b%radk(i,0 ,k) = b%radk(i,2, k)
        b%radk(i,1 ,k) = b%radk(i,2, k)
        b%radk(i,je,k) = b%radk(i,jl,k)
        b%radk(i,jb,k) = b%radk(i,jl,k)
     end do
  end do

! --- set coeff. of implicit residual smoothing -------------------------------

  if ( smoo .gt. 0.0d0) then
     cflrat = sqrt(1.d0+4.d0*smoo)
     do k=2,kl
        do j=2,jl
           do i=2,il
              reval  = b%radj(i,j,k)/b%radi(i,j,k)+b%radk(i,j,k)/b%radi(i,j,k)
              ex     = 0.25d0*((cflrat/(1.d0+0.0625d0*reval ))**2.0d0 - 1.0d0)
              ex     = min( smoo, ex )
              b%epsx(i,j,k) = max( 0.0d0, ex )

              reval  = b%radk(i,j,k)/b%radj(i,j,k)+b%radi(i,j,k)/b%radj(i,j,k)
              ex     = 0.25d0*((cflrat/(1.d0+0.0625d0*reval ))**2.0d0 - 1.0d0)
              ex     = min( smoo, ex )
              b%epsy(i,j,k) = max( 0.0d0, ex )

              reval  = b%radi(i,j,k)/b%radk(i,j,k)+b%radj(i,j,k)/b%radk(i,j,k)
              ex     = 0.25d0*((cflrat/(1.d0+0.0625d0*reval ))**2.0d0 - 1.0d0)
              ex     = min( smoo, ex )
              b%epsz(i,j,k) = max( 0.0d0, ex )
              !print *, i,j,k,b%epsx(i,j,k),b%epsy(i,j,k),b%epsz(i,j,k)
            enddo
        end do
     end do
  end if


  return
end subroutine timestep
