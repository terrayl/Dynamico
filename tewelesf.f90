subroutine tewelesf(vart, varf, plat, plon, ntt, ntf, nlat, nlon, tws)
! Estimate the Teweles-Wobus distance between the given atmospheric circulation 
! variable, say X, daily or monthly mean and all other daily or monthly means within the training period
! 
! INPUT->  vart(nxy,ntt)      : X array for all months/days of the training period
!          varf(nxy,ntf)      : X array for the month/day under scrutiny for the analogue search
!          plat(nlat)         : latitude array
!          plon(nlon)         : longitude array
!          ntt                : number of months/days in the training period 
!          ntf                : number of months/days under scrutiny for the analogue search
!          nlat               : number of latitudes
!          nlon               : number of longitudes
!
! OUTPUT-> tws(ntf,ntt)       : distance array
! 
! Note: this function uses a simple cos(lat) weighting in the gradient estimation
  implicit none
   integer      ntt, ntf, nlat, nlon
   integer      n, m, i, j
   real         vart(nlon,nlat,ntt), varf(nlon,nlat,ntf), plat(nlat), plon(nlon), tws(ntf,ntt)
   real         di(nlat), dvari(nlat), mdvari(nlat), dj(nlon), dvarj(nlon), mdvarj(nlon)
   real         gradi, gradj, Gi, Gj, pi
   
   ! Get Pi value.
   pi             = 4.*atan(1.)

   ! Weights (works for regular grids)
   ! 
   di             = cos(plat*pi/180) * 2 * (plon(2) - plon(1))
   dj             = 2 * abs(plat(2) - plat(1))

   !$OMP PARALLEL DO PRIVATE(gradi, gradj, Gi, Gj, dvari, dvarj, mdvari, mdvarj, m, i, j)
   do n=1,ntt
      do m=1,ntf

         gradi             = 0
         Gi                = 0

         do i=2,nlon-1
            ! Gradient différences in longitude
            dvari          = abs(     (vart(i+1,:,n) - vart(i-1,:,n))     - (varf(i+1,:,m) - varf(i-1,:,m))   )
            ! Weighting
            gradi          = gradi + sum(dvari / di)

            ! Maximum gradient
            mdvari         = max( abs( vart(i+1,:,n) - vart(i-1,:,n) ), abs( varf(i+1,:,m) - varf(i-1,:,m) )  )
            Gi             = Gi + sum(mdvari / di)
         end do


         gradj             = 0
         Gj                = 0

         do j=2,nlat-1
            !  Gradient différences in latitude
            dvarj          = abs(     (vart(:,j+1,n) - vart(:,j-1,n))     - (varf(:,j+1,m) - varf(:,j-1,m))   )
            ! Weighting
            gradj          = gradj + sum(dvarj / dj)

            ! Maximum gradient
            mdvarj         = max( abs( vart(:,j+1,n) - vart(:,j-1,n) ), abs( varf(:,j+1,m) - varf(:,j-1,m) )  )
            Gj             = Gj + sum(mdvarj / dj)
         end do
            ! Get distance
        tws(m,n)             = 100 * (gradi + gradj) / (Gi + Gj)

      end do

   end do
   !$OMP END PARALLEL DO

end subroutine tewelesf
