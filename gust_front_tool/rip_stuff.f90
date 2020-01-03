module rip_stuff_mod

use types_mod, only : wrf_met_data, derived_data, wrf_static_data
use namelist_and_constants, only : pi, missing_val

implicit none

private

! public subroutines
public :: rip_potential_temperature, rip_frontogensis, rip_gradient

!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!

subroutine rip_potential_temperature(grid)

   type(wrf_met_data), intent(inout) :: grid
   integer :: nx,ny
   real, dimension(:,:), allocatable :: gammam
   real :: rgas, rgasmd, cp, cpmd, gamma, gammamd , grav ! not all these are used.

   nx = grid%nx
   ny = grid%ny

   ! rip constants
   rgas=287.04  !J/K/kg
   rgasmd=.608   ! rgas_moist=rgas*(1.+rgasmd*qvp)
   cp=1004.     ! J/K/kg  Note: not using Bolton's value of 1005.7
   cpmd=.887   ! cp_moist=cp*(1.+cpmd*qvp)
   gamma=rgas/cp
   gammamd=rgasmd-cpmd  ! gamma_moist=gamma*(1.+gammamd*qvp)
   grav=9.81    ! m/s**2

   allocate( gammam(nx,ny))
   gammam = gamma*(1.+gammamd*grid%q2)
   grid%theta2 = grid%t2*(100000./grid%psfc)**gammam
   write(*,'(a,2(f10.3,1x))')'Min/max input vals of theta2',minval(grid%theta2),maxval(grid%theta2)
   deallocate(gammam)

end subroutine rip_potential_temperature

subroutine rip_derivative(grid_static,arr,dir,deriv) ! derivative on mass points

   type(wrf_static_data), intent(in)  :: grid_static
   real, dimension(grid_static%nx,grid_static%ny), intent(in) :: arr
   character(len=1), intent(in)       :: dir
   real, dimension(grid_static%nx,grid_static%ny), intent(inout) :: deriv

   integer :: i,j,dg,jp,jm, ip, im
   real    :: dsi, dssi

!  Modified from RIP/src/ddx.f
   dsi = 1./grid_static%dx  ! inverse horizonal grid spacing

   if ( dir.eq.'x') then     ! derivative WRT x
      do i=2,grid_static%nx-1  ! loop over x
         ip = i+1
         im = i-1
         dg=ip-im
         do j=2,grid_static%ny-1  ! loop over y
            dssi=grid_static%mapfac_u(i,j)*dsi/dg
            deriv(i,j)=dssi*(arr(ip,j)-arr(im,j))
         enddo
      enddo
   else if ( dir.eq.'y') then   ! derivative WRT y
      do j=2,grid_static%ny-1  ! loop over y
         jp = j+1
         jm = j-1
         dg=jp-jm
         do i=2,grid_static%nx-1  ! loop over x
            dssi=grid_static%mapfac_u(i,j)*dsi/dg
            deriv(i,j)=dssi*(arr(i,jp)-arr(i,jm))
         enddo
      enddo
   else
      write(*,*)'Error in rip_derivative. dir needs to be x or y. Exiting.'
      stop
   endif

   ! deal with boundaries...just set them equal to the adjacent row/column
   deriv(1,:)              = deriv(2,:)              ! West
   deriv(grid_static%nx,:) = deriv(grid_static%nx-1,:) ! East
   deriv(:,1)              = deriv(:,2)              ! South
   deriv(:,grid_static%ny) = deriv(:,grid_static%ny-1) ! North

end subroutine rip_derivative

subroutine rip_gradient(grid_static,field,grad_mag,grad_angle)

   type(wrf_static_data), intent(in) :: grid_static
   real, dimension(grid_static%nx,grid_static%ny), intent(in) :: field
   real, dimension(grid_static%nx,grid_static%ny), intent(inout) :: grad_mag
   real, dimension(grid_static%nx,grid_static%ny), intent(inout) :: grad_angle

   real, dimension(:,:), allocatable :: ddx, ddy

   ! initialize output to unphysical value
   grad_mag = missing_val

   allocate( ddx(grid_static%nx,grid_static%ny))
   allocate( ddy(grid_static%nx,grid_static%ny))

   ! get derivates at mass points
   call rip_derivative(grid_static,field,'x',ddx)
   call rip_derivative(grid_static,field,'y',ddy)

   ! get magnitude of gradient
   grad_mag = sqrt(ddx*ddx + ddy*ddy) 
   
   ! angle of the gradient
   grad_angle = atan2(ddy,ddx) ! returns arctan(ddy/ddx) in radians
   grad_angle = grad_angle * (180.0/pi) ! convert to degrees
   where( grad_angle .lt.0 ) grad_angle = grad_angle + 360.0 ! get positive values.
   
   deallocate(ddx,ddy)

end subroutine rip_gradient

subroutine rip_frontogensis(grid_static,grid,derived_fields)

   type(wrf_static_data), intent(in) :: grid_static
   type(wrf_met_data), intent(in)    :: grid
   type(derived_data), intent(inout) :: derived_fields

   real, dimension(:,:), allocatable :: dtheta_Dx, dtheta_Dy, grad_theta
   real, dimension(:,:), allocatable :: dudx, dudy, dvdx, dvdy
   real, dimension(:,:), allocatable :: pl3, pl4

   integer :: nxm, nym

   nxm = grid_static%nx
   nym = grid_static%ny

   ! allocate arrays
   allocate( dtheta_Dx(nxm,nym))
   allocate( dtheta_Dy(nxm,nym))
   allocate( grad_theta(nxm,nym))
   allocate( dudx(nxm,nym), dudy(nxm,nym))
   allocate( dvdx(nxm,nym), dvdy(nxm,nym))
   allocate( pl3(nxm,nym), pl4(nxm,nym))

   ! get derivates at mass points
   call rip_derivative(grid_static,grid%theta2,'x',dtheta_Dx)
   call rip_derivative(grid_static,grid%theta2,'y',dtheta_Dy)
   call rip_derivative(grid_static,grid%u10,'x',dudx)
   call rip_derivative(grid_static,grid%u10,'y',dudy)
   call rip_derivative(grid_static,grid%v10,'x',dvdx)
   call rip_derivative(grid_static,grid%v10,'y',dvdy)

   grad_theta =sqrt(dtheta_Dx*dtheta_Dx + dtheta_Dy*dtheta_Dy) ! magnitude of theta gradient

   ! load up output structure
   derived_fields%dtheta_Dx = dtheta_Dx * 100000. ! convert from K/m to K/100km
   derived_fields%dtheta_Dy = dtheta_Dy * 100000. 
   derived_fields%grad_theta = grad_theta * 100000.

   ! Frontogenesis divergence term from RIP/src/fields.f
   where (grad_theta.gt.0.)
      pl3=-.5/grad_theta*(dtheta_Dx*dtheta_Dx + dtheta_Dy*dtheta_Dy)*(dudx+dvdy)
   elsewhere
      pl3=0.
   endwhere
   derived_fields%fronto_div = pl3 * 3.6e8 ! convert from K/m/s to K/100km/hour

   ! Frontogeneis deformation term from RIP/src/fields.f
   where (grad_theta.gt.0.)
      pl4=-.5/grad_theta*(dtheta_Dx*dtheta_Dx - dtheta_Dy*dtheta_Dy)*(dudx-dvdy)
      pl4=pl4-1./grad_theta*dtheta_Dx*dtheta_Dy*(dudy + dvdx)
   elsewhere
      pl4=0.
   endwhere
   derived_fields%fronto_def = pl4*3.6e8 ! convert from K/m/s to K/100km/hour

   ! total frontogenesis
   derived_fields%fronto_tot = derived_fields%fronto_div + derived_fields%fronto_def

   ! clean-up
   deallocate(dtheta_Dx,dtheta_Dy,grad_theta)
   deallocate(dudx, dudy, dvdx, dvdy)
   deallocate(pl3,pl4)

end subroutine rip_frontogensis

end module rip_stuff_mod
