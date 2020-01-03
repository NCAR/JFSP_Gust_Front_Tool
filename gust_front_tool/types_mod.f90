module types_mod

!use netcdf

implicit none

private

! public stuff
public :: wrf_met_data, derived_data, wrf_static_data
public :: deallocate_met_data, deallocate_static_data, deallocate_derived_data
public :: allocate_met_data, allocate_derived_data

type wrf_met_data  
     integer                    :: nx   ! number of mass points
     integer                    :: ny   ! number of mass points
     real                       :: dx   ! horizontal grid spacing in meters
     real,                      dimension(:,:),    allocatable :: u10 ! 10-m zonal wind on mass points
     real,                      dimension(:,:),    allocatable :: v10 ! 10-m meridional wind on mass points
     real,                      dimension(:,:),    allocatable :: t2 ! 2-mtemperature (K) on mass points
     real,                      dimension(:,:),    allocatable :: theta2 ! 2-m potential temperature (K) on mass points
     real,                      dimension(:,:),    allocatable :: q2 ! 2-m water vapor mixing ratio (kg/kg) on mass points
     real,                      dimension(:,:),    allocatable :: psfc ! surface pressure (Pa) on mass points
     real,                      dimension(:,:),    allocatable :: col_max_refl ! column maximum reflectivity (dBz)
     real,                      dimension(:,:),    allocatable :: tot_precip ! total 1-h accumulated precip on mass points (mm)
end type wrf_met_data

type derived_data
     integer                    :: nx   ! number of mass points
     integer                    :: ny   ! number of mass points
     real                       :: dx   ! horizontal grid spacing in meters
     real,                      dimension(:,:),    allocatable :: fronto_div ! Miller frontogensis divergence  term (K/100km/hour)
     real,                      dimension(:,:),    allocatable :: fronto_def ! Miller frontogensis deformation term (K/100km/hour)
     real,                      dimension(:,:),    allocatable :: fronto_tot ! Sum of fronto_div + fronto_def
     real,                      dimension(:,:),    allocatable :: grad_fronto_tot ! magnitude of the gradient of fronto_tot
     real,                      dimension(:,:),    allocatable :: grad_fronto_tot_angle ! magnitude of the gradient of fronto_tot
     real,                      dimension(:,:),    allocatable :: vec_wind_diff_mag ! Magnitude of the vector wind difference between current hour and previous hour
     real,                      dimension(:,:),    allocatable :: vec_wind_diff_angle ! Angle bewteen the wind vectors for current hour and previous hour
     real,                      dimension(:,:),    allocatable :: dtheta_Dx
     real,                      dimension(:,:),    allocatable :: dtheta_Dy
     real,                      dimension(:,:),    allocatable :: grad_theta ! magnitude of potential temperature gradient
     real,                      dimension(:,:),    allocatable :: gust_front ! binary (yes/no) occurrence of gust front. 1s and 0s, but keep as real.
end type derived_data

type wrf_static_data
   integer                    :: nx   ! number of mass points
   integer                    :: ny   ! number of mass points
   real                       :: dx   ! horizontal grid spacing in meters
   real,                      dimension(:,:),    allocatable :: mapfac_m ! mapfactor at mass points
   real,                      dimension(:,:),    allocatable :: mapfac_u ! mapfactor at u-staggered points
   real,                      dimension(:,:),    allocatable :: mapfac_v ! mapfactor at v-staggered points
   real,                      dimension(:,:),    allocatable :: ter      ! model terrain height (m)
end type wrf_static_data

contains

! ----------------------
! allocation routines 
! ----------------------
subroutine allocate_derived_data(grid,nx,ny)
 type(derived_data), intent(inout) :: grid
 integer, intent(in)               :: nx,ny
 grid%nx = nx
 grid%ny = ny
 if(.not. allocated(grid%fronto_div))   allocate(grid%fronto_div(nx,ny))
 if(.not. allocated(grid%fronto_def))   allocate(grid%fronto_def(nx,ny))
 if(.not. allocated(grid%fronto_tot))   allocate(grid%fronto_tot(nx,ny))
 if(.not. allocated(grid%grad_fronto_tot))   allocate(grid%grad_fronto_tot(nx,ny))
 if(.not. allocated(grid%grad_fronto_tot_angle))   allocate(grid%grad_fronto_tot_angle(nx,ny))
 if(.not. allocated(grid%vec_wind_diff_mag))   allocate(grid%vec_wind_diff_mag(nx,ny))
 if(.not. allocated(grid%vec_wind_diff_angle))   allocate(grid%vec_wind_diff_angle(nx,ny))
 if(.not. allocated(grid%dtheta_Dx))       allocate(grid%dtheta_Dx(nx,ny))
 if(.not. allocated(grid%dtheta_Dy))       allocate(grid%dtheta_Dy(nx,ny))
 if(.not. allocated(grid%grad_theta))      allocate(grid%grad_theta(nx,ny))
 if(.not. allocated(grid%gust_front))      allocate(grid%gust_front(nx,ny))
end subroutine allocate_derived_data

subroutine allocate_met_data(grid,nx,ny,dx)
 type(wrf_met_data), intent(inout) :: grid
 integer, intent(in)               :: nx, ny 
 real,    intent(in)               :: dx
 grid%nx = nx
 grid%ny = ny
 grid%dx = dx
 if(.not. allocated(grid%u10))         allocate(grid%u10(nx,ny))
 if(.not. allocated(grid%v10))         allocate(grid%v10(nx,ny))
 if(.not. allocated(grid%t2))          allocate(grid%t2(nx,ny))
 if(.not. allocated(grid%theta2))      allocate(grid%theta2(nx,ny))
 if(.not. allocated(grid%q2))          allocate(grid%q2(nx,ny))
 if(.not. allocated(grid%psfc))        allocate(grid%psfc(nx,ny))
 if(.not. allocated(grid%col_max_refl))        allocate(grid%col_max_refl(nx,ny))
 if(.not. allocated(grid%tot_precip))  allocate(grid%tot_precip(nx,ny))
end subroutine allocate_met_data

! ----------------------
! deallocation routines
! ----------------------
subroutine deallocate_met_data(grid)
 type(wrf_met_data), intent(inout) :: grid
 if(allocated(grid%u10))          deallocate(grid%u10)
 if(allocated(grid%v10))          deallocate(grid%v10)
 if(allocated(grid%t2))           deallocate(grid%t2)
 if(allocated(grid%theta2))       deallocate(grid%theta2)
 if(allocated(grid%q2))           deallocate(grid%q2)
 if(allocated(grid%psfc))         deallocate(grid%psfc)
 if(allocated(grid%col_max_refl))    deallocate(grid%col_max_refl)
 if(allocated(grid%tot_precip))   deallocate(grid%tot_precip)
end subroutine deallocate_met_data

subroutine deallocate_derived_data(grid)
 type(derived_data), intent(inout) :: grid
 if(allocated(grid%fronto_div))   deallocate(grid%fronto_div)
 if(allocated(grid%fronto_def))   deallocate(grid%fronto_def)
 if(allocated(grid%fronto_tot))   deallocate(grid%fronto_tot)
 if(allocated(grid%grad_fronto_tot))   deallocate(grid%grad_fronto_tot)
 if(allocated(grid%grad_fronto_tot_angle))   deallocate(grid%grad_fronto_tot_angle)
 if(allocated(grid%vec_wind_diff_mag))   deallocate(grid%vec_wind_diff_mag)
 if(allocated(grid%vec_wind_diff_angle))   deallocate(grid%vec_wind_diff_angle)
 if(allocated(grid%dtheta_Dx))   deallocate(grid%dtheta_Dx)
 if(allocated(grid%dtheta_Dy))   deallocate(grid%dtheta_Dy)
 if(allocated(grid%grad_theta))   deallocate(grid%grad_theta)
 if(allocated(grid%gust_front))   deallocate(grid%gust_front)
end subroutine deallocate_derived_data

subroutine deallocate_static_data(grid)
 type(wrf_static_data), intent(inout) :: grid
 if(allocated(grid%mapfac_m))     deallocate(grid%mapfac_m)
 if(allocated(grid%mapfac_u))     deallocate(grid%mapfac_u)
 if(allocated(grid%mapfac_v))     deallocate(grid%mapfac_v)
 if(allocated(grid%ter))     deallocate(grid%ter)
end subroutine deallocate_static_data

end module types_mod
