module grib_io_mod

use rip_stuff_mod, only : rip_potential_temperature
use types_mod, only : wrf_met_data

use grib_mod   ! GRIB2 module/library

implicit none

private

! public subroutines
public :: fill_grids_grib

! variables visible by routines in this module but not public
type(gribfield) :: gfld
integer :: kr,kfo,irets,ierr,pds_template,gds_template
integer :: jids(200), jpds(200),jgds(200),kpds(200),kgds(200)
character (len = 500) :: outname,field_name

integer :: iunit = 0

!!!!!!!!!!!!!!!!!!!!!

contains

subroutine fill_grids_grib(fname,grid)
   character(len=*), intent(in)  :: fname
   type(wrf_met_data), intent(inout) :: grid

   real, dimension(:)  , allocatable :: input_data(:)
   integer                           :: nx, ny

   !----------------------------------------------------------------------------------
   ! load up the wrf_met_data structure with meteorological data variable by variable
   ! dimensions have already been allocated (all mass points)
   !----------------------------------------------------------------------------------

   ! input grid has nx,ny attached to it
   nx = grid%nx
   ny = grid%ny
   allocate(input_data(nx*ny))

   iunit=iunit + 1 ! get unique unit for each file that is read

   ! these stay the same for all fields
   jids=-9999
   jpds=-9999
   jgds=-9999 
   gds_template=30 ! lambert conformal ... could possibly just set to -9999 (or maybe -1)

   ! open the GRIB2 file for reading
   call baopenr(iunit,trim(adjustl(fname)),ierr)
   if ( ierr.ne.0 ) then
      write(*,fmt='(a)') 'Error opening ',trim(adjustl(fname))
      stop
   endif

   ! 2-m temperature (grid%t2)
   jpds(1)=0    ! category
   jpds(2)=0    ! parameter
   jpds(10)=103 ! level type (level 1)
   jpds(12)=2   ! level 1 value
   pds_template = 0  ! means template 4.0
   call getgb2(iunit,0,0,-1,jids,pds_template,jpds,gds_template,jgds,.true.,kfo,gfld,irets)
   call error_check(fname,'T2',irets)
   input_data = gfld%fld
   call convert_1d_to_2d(grid,input_data,grid%t2)
   write(*,'(a,2(f10.3,1x))')'Min/max input vals of T2 ',minval(grid%t2),maxval(grid%t2)

   ! 2-m water vapor mixing ratio (kg/kg)
   ! we're actually getting specific humidity from the GRIB2 file, so convert to mixing ratio
   jpds(1)=1    ! category
   jpds(2)=0    ! parameter
   jpds(10)=103 ! level type (level 1)
   jpds(12)=2   ! level 1 value
   pds_template = 0  ! means template 4.0
   call getgb2(iunit,0,0,-1,jids,pds_template,jpds,gds_template,jgds,.true.,kfo,gfld,irets)
   call error_check(fname,'Q2',irets)
   input_data = gfld%fld
   call convert_1d_to_2d(grid,input_data,grid%q2)
   grid%q2 = grid%q2 / ( 1.0 - grid%q2 ) ! specific humidity to mixing ratio
   write(*,'(a,2(f10.3,1x))')'Min/max input vals of Q2 ',minval(grid%q2),maxval(grid%q2)
   
   ! 10-m u-wind
   jpds(1)=2    ! category
   jpds(2)=2    ! parameter
   jpds(10)=103 ! level type (level 1)
   jpds(12)=10  ! level 1 value
   pds_template = 0  ! means template 4.0
   call getgb2(iunit,0,0,-1,jids,pds_template,jpds,gds_template,jgds,.true.,kfo,gfld,irets)
   call error_check(fname,'U10',irets)
   input_data = gfld%fld
   call convert_1d_to_2d(grid,input_data,grid%u10)
   write(*,'(a,2(f10.3,1x))')'Min/max input vals of U10 ',minval(grid%u10),maxval(grid%u10)

   ! 10-m v-wind
   jpds(1)=2    ! category
   jpds(2)=3    ! parameter
   jpds(10)=103 ! level type (level 1)
   jpds(12)=10  ! level 1 value
   pds_template = 0  ! means template 4.0
   call getgb2(iunit,0,0,-1,jids,pds_template,jpds,gds_template,jgds,.true.,kfo,gfld,irets)
   call error_check(fname,'V10',irets)
   input_data = gfld%fld
   call convert_1d_to_2d(grid,input_data,grid%v10)
   write(*,'(a,2(f10.3,1x))')'Min/max input vals of V10 ',minval(grid%v10),maxval(grid%v10)

   ! surface pressure (Pa)
   jpds(1)=3    ! category
   jpds(2)=0    ! parameter
   jpds(10)=1 ! level type (level 1)
   jpds(12)=-9999   ! level 1 value
   pds_template = 0  ! means template 4.0
   call getgb2(iunit,0,0,-1,jids,pds_template,jpds,gds_template,jgds,.true.,kfo,gfld,irets)
   call error_check(fname,'PSFC',irets)
   input_data = gfld%fld
   call convert_1d_to_2d(grid,input_data,grid%psfc)
   write(*,'(a,2(f10.3,1x))')'Min/max input vals of PSFC ',minval(grid%psfc),maxval(grid%psfc)

   ! column maximum reflectivity (dBZ)
   jpds(1)=16   ! category
   jpds(2)=196  ! parameter
   jpds(10)=10 ! level type (level 1)
   jpds(12)=-9999   ! level 1 value
   pds_template = 0  ! means template 4.0
   call getgb2(iunit,0,0,-1,jids,pds_template,jpds,gds_template,jgds,.true.,kfo,gfld,irets)
   if ( irets.ne.0 ) then ! if there's an error try different GRIB2 parameters for URMA
      jpds(1)=10   ! category
      jpds(2)=0  ! parameter
      jpds(10)=-9999 ! level type (level 1)
      jpds(12)=-9999   ! level 1 value
      pds_template = 0  ! means template 4.0
      call getgb2(iunit,0,0,-1,jids,pds_template,jpds,gds_template,jgds,.true.,kfo,gfld,irets)
      call error_check(fname,'REFD_MAX',irets)
   endif
   input_data = gfld%fld
   call convert_1d_to_2d(grid,input_data,grid%col_max_refl)
   write(*,'(a,2(f10.3,1x))')'Min/max input vals of REFD_MAX ',minval(grid%col_max_refl),maxval(grid%col_max_refl)

   ! get 2-m potential temperature (K) on mass points using RIP formula
   call rip_potential_temperature(grid)

   ! total 1-hr precipitation (mm)
   ! this should go last, since we're assigning jpds(27) a value
   jpds(1)=1    ! category
   jpds(2)=8    ! parameter
   jpds(10)=1 ! level type (level 1)
   jpds(12)=-9999   ! level 1 value
   pds_template = 8  ! means template 4.8 for accumulated quantity
   jpds(27) = 1 ! indicates 1-h accumulation, incase there are multiple precip records in file
   call getgb2(iunit,0,0,-1,jids,pds_template,jpds,gds_template,jgds,.true.,kfo,gfld,irets)
   if ( irets.ne.0 ) then ! if there's an error try different GRIB2 parameters for URMA
      jpds(1)=6   ! category
      jpds(2)=9  ! parameter
      jpds(10)=-9999 ! level type (level 1)
      jpds(12)=-9999   ! level 1 value
      pds_template = 0  ! MRMS is template 4.0
      jpds(27) = -9999 ! indicates 1-h accumulation, incase there are multiple precip records in file
      call getgb2(iunit,0,0,-1,jids,pds_template,jpds,gds_template,jgds,.true.,kfo,gfld,irets)
      call error_check(fname,'tot precip',irets)
   endif
   input_data = gfld%fld
   call convert_1d_to_2d(grid,input_data,grid%tot_precip)
   write(*,'(a,2(f10.3,1x))')'Min/max input vals of tot_precip ',minval(grid%tot_precip),maxval(grid%tot_precip)

   ! done getting variales.  now close GRIB2 file
   call baclose(iunit,ierr)

   ! clean-up
   call gf_free(gfld) ! Deallocate grib2 structure
   deallocate(input_data)

end subroutine fill_grids_grib

subroutine convert_1d_to_2d(grid,input_data_1d,input_data_2d)
   type(wrf_met_data), intent(in) :: grid
   real,dimension(grid%nx*grid%ny),intent(in) :: input_data_1d
   real,dimension(grid%nx,grid%ny),intent(inout) :: input_data_2d

   integer :: i,j,kk

   ! make 1d input data 2d.
   kk = 0
   do j = 1,grid%ny
      do i = 1,grid%nx
         kk = kk + 1
         input_data_2d(i,j) = input_data_1d(kk)
      enddo
   enddo
end subroutine convert_1d_to_2d

subroutine error_check(fname,variable,irets)
   character(len=*), intent(in)  :: fname, variable
   integer, intent(in)           :: irets
   if ( irets.ne.0 ) then
      write(*,fmt='(a)') 'Error reading '//trim(adjustl(variable))//'from '//trim(adjustl(fname))
      stop
   endif
end subroutine error_check

!!!!!!!!!!!!!!!!!!!!!!!

end module grib_io_mod
