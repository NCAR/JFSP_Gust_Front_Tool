module netcdf_io_mod

use namelist_and_constants, only : fname_grid_info

use rip_stuff_mod, only :  rip_potential_temperature
use types_mod, only : wrf_met_data, wrf_static_data

use netcdf

implicit none

private

! public subroutines
public :: fill_grids_netcdf, get_static_data
public :: init_netcdf_output, output_netcdf_variable, finalize_netcdf_output

! variables visible by routines in this module but not public
integer :: ncstatus, ncvarid 
integer :: ncid_out
integer, dimension(2) :: id2

!!!!!!!!!!!!!!!!!!!!!
contains

subroutine get_static_data(grid)

   type(wrf_static_data), intent(inout) :: grid

   integer  :: ncfileid_static
   character(len=500)     :: netcdf_var
   integer  :: nx, ny

   !----------------------------------------------------------------------------------
   ! get static information about the grid, like map factors, etc
   !----------------------------------------------------------------------------------
   ! try opening fname_grid_info, which has static grid information
   ncstatus = nf90_open(path=trim(adjustl(fname_grid_info)),mode=nf90_nowrite,ncid=ncfileid_static)
   if ( ncstatus.ne.0 ) then
      write(0,'(a)')' Problem opening '//trim(adjustl(fname_grid_info))
      stop
   endif

   ! horizontal grid spacing
   ncstatus = nf90_get_att(ncfileid_static,nf90_global,'DX',grid%dx) ! get dx from netcdf file global attribute
   if ( ncstatus /= 0 ) then
      write(0,*) 'Error getting global attribute DX from netCDF file'
      stop
   else
      write(0,fmt='(a,f8.3)') 'Horizontal grid spacing is ',grid%dx
      write(*,*) ''
   endif

   ! get mapfactor on mass point
   netcdf_var = 'MAPFAC_M'
   call get_var_dims_netcdf(ncfileid_static,trim(adjustl(fname_grid_info)),netcdf_var,nx,ny)
   allocate( grid%mapfac_m(nx,ny))
   call get_netcdf_var_2d_real(ncfileid_static,fname_grid_info,trim(adjustl(netcdf_var)),nx,ny,grid%mapfac_m)

   ! mapfac_m is on mass points, so assign these to the data structure
   grid%nx = nx
   grid%ny = ny

   ! get mapfactor on u-staggered points
   netcdf_var = 'MAPFAC_U'
   call get_var_dims_netcdf(ncfileid_static,trim(adjustl(fname_grid_info)),netcdf_var,nx,ny)
   allocate( grid%mapfac_u(nx,ny))
   call get_netcdf_var_2d_real(ncfileid_static,fname_grid_info,trim(adjustl(netcdf_var)),nx,ny,grid%mapfac_u)

   ! get mapfactor on v-staggered points
   netcdf_var = 'MAPFAC_V'
   call get_var_dims_netcdf(ncfileid_static,trim(adjustl(fname_grid_info)),netcdf_var,nx,ny)
   allocate( grid%mapfac_v(nx,ny))
   call get_netcdf_var_2d_real(ncfileid_static,fname_grid_info,trim(adjustl(netcdf_var)),nx,ny,grid%mapfac_v)

   ! get model terrain height
   netcdf_var = 'HGT'
   call get_var_dims_netcdf(ncfileid_static,trim(adjustl(fname_grid_info)),netcdf_var,nx,ny)
   allocate( grid%ter(nx,ny))
   call get_netcdf_var_2d_real(ncfileid_static,fname_grid_info,trim(adjustl(netcdf_var)),nx,ny,grid%ter)

   ncstatus =  nf90_close(ncfileid_static) ! close file...done with static file

end subroutine get_static_data

subroutine fill_grids_netcdf(ncfileid,fname,grid)
   integer, intent(in)           :: ncfileid
   character(len=*), intent(in)  :: fname
   type(wrf_met_data), intent(inout) :: grid

   character(len=500)     :: netcdf_var
   integer                :: nx, ny
   real, dimension(:,:), allocatable :: work2d

   !----------------------------------------------------------------------------------
   ! load up the wrf_met_data structure with meteorological data variable by variable
   ! dimensions have already been allocated (all mass points)
   !----------------------------------------------------------------------------------

   ! input grid has nx,ny attached to it
   nx = grid%nx
   ny = grid%ny

   ! 2-m temperature
   netcdf_var = 'T2'
   call get_netcdf_var_2d_real(ncfileid,fname,trim(adjustl(netcdf_var)),nx,ny,grid%t2)

   ! 2-m water vapor mixing ratio (kg/kg)
   netcdf_var = 'Q2'
   call get_netcdf_var_2d_real(ncfileid,fname,trim(adjustl(netcdf_var)),nx,ny,grid%q2)

   ! 10-m u-wind
   netcdf_var = 'U10'
   call get_netcdf_var_2d_real(ncfileid,fname,trim(adjustl(netcdf_var)),nx,ny,grid%u10)

   ! 10-m v-wind
   netcdf_var = 'V10'
   call get_netcdf_var_2d_real(ncfileid,fname,trim(adjustl(netcdf_var)),nx,ny,grid%v10)

   ! surface pressure (Pa)
   netcdf_var = 'PSFC'
   call get_netcdf_var_2d_real(ncfileid,fname,trim(adjustl(netcdf_var)),nx,ny,grid%psfc)

   ! column maximum reflectivity (dBz)
   netcdf_var = 'REFD_MAX'
   call get_netcdf_var_2d_real(ncfileid,fname,trim(adjustl(netcdf_var)),nx,ny,grid%col_max_refl)

   ! get 2-m potential temperature (K) on mass points using RIP formula
   call rip_potential_temperature(grid)

   ! total 1-h precipitation (mm), which is really the sum of PREC_ACC_NC and PREC_ACC_C
   ! quite possible that PREC_ACC_C isn't in the file...if so get_netcdf_var_2d_real should return array of 0, which is ok
   netcdf_var = 'PREC_ACC_NC'
   call get_netcdf_var_2d_real(ncfileid,fname,trim(adjustl(netcdf_var)),nx,ny,grid%tot_precip)

   netcdf_var = 'PREC_ACC_C'
   allocate( work2d(nx,ny))
   call get_netcdf_var_2d_real(ncfileid,fname,trim(adjustl(netcdf_var)),nx,ny,work2d)
   grid%tot_precip = grid%tot_precip + work2d
   deallocate(work2d)

end subroutine fill_grids_netcdf

subroutine get_var_dims_netcdf(ncfileid,fname,netcdf_var,nx_out,ny_out)
   implicit none

   integer, intent(in)           :: ncfileid
   character(len=*), intent(in)  :: fname
   character(len=*), intent(inout)  :: netcdf_var 
   integer, intent(inout) :: nx_out, ny_out

   integer, dimension(4) :: dims
   integer :: ndims, i,xtype,natts,dimids(10)
   character (len=80) :: varname

   ! get dimensions for this variable
   ncstatus = nf90_inq_varid( ncfileid, trim(adjustl(netcdf_var)), ncvarid)
   if ( ncstatus /= 0 ) then
      if ( trim(adjustl(netcdf_var)).eq.'HGT') then  ! try another variable for HGT
         netcdf_var = 'HGT_M'
         ncstatus = nf90_inq_varid( ncfileid, trim(adjustl(netcdf_var)), ncvarid)
      endif
      if ( ncstatus /= 0 ) then
         write(*,'(a)') 'Variable '//trim(adjustl(netcdf_var))//' not found in '//trim(adjustl(fname))//'.  Stop.'
         stop
      endif
   endif
   ncstatus = nf90_Inquire_Variable(ncfileid, ncvarid, varname, xtype, ndims, dimids, natts)
   do i=1,ndims
      ncstatus = nf90_inquire_dimension( ncfileid, dimids(i), len=dims(i) )
   enddo
   nx_out = dims(1)
   ny_out = dims(2)
   write(*,fmt='( 2(a,1x),2i5)') trim(adjustl(netcdf_var)),trim(adjustl(fname)),nx_out,ny_out

end subroutine get_var_dims_netcdf

subroutine get_netcdf_var_2d_real(fileid,fname,variable,dim1,dim2,output)

   integer, intent(in) :: fileid, dim1,dim2
   character(len=*), intent(in) :: fname, variable
   real, intent(inout), dimension(dim1,dim2)  :: output

   integer :: istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   if ( ncstatus /= 0 ) then
      write(*,fmt='(a)') 'Error inq_varid '//trim(adjustl(variable))//' from '//trim(adjustl(fname))
      istatus = istatus + ncstatus
   endif
   ncstatus = nf90_get_var(fileid,ncvarid,output)
   if ( ncstatus /= 0 ) then
      write(*,fmt='(a)') 'Error inq_get_var '//trim(adjustl(variable))//' from '//trim(adjustl(fname))
      istatus = istatus + ncstatus
   endif

   if ( istatus /= 0 ) then
      output = 0.0  ! Set to zero if error
   else
      write(*,'(a,2(f10.3,1x))')'Min/max input vals of '//trim(adjustl(variable)),minval(output),maxval(output)
   endif

end subroutine get_netcdf_var_2d_real

subroutine init_netcdf_output(grid,outname)

   type(wrf_static_data), intent(in) :: grid
   character(len=*), intent(in) :: outname

   integer :: nx_id, ny_id

   ! ---------------------------------------------
   ! open output file ... opened in "define mode" 
   ! ---------------------------------------------
   ncstatus = nf90_create(trim(adjustl(outname)), NF90_CLOBBER, ncid_out)
   if ( ncstatus.ne.0 ) write(0,'(a)')'Error opening '//trim(adjustl(outname))//' for writing.'

   ! ---------------------------------------------------------------------------
   ! declare the dimensions of the output file ... file already in "define mode"
   ! ---------------------------------------------------------------------------
   ncstatus = nf90_def_dim(ncid_out, 'east_west',  grid%nx, nx_id)
   ncstatus = nf90_def_dim(ncid_out, 'north_south',grid%ny ,ny_id)
   id2(1) = nx_id
   id2(2) = ny_id

   ! -------------------
   ! global attributes 
   ! -------------------
   ncstatus = nf90_put_att(ncid_out, NF90_GLOBAL, 'history', 'test file for JFSP gust front tool')

   ! --------------------
   ! leave define mode
   ! --------------------
   ncstatus = nf90_enddef(ncid_out)

end subroutine init_netcdf_output

subroutine output_netcdf_variable(nx,ny,output_data,netcdf_var,units,long_name)

   integer, intent(in) :: nx, ny
   real, dimension(nx,ny), intent(in) :: output_data
   character(len=*), intent(in) :: netcdf_var, units, long_name
   
   integer :: var_id 

   ! put file back into define mode so dimensions/attributes for this variable can be declared
   ncstatus = nf90_redef(ncid_out)

   ! define variable, its dimensions, and attributes
   ncstatus = nf90_def_var(ncid_out, trim(adjustl(netcdf_var)), NF90_FLOAT, id2, var_id)
   ncstatus = nf90_put_att(ncid_out, var_id, 'units', trim(adjustl(units)) )
   ncstatus = nf90_put_att(ncid_out, var_id, 'long_name', trim(adjustl(long_name)) )

   ! leave define mode so data can be output
   ncstatus = nf90_enddef(ncid_out) 

   ! now actually output the data
   ncstatus = nf90_put_var(ncid_out, var_id, output_data)

end subroutine output_netcdf_variable

subroutine finalize_netcdf_output(outname)
   character(len=*), intent(in) :: outname
   ! close file. done writing
   ncstatus = nf90_close(ncid_out)
   if ( ncstatus .ne. 0 ) write(0,*)'Error closing '//trim(adjustl(outname))
end subroutine finalize_netcdf_output

!!!!!!!!!!!!!!!!!!!!!!!

end module netcdf_io_mod
