module namelist_and_constants

use netcdf

implicit none

private

! public subroutines
public ::  read_namelist

! public namelist variables
public :: ens_size
public :: fname_grid_info
public :: search_radius_precip, search_radius_refl
public :: precip_thresh, refl_thresh
public :: wind_diff_min_mag, wind_diff_min_ang
public :: gaussian_length_scale
public :: nms_window_size, num_binary_dilations
public :: min_gf_size

! other public variables
public :: missing_val, pi
public :: fnames_model, fnames_model_prev_hour
public :: fnames_output

! declare namelist variables
integer, parameter :: nmax = 300
integer :: ens_size = 1
character (len = 500) :: fname_grid_info =  ''
character (len = 500) :: input_file_list = '' ! text file containing files to process. not a public variable outside this module.
character (len = 500) :: input_file_list_prev_hour = '' ! text file containing files to process for previous forecast hour. not a public variable outside this module.
character (len = 500) :: output_file_list = ''    ! text file containing output file names.  not a public variable outside this module.
real    :: search_radius_precip = 0  ! radius (km) to search about each grid point for precip
real    :: search_radius_refl = 0  ! radius (km) to search about each grid point for column-max reflectivity
real    :: precip_thresh = 0       ! minimum precipitation to search for in the neighborhood (search_radius_precip) to define event
real    :: refl_thresh   = 0       ! minimum col-max reflectivity to search for in the neighborhood (search_radius_refl) to define event
real    :: wind_diff_min_mag = 0       ! minimum threshold (m/s) of 1-h vector wind difference required for a gust front to be present
real    :: wind_diff_min_ang = 0       ! minimum angle (degrees) of 1-h wind direction shift   required for a gust front to be present
real    :: gaussian_length_scale = 0 ! gaussian length scale for smoothing (km).
integer :: nms_window_size = 0     ! window size for NMS. an integer that is the halfwith of rectange. total number of points is 2*nms_window_size. Set to < 0 to deactivate NMS
integer :: num_binary_dilations = 0 ! number of binary dilations applied to NMS. only active if nms_window_size > 0
integer :: min_gf_size = 0 ! minimum size of the gust front (number of contiguous pixels)

! declare other public variables. these variables also visible to this module
character (len = 500) :: fnames_model(nmax) = ''
character (len = 500) :: fnames_model_prev_hour(nmax) = ''
character (len = 500) :: fnames_output(nmax) = ''
real, parameter :: missing_val = -999999.0
real, parameter :: pi=3.1415926535

namelist /share/ ens_size,  &
                 input_file_list, input_file_list_prev_hour, output_file_list, &
                 fname_grid_info, &
                 search_radius_precip, search_radius_refl, wind_diff_min_mag, wind_diff_min_ang, &
                 precip_thresh, refl_thresh, gaussian_length_scale, nms_window_size, &
                 num_binary_dilations, min_gf_size

!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!

subroutine read_namelist

   integer :: iunit = 912
   integer :: i,j

   open(file='input.nml',unit=iunit, status = 'old')
   read(iunit,nml=share)
   close(iunit)

   ! read the text files containing files to process
   open(unit=57,file=trim(adjustl(input_file_list)),status='old')
   do j = 1,ens_size
      read(57,fmt='(a)')fnames_model(j)
   enddo
   close(57) ! close file

   open(unit=58,file=trim(adjustl(input_file_list_prev_hour)),status='old')
   do j = 1,ens_size
      read(58,fmt='(a)')fnames_model_prev_hour(j)
   enddo
   close(58) ! close file

   ! read the text file containing output file names
   open(unit=59,file=trim(adjustl(output_file_list)),status='old')
   do j = 1,ens_size
      read(59,fmt='(a)')fnames_output(j)
   enddo
   close(59) ! close file

   ! Sanity check regarding number of model files to read
   j = 0
   do i = 1,nmax
      if ( fnames_model(i).eq.'' ) cycle
      j = j + 1
   enddo
   if ( j.ne.ens_size ) then
      write(*,*)'Mismatch regaring ens_size and fnames_model. Change something and try again.'
      write(*,*) j,ens_size
      stop(3)
   endif

   ! Sanity check regarding number of files to write
   j = 0
   do i = 1,nmax
      if ( fnames_output(i).eq.'' ) cycle
      j = j + 1
   enddo
   if ( j.ne.ens_size ) then
      write(*,*)'Mismatch regaring ens_size and fnames_output. Change something and try again.'
      write(*,*) j,ens_size
      stop(3)
   endif

end subroutine read_namelist

end module namelist_and_constants

