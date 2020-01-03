program gf_tool

use namelist_and_constants, only : read_namelist, &
                                   fnames_model, fnames_model_prev_hour, fnames_output, &
                                   fname_grid_info, ens_size, search_radius_precip, search_radius_refl, wind_diff_min_mag, wind_diff_min_ang, &
                                   precip_thresh, refl_thresh, gaussian_length_scale, nms_window_size, num_binary_dilations, &
                                   min_gf_size, pi
use rip_stuff_mod, only : rip_frontogensis, rip_gradient
use types_mod, only : wrf_met_data, derived_data, wrf_static_data, &
                      allocate_met_data, allocate_derived_data, &
                      deallocate_met_data, deallocate_static_data, deallocate_derived_data
use netcdf_io_mod, only : fill_grids_netcdf, get_static_data, &
                          init_netcdf_output, output_netcdf_variable, finalize_netcdf_output
use grib_io_mod, only : fill_grids_grib
use neighborhood_stuff_mod, only : check_neighborhood_for_event, calc_gaussian_weights, apply_gaussian_filter, deallocate_gaussian_weights, &
                                    apply_nms, apply_binary_dilation, apply_binary_closing, apply_size_threshold
use netcdf

implicit none

type(wrf_met_data), dimension(:)  , allocatable :: input_wrf_data, input_wrf_data_prev_hr
type(derived_data), dimension(:)  , allocatable :: derived_fields
type(wrf_static_data) :: grid_static

logical, dimension(:,:)  , allocatable :: event_in_neighborhood
logical, dimension(:)    , allocatable :: is_missing
real, dimension(:,:)     , allocatable :: dummy_output_2d
character(len=3) :: cie  
character(len=500) :: outname

integer :: i, j, ncstatus, ncfileid
integer :: nxm, nym ! number of x and y mass points

!!!!!!!!!!!!!!!!!!!!!!!

! read namelist
call read_namelist

! get static grid information from a netcdf file
call get_static_data(grid_static)
nxm = grid_static%nx  ! local variables with mass point dimensions
nym = grid_static%ny

! calculate gaussian weights for smoothing
if ( gaussian_length_scale.gt.0 ) call calc_gaussian_weights(grid_static,gaussian_length_scale)

! ---------------------------------------------------------------------
! Load up wrf_met_data structure with all the input fields we need
! ---------------------------------------------------------------------
! keep track of whether valid pairs of files are there for each member
allocate(is_missing(ens_size))
allocate(input_wrf_data(ens_size), input_wrf_data_prev_hr(ens_size) )  
allocate(dummy_output_2d(nxm,nym))

do i = 1,ens_size

   ! need files to be present both for this time and the previous forecast hour
   call file_check(fnames_model(i),fnames_model_prev_hour(i),is_missing(i)) ! output is is_missing(i)
   if ( is_missing(i)) cycle

   ! do this time for the ith member
   call allocate_met_data( input_wrf_data(i),nxm,nym,grid_static%dx) ! allocate all data arrays in structure
   ncstatus = nf90_open(path=trim(adjustl(fnames_model(i))),mode=nf90_nowrite,ncid=ncfileid) ! see if it's a netcdf file
   if ( ncstatus.eq.0 ) then 
      call fill_grids_netcdf( ncfileid, trim(adjustl(fnames_model(i))), input_wrf_data(i) )
      ncstatus = nf90_close(ncfileid) ! close file
   else 
      call fill_grids_grib( trim(adjustl(fnames_model(i))), input_wrf_data(i) )
   endif

   ! do the previous forecast hour for the ith member (same as above)
   call allocate_met_data( input_wrf_data_prev_hr(i),nxm,nym,grid_static%dx) 
   ncstatus = nf90_open(path=trim(adjustl(fnames_model_prev_hour(i))),mode=nf90_nowrite,ncid=ncfileid)
   if ( ncstatus.eq.0 ) then 
      call fill_grids_netcdf( ncfileid, trim(adjustl(fnames_model_prev_hour(i))), input_wrf_data_prev_hr(i) )
      ncstatus =  nf90_close(ncfileid) 
   else 
      call fill_grids_grib( trim(adjustl(fnames_model_prev_hour(i))), input_wrf_data_prev_hr(i) ) 
   endif

enddo

! Make sure that at least some good data are there.
if ( all(is_missing)) call fatal_error_message('No valid pairs of input files for any ensemble member')

! ---------------------------------------------------------------------
! Compute derived quantities, including differences across times
! ---------------------------------------------------------------------
! allocate data structure for derived fields. only needed at one time for each member.
allocate(derived_fields(ens_size))
allocate(event_in_neighborhood(nxm,nym))
do i = 1,ens_size 
   if ( .not. is_missing(i)) then
      call allocate_derived_data(derived_fields(i),nxm,nym) ! derived fields on mass points
      call rip_frontogensis(grid_static,input_wrf_data(i),derived_fields(i)) ! fills up derived_fields(i) data structure with lots of fields

      ! get magnitude of the 1-h vector wind difference
      derived_fields(i)%vec_wind_diff_mag = sqrt(  (input_wrf_data(i)%u10 - input_wrf_data_prev_hr(i)%u10)**2 + & 
                                               (input_wrf_data(i)%v10 - input_wrf_data_prev_hr(i)%v10)**2   )

      ! use the dot product to get the angle between the two wind vectors
      derived_fields(i)%vec_wind_diff_angle = &
                    acos( &
                       (input_wrf_data(i)%u10 * input_wrf_data_prev_hr(i)%u10 + input_wrf_data(i)%v10 * input_wrf_data_prev_hr(i)%v10) / &
                       ( sqrt( input_wrf_data(i)%u10 * input_wrf_data(i)%u10 + input_wrf_data(i)%v10 * input_wrf_data(i)%v10 ) * &
                         sqrt( input_wrf_data_prev_hr(i)%u10 * input_wrf_data_prev_hr(i)%u10 + input_wrf_data_prev_hr(i)%v10 * input_wrf_data_prev_hr(i)%v10) ) &
                    )
      derived_fields(i)%vec_wind_diff_angle = derived_fields(i)%vec_wind_diff_angle * (180.0/pi) ! convert to degrees

      ! smooth the frontogenesis field. needs to be applied here before we start setting points = 0
      ! where thresholds aren't met.
      if ( gaussian_length_scale.gt.0 ) then
         call apply_gaussian_filter(grid_static,gaussian_length_scale,derived_fields(i)%fronto_tot,dummy_output_2d) ! smoothed field is in dummy_output_2d
         derived_fields(i)%fronto_tot = dummy_output_2d
      endif

      ! see if precipition exceeding precip_thresh occurs within search_radius_precip (km) of each point
      if ( search_radius_precip .gt. 0 ) then
         ! get rid of grid-point storms (contiguous precip areas >= 7 grid points)
         call apply_size_threshold(grid_static,7,precip_thresh,input_wrf_data(i)%tot_precip) ! input_wrf_data(i)%tot_precip gets updated

         call check_neighborhood_for_event(grid_static,input_wrf_data(i)%tot_precip,event_in_neighborhood,search_radius_precip,precip_thresh) ! output is event_in_neighborhood
         where (.not. event_in_neighborhood ) derived_fields(i)%fronto_tot =  0.0
      endif

      ! see if column-max reflectivity exceeding refl_thresh occurs within search_radius_precip (km) of each point
      if ( search_radius_refl .gt. 0 ) then
         call check_neighborhood_for_event(grid_static,input_wrf_data(i)%col_max_refl,event_in_neighborhood,search_radius_refl,refl_thresh) ! output is event_in_neighborhood
         where (.not. event_in_neighborhood ) derived_fields(i)%fronto_tot = 0.0
      endif

      ! get rid of negative values of frontogensis
     !where ( derived_fields(i)%fronto_tot .le. 0 ) derived_fields(i)%fronto_tot = 0.0

      ! apply minmum 1-h vector wind magnitude threshold
      write(*,*)'removing gust fronts where 1-h vector wind diff magnitude < ',wind_diff_min_mag
      where ( derived_fields(i)%vec_wind_diff_mag .lt. wind_diff_min_mag ) derived_fields(i)%fronto_tot = 0.0

      ! apply minmum 1-h wind direction shift threshold
      write(*,fmt='(a,f6.2,a10)')'removing gust fronts where 1-h wind direction shifted < ',wind_diff_min_ang,' degrees'
      where ( derived_fields(i)%vec_wind_diff_angle .lt. wind_diff_min_ang ) derived_fields(i)%fronto_tot = 0.0

      ! get gradient of frontogenesis.  needed for NMS purposes.
      call rip_gradient(grid_static,derived_fields(i)%fronto_tot,derived_fields(i)%grad_fronto_tot,derived_fields(i)%grad_fronto_tot_angle)
      derived_fields(i)%grad_fronto_tot = derived_fields(i)%grad_fronto_tot * 1000.0 ! get into K/(h*100km^2)

      ! apply NMS ... need to smooth frontogenesis field before applying.
      if ( nms_window_size .gt. 0 ) then 
         ! applying NMS on the gradient of frontogenesis can yield double line strucures, one each at the leading and trailing edge of the boundary
         ! NMS is acting properly in this case, but this result is undesirable.  We just want to get the leading edge and want to ignore trailing edge.
         ! To fix the double line problem, apply NMS to the gradient of the theta field, but only apply NMS (i.e., mask NMS) in areas where
         ! the gradient of the frontogenesis > 0. 
         ! Output of NMS is a binary field of 1s and 0s indicating where a gust front is present

         call apply_nms(grid_static,nms_window_size,derived_fields(i)%grad_theta,derived_fields(i)%grad_fronto_tot,0.0,derived_fields(i)%gust_front)

         if ( num_binary_dilations .gt. 0 ) call apply_binary_dilation(grid_static,num_binary_dilations,derived_fields(i)%gust_front)

        !call apply_binary_closing(grid_static,1,derived_fields(i)%gust_front)

         if ( min_gf_size .gt. 0 ) call apply_size_threshold(grid_static,min_gf_size,0.5,derived_fields(i)%gust_front) ! remove gust fronts < min_gf_size
      endif

   endif 
enddo

! ---------------------------------------------------------------------
! Time to write ... output a netcdf file, one for each member
! ---------------------------------------------------------------------
do i = 1,ens_size
   if ( .not. is_missing(i)) then
    
      outname = fnames_output(i)

      ! open netcdf file for writing; define filename, dimensions, and global attributes
      call init_netcdf_output(grid_static,outname)

      ! write variables to netcdf file
      call output_netcdf_variable(nxm,nym,input_wrf_data(i)%theta2,'TH2','K','2-m potential temperature')
      call output_netcdf_variable(nxm,nym,input_wrf_data(i)%u10,'U10','m/s','10-m zonal wind')
      call output_netcdf_variable(nxm,nym,input_wrf_data(i)%v10,'V10','m/s','10-m meridional wind')
      call output_netcdf_variable(nxm,nym,input_wrf_data(i)%tot_precip,'tot_precip','mm','total accumulated precip')
      call output_netcdf_variable(nxm,nym,input_wrf_data(i)%col_max_refl,'REFL_MAX','dBZ','column-max reflectivity')
      !call output_netcdf_variable(nxm,nym,derived_fields(i)%dtheta_Dx,'dthetaDx','K/100km','dthetaDx')
      !call output_netcdf_variable(nxm,nym,derived_fields(i)%dtheta_Dy,'dthetaDy','K/100km','dthetaDy')
      call output_netcdf_variable(nxm,nym,derived_fields(i)%grad_theta,'gradtheta','K/100km','gradtheta')
      call output_netcdf_variable(nxm,nym,derived_fields(i)%fronto_div,'FRONTdiv','K/100km/h','Frontogenesis divergence term')
      call output_netcdf_variable(nxm,nym,derived_fields(i)%fronto_def,'FRONTdef','K/100km/h','Frontogensis deformation term')
      call output_netcdf_variable(nxm,nym,derived_fields(i)%fronto_tot,'FRONTtot','K/100km/h','Total frontogenesis')
      call output_netcdf_variable(nxm,nym,derived_fields(i)%grad_fronto_tot,'GradFRONTtot','K/100km^2/h','Magnitude of gradient of total frontogenesis')
      call output_netcdf_variable(nxm,nym,derived_fields(i)%grad_fronto_tot_angle,'GradFRONTtotAng','degees','Angle of the gradient of total frontogenesis')
      call output_netcdf_variable(nxm,nym,derived_fields(i)%vec_wind_diff_mag,'VEC_WIND_DIFF_MAG','m/s','1-hr vector wind diff mag')
      call output_netcdf_variable(nxm,nym,derived_fields(i)%vec_wind_diff_angle,'VEC_WIND_DIFF_ANGLE','degrees','1-hr vector wind diff angle')
      call output_netcdf_variable(nxm,nym,derived_fields(i)%gust_front,'binary_gust_front','','Gust front')

      ! close netcdf file
      call finalize_netcdf_output(outname)
   endif
enddo

! -----------------------------------------------------------------------------------
! clean-up ... in the deallocation routines, it makes sure the variable is allocated
!              so no need to make sure every member actually got allocated by using
!              is_missing variable
! -----------------------------------------------------------------------------------
call deallocate_static_data(grid_static)
do i = 1,ens_size
   call deallocate_met_data(input_wrf_data(i))
   call deallocate_met_data(input_wrf_data_prev_hr(i))
   call deallocate_derived_data(derived_fields(i))
enddo
call deallocate_gaussian_weights
deallocate(input_wrf_data, input_wrf_data_prev_hr)
deallocate(derived_fields)
deallocate(is_missing, event_in_neighborhood)
deallocate(dummy_output_2d)

! --------
! all done
! --------
write(*,*)'All done!'

!!!!!!!

contains

subroutine file_check(fname1,fname2,missing)
   character(len=*),intent(in) :: fname1, fname2
   logical , intent(inout)     :: missing
   logical :: lexist1,lexist2

   inquire(file=trim(adjustl(fname1)), exist=lexist1)
   inquire(file=trim(adjustl(fname2)), exist=lexist2)

   missing = .false. ! assume file is there
   if ( (.not. lexist1) .or. (.not. lexist2) ) then
      missing = .true.
      write(0,'(a)')'missing file. skipping '//trim(adjustl(fname1))
   endif

end subroutine file_check

subroutine fatal_error_message(text)
   character(len=*),intent(in) :: text
   write(*,*) trim(adjustl(text))
   stop
end subroutine fatal_error_message

end program gf_tool
