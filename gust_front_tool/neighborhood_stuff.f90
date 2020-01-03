module neighborhood_stuff_mod

use types_mod, only : wrf_static_data
use namelist_and_constants, only : missing_val, pi

implicit none

private

! public subroutines
public :: check_neighborhood_for_event
public :: calc_gaussian_weights, apply_gaussian_filter, deallocate_gaussian_weights
public :: apply_nms, apply_binary_dilation, apply_binary_closing
public :: apply_size_threshold

! define variables visible to all subroutines in this module but are not public 
real, allocatable, dimension(:) :: weights

!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!

subroutine check_neighborhood_for_event(grid,field,event_in_neighborhood,roi_in,threshold)

   type(wrf_static_data), intent(in) :: grid
   real, dimension(grid%nx,grid%ny), intent(in) :: field
   logical, dimension(grid%nx,grid%ny), intent(inout) :: event_in_neighborhood
   real, intent(in) :: roi_in, threshold ! roi_in in km

   real :: roi, distsq, dist, dist2
   integer :: x,y,ng
   integer :: xw,xe,ys,yn,x1,y1,xdiff,ydiff

   roi = roi_in * 1000. ! convert km to m
   ng = int(roi/grid%dx)     ! grid%dx in m
   write(*,*)'dx, ng neigh =',grid%dx, ng

   event_in_neighborhood = .false. ! initialize to false

   do y = 1,grid%ny
      do x = 1,grid%nx

         xw = x-ng
         xe = x+ng
         ys = y-ng
         yn = y+ng

         !  check to see if point (x1,y1) is within bounds of neighborhood and within the model domain. If so, keep going.
         if(xw.ge.1 .and. xe.le.grid%nx .and. ys.ge.1 .and. yn.le.grid%ny) then 
            xloop: do x1 = xw,xe
            yloop:   do y1 = ys,yn
                  xdiff = x-x1
                  ydiff = y-y1
                  distsq = real(xdiff)**2+real(ydiff)**2
                  dist = distsq**0.5
                  dist2 = dist*grid%dx
                  if(dist2.le.roi) then ! If true, then (x1,y1) is in the neighborhood
                     if ( field(x1,y1).gt.threshold) then ! just look for one point in the neighborhood where condition is met
                        event_in_neighborhood(x,y) = .true.
                        exit xloop  ! if we found a point meeting condition, no need to continue searching the neighborhood
                     endif
                  endif ! dist2.le.roi
               enddo yloop !yloop
            enddo xloop ! xloop
         endif ! xw, xe, ys, yn
      enddo !nx
   enddo ! ny

end subroutine check_neighborhood_for_event

subroutine calc_gaussian_weights(grid,roi_in)

   type(wrf_static_data), intent(in) :: grid
   real, intent(in) :: roi_in ! roi_in in km

   real :: roi, distsq, sigma, sigmasq, amp, sqng
   integer :: x,y,ng,ngn,nw
   integer :: xw,xe,ys,yn,x1,y1,xdiff,ydiff

   roi = roi_in * 1000. ! convert km to m

   sigma=roi/grid%dx ! roi is in m and grid%dx is meters
   sigmasq=sigma*sigma
   amp = 1./(2.*pi*sigmasq)

   ng = int(4.*sigma) ! truncate at 4*sigma...only search this size box.
   ngn = (-1)*ng
   sqng = real(ng*ng)

   allocate( weights( (2*ng+1)**2 ) ) ! variable with gaussian weights
   weights = 0.

   nw=0
   do x = ngn,ng
      do y = ngn,ng
         nw=nw+1
         distsq = real(x*x)+real(y*y)
         if(distsq .le. sqng )then ! this makes a circular weighting. points outisde sqng default to 0.
            weights(nw)=amp*EXP(-0.5*distsq/sigmasq)
         endif
      enddo
   enddo

   write(*,*)' size(weights)= ',size(weights)
   write(*,*)' maxval(weights)= ',maxval(weights)
   write(*,*)' weights= ',weights

end subroutine calc_gaussian_weights

subroutine apply_gaussian_filter(grid,roi_in,input_field,output_field)

   type(wrf_static_data), intent(in) :: grid
   real, intent(in) :: roi_in ! roi_in in km
   real, dimension(grid%nx,grid%ny), intent(in) :: input_field
   real, dimension(grid%nx,grid%ny), intent(inout) :: output_field

   real :: roi, sigma
   integer :: x,y,ng,nw
   integer :: xw,xe,ys,yn,x1,y1

   roi = roi_in * 1000. ! convert km to m

   sigma=roi/grid%dx ! roi is in m and grid%dx is meters
   ng = int(4.*sigma)
   write(*,*)'dx, ng gaus =',grid%dx, ng

   ! initialize output to zero
   output_field = 0.0

   ! Now apply the weights to the input
   do y = 1, grid%ny
      do x = 1,grid%nx
         xw = x-ng
         xe = x+ng
         ys = y-ng
         yn = y+ng
         nw = 0
         ! no need to do a distance check because weights will just be 0 if point is outisde the desired radius
         do y1 = ys,yn
            do x1 = xw,xe
               nw = nw+1
               if(x1.ge.1 .and. x1.le.grid%nx .and. y1.ge.1 .and. y1.le.grid%ny) then
                  output_field(x1,y1)=output_field(x1,y1) + weights(nw)*input_field(x,y)
               endif
            enddo
         enddo
      enddo
   enddo
end subroutine apply_gaussian_filter

subroutine deallocate_gaussian_weights
   if ( allocated(weights)) deallocate(weights)
end subroutine deallocate_gaussian_weights

subroutine apply_nms(grid,window_size,input_field,mask_field,mask_thresh,output_field)

   type(wrf_static_data), intent(in) :: grid
   integer, intent(in) :: window_size
   real, dimension(grid%nx,grid%ny), intent(in) :: input_field
   real, dimension(grid%nx,grid%ny), intent(in) :: mask_field
   real, intent(in)    ::  mask_thresh
   real, dimension(grid%nx,grid%ny), intent(inout) :: output_field

   integer :: x,y,j,k
   integer :: x1,y1
   real :: processed = 0.0
   real, allocatable, dimension(:,:,:) :: directional_output

   ! initialize output to zero
   output_field = 0.0

   ! window size is a number of grid boxes symmetrically about a central grid box.
   ! essentially square windows with total number of boxes to search is 2*window_size + 1
   write(*,*)'Total window size for NMS = ',(2*window_size + 1)

   ! allocate output for each direction. 4 angles, hence, first dimension size of 4
   ! then initialize to missing
   allocate( directional_output(4,grid%nx,grid%ny))
   directional_output = missing_val

   ! Now apply the weights to the input
   do y = 1, grid%ny
      do x = 1,grid%nx
         if ( mask_field(x,y) .gt. mask_thresh ) then ! don't need to apply where mask says not to

            ! loop over the 4 possible directions
            ! algorithm from the paper Clark et al. (2015) cited
            ! Sun, C., and P. Vallotton, 2009: Fast linear feature detection using multiple directional non-maximum suppression. J. Microsc., 
            ! 234, 147–157, doi:https://doi.org/10.1111/j.1365-2818.2009.03156.x
            do k = 1,4 
               if ( directional_output(k,x,y) .ne. processed ) then
                  do j = -1*window_size, window_size
                     if ( k.eq.1) then ! 0 degrees (constant y)
                        x1 = x + j
                        y1 = y 
                     endif
                     if ( k.eq.2) then ! 90 degrees (constant x)
                        x1 = x 
                        y1 = y + j
                     endif
                     if ( k.eq.3) then ! 45 degrees
                        x1 = x + j
                        y1 = y + j
                     endif
                     if ( k.eq.4) then ! 135 degrees
                        x1 = x + j
                        y1 = y - j
                     endif
                     if ( x1.eq.x .and. y1.eq.y) cycle ! don't check current pixel
                     if(x1.ge.1 .and. x1.le.grid%nx .and. y1.ge.1 .and. y1.le.grid%ny) then ! make sure point is in domain
                        if ( input_field(x,y) .le. input_field(x1,y1)) then
                           goto 88 ! don't fill directional_output(k,x,y)
                        else
                           directional_output(k,x1,y1) = processed
                        endif
                     endif
                  enddo ! loop over j ( window size)
                  directional_output(k,x,y) = input_field(x,y)
               endif
 88            continue
            enddo ! loop over k
            ! end testing

            ! now just take the union of all the directional output
            ! essentially, we're just keeping all the points where in any of the directions
            !    directional_output(?,x,y) = input_field(x,y)
            if ( any(directional_output(:,x,y) .eq. input_field(x,y)) ) then
               output_field(x,y) = 1.0
            else
               output_field(x,y) = 0.0
            endif

         endif ! check for minimum threshold
      enddo ! loop over x
   enddo ! loop over y

   ! clean-up
   deallocate(directional_output)
end subroutine apply_nms

subroutine apply_binary_dilation(grid,num_dilations,field)
   type(wrf_static_data), intent(in) :: grid
   integer, intent(in) :: num_dilations 
   real, dimension(grid%nx,grid%ny), intent(inout) :: field

   integer :: k,x,y,ng
   integer :: xw,xe,ys,yn,x1,y1
   real, dimension(grid%nx,grid%ny) :: temp_data
   
   ng = 1 ! 3 x 3 kernel

   ! dilation expands lines to larger areas.

   ! initialize temporary data to input data
   temp_data = field

   write(*,*)'Doing ',num_dilations, ' binary dilations'

   do k = 1,num_dilations
      do y = 1, grid%ny
         do x = 1,grid%nx
            if ( field(x,y) .lt. 0.5 ) then ! only process pixels with zeroes. we will reassign the point to 1 if there 1s ANYWHERE in the neighborhood. this is the point of dilation.
               xw = x-ng
               xe = x+ng
               ys = y-ng
               yn = y+ng
               yloop2 : do y1 = ys,yn
                  do x1 = xw,xe
                     if(x1.ge.1 .and. x1.le.grid%nx .and. y1.ge.1 .and. y1.le.grid%ny) then
                        if ( field(x1,y1) .gt. 0.5 ) then ! if any value in the neighborhood is 1, reassign 0 with 1 and go to next point
                           temp_data(x,y) = 1.0
                           exit yloop2 ! if we found a point meeting condition, no need to continue searching the neighborhood
                        endif
                     endif
                  enddo ! x1
               enddo yloop2 ! y1
            endif
         enddo ! x
      enddo ! y
      field = temp_data ! overwrite input
   enddo ! loop over k
end subroutine apply_binary_dilation

subroutine apply_binary_erosion(grid,num_erosions,field)
   type(wrf_static_data), intent(in) :: grid
   integer, intent(in) :: num_erosions
   real, dimension(grid%nx,grid%ny), intent(inout) :: field

   integer :: k,x,y,ng
   integer :: xw,xe,ys,yn,x1,y1
   real, dimension(grid%nx,grid%ny) :: temp_data
   
   ng = 1 ! 3 x 3 kernel

   ! erosion shrinks lines to smaller areas.

   ! initialize temporary data to input data
   temp_data = field

   write(*,*)'Doing ',num_erosions, ' binary erosions'

   do k = 1,num_erosions
      do y = 1, grid%ny
         do x = 1,grid%nx
            if ( field(x,y) .gt. 0.5 ) then ! only process pixels with ones. we will reassign the point to 0 if there are 0s ANYWHERE in the neighborhood. this is the point of erosion.
               xw = x-ng
               xe = x+ng
               ys = y-ng
               yn = y+ng
               yloop3 : do y1 = ys,yn
                  do x1 = xw,xe
                     if(x1.ge.1 .and. x1.le.grid%nx .and. y1.ge.1 .and. y1.le.grid%ny) then
                        if ( field(x1,y1) .lt. 0.5 ) then ! if any value in the neighborhood is 0, reassign 1 with 0 and go to next point
                           temp_data(x,y) = 0.0
                           exit yloop3 ! if we found a point meeting condition, no need to continue searching the neighborhood
                        endif
                     endif
                  enddo ! x1
               enddo yloop3 ! y1
            endif
         enddo ! x
      enddo ! y
      field = temp_data ! overwrite input
   enddo ! loop over k
end subroutine apply_binary_erosion

subroutine apply_binary_closing(grid,num_iterations,field)
   type(wrf_static_data), intent(in) :: grid
   integer, intent(in) :: num_iterations
   real, dimension(grid%nx,grid%ny), intent(inout) :: field

   integer :: k

   write(*,*)'Doing binary closing with ',num_iterations, ' iterations'

   do k = 1,num_iterations
      call apply_binary_dilation(grid,1,field)
      call apply_binary_erosion(grid,1,field)
   enddo

end subroutine apply_binary_closing

subroutine apply_size_threshold(grid,min_size,threshold,field)
   type(wrf_static_data), intent(in) :: grid
   integer, intent(in) :: min_size
   real,    intent(in) :: threshold
   real, dimension(grid%nx,grid%ny), intent(inout) :: field

   integer :: x,y,xy,k
   integer :: xw,xe,ys,yn,x1,y1
   integer :: ifound,ifound2
   integer :: xfound(100000),yfound(100000)
   real, dimension(grid%nx,grid%ny) :: temp_data

   write(*,*)' Removing gust fronts smaller than ',min_size,' points.'

   xy = grid%nx * grid%ny

   ! initialize temporary data to input data because this field gets changed.
   temp_data = field

   do x = 1,grid%nx
      do y = 1,grid%ny
         if ( temp_data(x,y).ge. threshold ) then
            ifound = 0
            xfound(:) = 0
            yfound(:) = 0
            ! now...this is the tricky part...
            ! need to zero out all contiguous grid pts so
            ! the storm is only counted once.
            do x1=x-1,x+1
               do y1=y-1,y+1
                  if(x1.lt.1.or.y1.lt.1) cycle
                  if(x1.gt.grid%nx.or.y1.gt.grid%ny) cycle
                  if(temp_data(x1,y1).ge. threshold) then
                     ifound = ifound + 1
                     xfound(ifound) = x1
                     yfound(ifound) = y1
                     temp_data(x1,y1) = -999. ! make it small...
                  endif
               enddo
            enddo
 600        ifound2 = ifound
            if(ifound.gt.xy) stop 'Array is too small...99999'
            do k=1,ifound
               do x1=xfound(k)-1,xfound(k)+1
                  do y1=yfound(k)-1,yfound(k)+1
                     if(x1.lt.1.or.y1.lt.1) cycle
                     if(x1.gt.grid%nx.or.y1.gt.grid%ny) cycle
                     if(temp_data(x1,y1).ge. threshold ) then
                        ifound = ifound + 1
                        xfound(ifound) = x1
                        yfound(ifound) = y1
                        temp_data(x1,y1) = -999. ! make it small...
                     endif
                  enddo
               enddo
            enddo
            if(ifound.gt.ifound2) goto 600 ! do it again until the "storm"is wiped out.

            ! ifound is the number of contiguous grid pts for this entity
            ! if the number of points is too small, then get rid of the
            !  values in the output. arrays xfound,yfound have the x/y points associated
            !  with the current contiguous feature
            if ( ifound .lt. min_size ) then
               do k = 1,ifound
                  field(xfound(k),yfound(k)) = 0.0
               enddo
            endif

         endif ! fcst_data.ge.0.5
      end do ! y
   end do ! x

end subroutine apply_size_threshold

!!!!
end module neighborhood_stuff_mod
