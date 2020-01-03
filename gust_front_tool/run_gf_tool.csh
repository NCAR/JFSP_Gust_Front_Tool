#!/bin/csh
#  run_gf_tool.csh   12/23/19
#   Usage: run_gf_tool.csh

echo =============  Begin run_gf_tool.csh `date` =============

set echo
setenv NCARGFT /Fullpath/to/ncar_gft
set gf_tool_exec = ${NCARGFT}/gf_tool.x  # Path name to the gust front tool executable

set start_init = 2017090800   # initializations to process
set final_init = $start_init
set inc_init = 24             # hour increment

set fhrs = ( 24 ) # forecast hours to process for each initialization

set model = 0   # Type of model output  0 = ncar ensemble  1 = HRRR   2 = URMA
set storage_top = . # Where to put the output (directory will be created if it doesn't exist)
set workdir = $PWD  # Where the code will run (directory will be created if it doesn't exist)

set file_top_dir = /glade/collections/rda/data/ds300.0 # Where the model output files are located

set clean_files = true # if true files that are copied to the working directory will be removed. 
                       # To clean-up and save space.

# Namelist options for the gust-front tool. Values here will overwrite the input.nml file.
set search_radius_precip = 40. # distance (km) within which to search for precip >= precip_thresh
set precip_thresh = 0.10
set search_radius_refl = 40. # distance (km) within which to search for column-max reflectivity >= refl_thresh
set refl_thresh = 10.0
set wind_diff_min_mag = 3.0   # minimum threshold (m/s) of 1-h vector wind difference required for a gust front to be present
set wind_diff_min_ang = 10.0  # minimum wind direction shift of 1-h wind for a gust front to be present.
set gaussian_length_scale = 6 #distance (km) for Gaussian filter to smooth fields. No smoothing if $gaussian_length_scale < 0
set nms_window_size = 10      # number of grid squares for NMS. Value is halfwidth of a rectange centered about the ith point.
set num_binary_dilations = 0  # number of binary dilations to apply to NMS.  Only active if nms_window_size > 0
set min_gf_size = 10          # minimum number of points to count as gust front. Only active if nms_window_siave > 0

set fname_grid_info = ./rt_ensemble_wrfinput_d02 # a NETCDF file on the model's grid

set ens_members = `seq 1 1 10` # The ensemble members to process. Unix command. See 'man seq' for details. 
                               # This command processes all 10 members of the NCAR ensemble

if ( $model == 1 ) then
  set fname_grid_info = ${NCARGFT}/con/hrrr_geo.nc
  set ens_members = 1
endif

# ------------------------------------------------------------------------------------------------
#### Don't need much changing below here other than pointing to the correct directory strucutre
# ------------------------------------------------------------------------------------------------

set ens_size = $#ens_members

mkdir -p $workdir
cd $workdir

#-----------

set init = $start_init
while ( $init <= $final_init )

   set yyyy  = `echo "$init" | cut -c 1-4`
   set sinit = `echo "$init" | cut -c 1-8`
   set hh    = `echo "$init" | cut -c 9-10`
   set hrym  = `echo "$init" | cut -c 3-8`

   set storage = ${storage_top} #/${init}
   mkdir -p $storage

   if ( ( $model == 1 ) then
     set file_dir = ${file_top_dir}
   else
     set file_dir = ${file_top_dir}/${yyyy}/${sinit}
   endif

   echo "Processing ${init}"

   foreach fhr ( $fhrs )

      if ( $fhr == 0 ) continue # don't run tool for forecast hour 0

      set fhr_prev = `expr $fhr - 1` # also need the previous forecast hour
      if ( $model == 1 ) then
        echo processing hrrr output
        set f2 = `printf %02d $fhr`  # get 3-digit forecast hours
        set f2_prev = `printf %02d $fhr_prev`
      endif
      set f3 = `printf %03d $fhr`  # get 3-digit forecast hours
      set f3_prev = `printf %03d $fhr_prev`

      rm -f ./fnames_ensemble.txt
      rm -f ./fnames_ensemble_prev.txt
      rm -f ./fnames_output.txt

      foreach mem ( $ens_members ) 
	 foreach f ( $f3 $f3_prev )  # need to get files for each member for this forecast hour and the previous one.
	    set this_fname = ${file_dir}/diags_d02_${init}_mem_${mem}_f${f}.nc.gz
            if ( $model == 1 ) then
              set hf = `echo $f | cut -c2-3`
              set this_fname = ${file_dir}/hrrr_${hrym}_${hh}_${hf}.grb2
            endif
	    set local_file = `basename $this_fname`    # the filename without the directory and suffix
	    set local_file_unzipped = `basename $this_fname .gz`
	    if ( ! -e $local_file_unzipped ) then # check to see if file is already there and unzipped. if not, unzip it.
	       if ( -e $this_fname ) then
                  if ( $model == 1 ) then
                    ln -s $this_fname ${local_file}
                  else
		    cp $this_fname ./${local_file}      
		    gunzip -f $local_file  # once this is done you'll have $local_file_unzipped
                  endif
	       else
		  echo "$this_fname is missing."
	       endif
	    endif
	    # Build file containing list of ensemble names. Doesn't matter if they're missing.
	    # The tool will handle missing files, but it's important to have one-to-one correspondence
	    # bewteen files/members in the two files for the different forecast hours
	    if ( $f == $f3 )      echo "${local_file_unzipped}" >>! ./fnames_ensemble.txt  # One filename per line
	    if ( $f == $f3_prev ) echo "${local_file_unzipped}" >>! ./fnames_ensemble_prev.txt  # One filename per line
	 end
	 # Build file of output names
	 echo "${storage}/gf_tool_mem${mem}_${init}_f${f3}.nc" >>! ./fnames_output.txt
      end # loop over members

      # Fill the namelist file for the gust-front tool
      rm -f ./input.nml
      cat > ./input.nml << EOF
&share
ens_size=${ens_size},
input_file_list='./fnames_ensemble.txt',
input_file_list_prev_hour='./fnames_ensemble_prev.txt',
output_file_list='./fnames_output.txt'
fname_grid_info='${fname_grid_info}',
search_radius_precip=${search_radius_precip},     ! distance (km) within which to search for precip >= precip_thresh
precip_thresh = ${precip_thresh}
search_radius_refl=${search_radius_refl}       ! distance (km) within which to search for column-max reflectivity >= refl_thresh
refl_thresh = ${refl_thresh},
wind_diff_min_mag = ${wind_diff_min_mag},
wind_diff_min_ang = ${wind_diff_min_ang},
gaussian_length_scale=${gaussian_length_scale}
nms_window_size=${nms_window_size}
num_binary_dilations=${num_binary_dilations}
min_gf_size=${min_gf_size}
/
EOF

      # Run the gust-front tool!
#     ln -sf $gf_tool_exec .
      $gf_tool_exec # >& log_gf_tool.out 

      # clean-up 
      if ( $clean_files == true || $clean_files == .true. ) then
	 rm -f ./diags_d02_${init}_mem*f${f3}.nc
	 rm -f ./diags_d02_${init}_mem*f${f3_prev}.nc
      endif

      rm -f ./core*

   end # loop over forecast hours

   # Go to next initialization
   set init = `date +%Y%m%d%H -d "$sinit $hh + ${inc_init} hour"`

end # loop over inits

echo =============  Completed run_gf_tool.csh `date` =============
exit
