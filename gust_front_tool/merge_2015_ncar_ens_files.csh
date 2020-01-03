#!/bin/csh

# For NCAR ensemble data before November-ish 2015, there is no PSFC variable in the netCDF "diag" files
# But, the gust-front tool needs PSFC
# So, get the surface pressure from the GRIB file and stuff it into the netCDF file

set time_exec =    /glade/u/home/schwartz/utils/da_advance_time.exe # Utility to manipulate/advance date/time

set start_init = 2015073000   # initializations to process
set final_init = $start_init
set inc_init = 24

set fhrs = ( `seq 1 1 36` ) # forecast hours to process for each initialization

set workdir = /glade/scratch/schwartz  # Where the code will run (directory will be created if it doesn't exist)
set storage_top = ${workdir} # Where to put the output (directory will be created if it doesn't exist)

set file_top_dir = /glade/collections/rda/data/ds300.0 # Where files are located

set fname_grid_info = /glade/p/mmm/parc/schwartz/rt_ensemble/wrfinput_d02 # a NETCDF file on the model's grid

set ens_members = `seq 1 1 10` # The ensemble members to process 

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

   set storage = ${storage_top} #/${init}
   mkdir -p $storage

   set file_dir = ${file_top_dir}/${yyyy}/${sinit}

   echo "Processing ${init}"

   foreach fhr ( $fhrs )

      if ( $fhr == 0 ) continue # don't run tool for forecast hour 0

      set f3 = `printf %03d $fhr`  # get 3-digit forecast hours

      foreach mem ( $ens_members ) 
	 set fname_diag   = ${file_dir}/diags_d02_${init}_mem_${mem}_f${f3}.nc.gz
	 set fname_grib   = ${file_dir}/ncar_3km_${init}_mem${mem}_f${f3}.grb.gz
	 foreach this_fname ( $fname_diag   $fname_grib )
	    set local_file = `basename $this_fname`
	    set local_file_unzipped = `basename $this_fname .gz`
	    if ( ! -e $local_file_unzipped ) then # check to see if file is already there and unzipped. if not, unzip it.
	       if ( -e $this_fname ) then
		  cp $this_fname ./${local_file}      
		  gunzip -f $local_file  # once this is done you'll have $local_file_unzipped
	       else
		  echo "$this_fname is missing."
	       endif
	    endif
	 end

	 # Convert the GRIB file to netCDF
	 # Only need variable PRES_GDS3_SFC
	 set fname_grib =   ./ncar_3km_${init}_mem${mem}_f${f3}.grb
	 set fname_netcdf = ./ncar_3km_${init}_mem${mem}_f${f3}.nc
	 ncl_convert2nc $fname_grib -v PRES_GDS3_SFC # Output is $fname_netcdf

	 # Now rename dimensions and then add a new dimension for Time
	 ncrename -v PRES_GDS3_SFC,PSFC    $fname_netcdf
	 ncrename -O -d g3_x_0,south_north $fname_netcdf
	 ncrename -O -d g3_y_1,west_east   $fname_netcdf
	 ncecat   -O -h $fname_netcdf      $fname_netcdf # Add a record variable, default name is "record"
	 ncrename -O -d record,Time        $fname_netcdf # Change name of "record" to "Time"

	 # Finally, stuff the field into the diag file
	 ncks -A $fname_netcdf ./diags_d02_${init}_mem_${mem}_f${f3}.nc

      end # loop over members

   end # loop over forecast hours

   # Go to next initialization
   set init = `$time_exec $init $inc_init`

end # loop over inits

exit
