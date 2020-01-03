#!/bin/csh

set time_exec =    /glade/u/home/schwartz/utils/da_advance_time.exe
set wgrib2_ip =    /glade/work/wrfrt/rt_ensemble_code/grib2/wgrib2/wgrib2 # This was compiled on yellowstone but works on cheyenne.

set start_init = 2017070323     # analysis times to process
set final_init = $start_init #2017070401
set inc_init = 1

set storage_top = . #$PWD/storage
set workdir = $PWD

set file_dir_urma   = /gpfs/fs1/p/mmm/parc/schwartz/jfsp/URMA     # Where files are located
set file_dir_precip = /gpfs/fs1/p/mmm/parc/sobash/MRMS
set file_dir_refl   = $file_dir_precip

set wgrib2_ip_gds = "lambert:-101.0:32.0:46.0:32.0 -120.811:1580:3000.0 23.159:985:3000.0"   # NCAR 3-km ensemble grid

# ------------------------------------------------------------------------------------------------
#### Don't need much changing below here other than pointing to the correct directory strucutre
# ------------------------------------------------------------------------------------------------

mkdir -p $workdir
cd $workdir

#-----------

set init = $start_init
while ( $init <= $final_init )

   set yyyy  = `echo "$init" | cut -c 1-4`
   set sinit = `echo "$init" | cut -c 1-8`
   set ii    = `echo "$init" | cut -c 9-10`

   set storage = ${storage_top}
   mkdir -p $storage

   # set filename of urma files to interpolate
   # these are full paths.
   # variables will be renamed to point to local paths
   set fname_urma   = ${file_dir_urma}/${init}/urma2p5.t${ii}z.2dvaranl_ndfd.grb2
   set fname_precip = ${file_dir_precip}/${sinit}/MRMS_GaugeCorr_QPE_01H_00.00_${sinit}-${ii}0000.grib2.gz
   set fname_refl   = `ls ${file_dir_refl}/${sinit}/MRMS_MergedReflectivityQCComposite_00.50_${sinit}-${ii}*.grib2.gz`

   # copy urma file to working directory
   cp $fname_urma . 
   set fname_urma = `basename $fname_urma`

   # copy precip file to working directory
   set base = `basename $fname_precip`
   cp $fname_precip ./${base}
   gunzip -f $base
   set fname_precip = `basename $fname_precip .gz`

   # copy refl file to working directory
   set base = `basename $fname_refl`
   cp $fname_refl ./${base}
   gunzip -f $base
   set fname_refl = `basename $fname_refl .gz`

   @ err = 0 
   if ( -e $fname_urma && -e $fname_precip && -e $fname_refl ) then

      echo "combining $fname_urma \n $fname_precip \n $fname_refl "

      # bilinear interpolation for everything except precip
      set outname_urma = interp_urma_met_${init}.grb2
      $wgrib2_ip ${fname_urma} -new_grid_winds grid -new_grid_interpolation bilinear -set_grib_type same -new_grid ${wgrib2_ip_gds}  ${outname_urma}
      if ( $status != 0 ) then
	 @ err ++
	 echo "error interpolating ${fname_urma}"
      endif

      # budget interpolation for precip
      set outname_precip = interp_urma_precip_${init}.grb2
      $wgrib2_ip ${fname_precip} -new_grid_winds grid -new_grid_interpolation budget -set_grib_type same -new_grid ${wgrib2_ip_gds}  ${outname_precip}
      if ( $status != 0 ) then
	 @ err ++
	 echo "error interpolating ${fname_precip}"
      endif

      # bilinear interpolation for reflectivity
      set outname_refl = interp_mrms_cref_${init}.grb2
      $wgrib2_ip ${fname_refl} -new_grid_winds grid -new_grid_interpolation bilinear -set_grib_type same -new_grid ${wgrib2_ip_gds}  ${outname_refl}
      if ( $status != 0 ) then
	 @ err ++
	 echo "error interpolating ${fname_refl}"
      endif

      # merge all the files together
      set final_outname = ${storage}/interp_urma_all_${init}.grb2
      rm -f $final_outname
      if ( -e $outname_urma && -e $outname_precip && -e $outname_refl && $err == 0 ) then
	 cat $outname_urma $outname_precip $outname_refl > $final_outname

	 # clean-up
	 rm -f $fname_urma $fname_precip $fname_refl
	 rm -f $outname_urma $outname_precip $outname_refl
      else
	 echo "Error for ${init}"
      endif
   else
      echo "Missing files for ${init}"
   endif

   # Go to next initialization
   set init = `$time_exec $init $inc_init`

end # loop over inits

exit
