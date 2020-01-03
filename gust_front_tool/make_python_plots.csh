#!/bin/csh

echo =============  Begin $0 `date` =============

set echo
setenv NCARGFT /Fullpath/to/ncar_gft
set start_init = 2017072500 #2017070300
set final_init = $start_init
set inc_init = 24 # hours

set model = 0   # Type of model output  0 = ncar ensemble  1 = HRRR   2 = URMA

set fhrs = (`seq 21 1 27`) #`seq 24 1 28` #( 18 21 24 27 ) #`seq 23 1 27`

set storage_top = $PWD # where PNGs are stored
set workdir = $PWD # where this script is run

set top_dir = /gpfs/fs1/scratch/schwartz/jfsp # where netCDF files from output of GF tool are located

# Center lat/lon for plot.
# Currently, the "delta" about the center point is 4.0 degrees.
#   To change the "delta", look for variable "delta" in webplot.py
set center_lat = 42.28   #45.49
set center_lon = -123.95 #-115.23

if ( $model == 0 ) then
  setenv ENS_SIZE 10   # Ensemble size
else 
  setenv ENS_SIZE 1    # Ensemble size
endif

set is_urma = False # True or False, capitalize first letter

# -----------------------------------------------------------------------------------------------------------------------
#### Don't need changing below here other than adding/tweaking/removing graphics products and how they're submitted #####
# -----------------------------------------------------------------------------------------------------------------------

mkdir -p $workdir
cd $workdir
setenv EXP_DIR $workdir  # Environmental variable needed for python scripts

if ( $model == 0 ) then
  ln -sf ${NCARGFT}/con/latlon_ens.nc latlon.nc
else ( $model == 1 ) then
  ln -sf ${NCARGFT}/con/latlon_hrrr.nc latlon.nc
else
  echo Need to hardwire file for URMA
  exit (1)
endif

#---

set init = $start_init
while ( $init <= $final_init )

   set yyyy  = `echo "$init" | cut -c 1-4`
   set sinit = `echo "$init" | cut -c 1-8`
   set hh    = `echo "$init" | cut -c 9-10`

   set storage = ${storage_top}/${init} # Graphics will be stored here
   mkdir -p $storage

 # set file_dir = ${top_dir}/${yyyy}/${sinit}
   set file_dir = ${top_dir}

  # loop over fhrs to make plots for individual forecast hours
  # first link the files to the working directory
   foreach fhr ( $fhrs )
      set f3 = `printf %03d $fhr`
      if ( $model == 2 ) then
         ln -sf ${file_dir}/gf_tool_urma_${init}_f${f3}.nc ./gf_tool_mem1_${init}_f${f3}.nc
      else
         ln -sf ${file_dir}/gf_tool_mem*_${init}_f${f3}.nc .
      endif
    #${file_dir}/diags_d02_${init}_mem*f${f3}.nc.gz .
    #gunzip -f diags_d02_${init}_mem*f${f3}.nc.gz

      # make plots
      ./make_webplot.py -d=${init} -tr=${fhr},${fhr} -f=binary-gf_neprob_0.5                  -clat=${center_lat}  -clon=${center_lon} -t='Probability of gust front within 40 km'
      ./make_webplot.py -d=${init} -tr=${fhr},${fhr} -f=binary-gf_paintball_0.5               -clat=${center_lat}  -clon=${center_lon} -t='Gust front paintball plot'
      ./make_webplot.py -d=${init} -tr=${fhr},${fhr} -f=cref_stamp   -c=binary-gf_stamp -b=wind10m_stamp  -clat=${center_lat}  -clon=${center_lon} -t='Comp. refl (dBz), 10-m wind (kts), gust fronts (contour)'
      ./make_webplot.py -d=${init} -tr=${fhr},${fhr} -f=precip_stamp -c=binary-gf_stamp -b=wind10m_stamp  -clat=${center_lat}  -clon=${center_lon} -t='1-h precip (in), 10-m wind (kts), gust fronts (contour)'
   end

   # Move PNG files to storage
   mv ./*f???*.png $storage

   # Clean-up
#  rm -f ./gf_tool_mem*_${init}_f*.nc .
   rm -f ./core*

   # Go on to the next initialization time
   set init = `date +%Y%m%d%H -d "$sinit $hh + ${inc_init} hour"`

end # Loop over initial times

exit
