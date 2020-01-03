;
;  Script to make panel-plots of WRF sfc-fields.  Paneling allows for easy
;   visual comparison of multiple experiments
; 
; Need to set paths to the data files to plot, the advance_time script, and map_info_file below
; Need to also set the start time and increment hours, and latitude longitude bounds below.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
 
begin
system("date")

;;;;;;;;;;;;;;;;;;;;;;;

EXEC="/glade/work/ncar_gft/advance_time.csh" ; script in gft source directory to advance date/time

; path to directories holding files of interest
parent_dir = "/glade/scratch/schwartz/jfsp/"
expt1 = parent_dir
; change the following line and the plot_dims line to set the number of experiments to plot
; currently set to plot 3 experiments
experiments = (/ expt1, expt1, expt1  /) ; , expt1, expt1, expt1 /)   ; this is what is ultimately used to point to directories
expt_prefix = "mem"+ispan(1,3,1)
expt_labels = expt_prefix

start_init = "2017070300"   ; forecast initializations to process
final_init = start_init 
inc_init   = 24

fhrs = (/24/) ; ispan(43,58,1)  ; forecast hours you want to process

variables = (/ "binary_gust_front" ,"tot_precip", "FRONTtot", "VEC_WIND_DIFF_MAG", "gradtheta" /) ; , "GradFRONTtot" /)  ; fields you want to plot
legend_labels = variables  ; labels for the variables

; only 2-d variables can be plotted (okay if a time dimension...code will check for that)
variables_to_contour_only = (/ "binary_gust_front" /)    ; these variables will not be colored filled--just drawn with contours

; overlay wind barbs?
overlay_winds = True   ; True or False (capital T or F)

; overlay binary gust front?
overlay_binary_gust_front = True ; True or False (capital T or F)

; a netcdf file with lat/lons of the grid of the model/data you're plotting
map_info_file = "/glade/work/ncar_gftp/con/ens_mapfac.nc   ; ncar ensemble
;map_info_file = "/glade/work/ncar_gftp/con/hrrr_geo.nc    ; ncep HRRR

user_specified_lat_lon_corners = True ; True to show domain over area specified below...lambert lats/lons controlled by "map_info_from_file"
;  Washington DC
;  coords_SW = (/ 38,-80  /)  ; (lat,lon) of SW corner
;  coords_NE = (/ 40,-76  /)      ; (lat,lon) of NE corner

   ; Rockies
;  coords_SW = (/ 38,-110 /)  ; (lat,lon) of SW corner
;  coords_NE = (/ 44,-100 /)      ; (lat,lon) of NE corner
   
   ; NM 2017070300_f024
;  coords_SW = (/ 30,-110  /)  ; (lat,lon) of SW corner
;  coords_NE = (/ 35,-105  /)      ; (lat,lon) of NE corner

   ; TX panhandle 2017070300_f024
   coords_SW = (/ 30,-105  /)  ; (lat,lon) of SW corner
   coords_NE = (/ 35,-100  /)      ; (lat,lon) of NE corner

storage = "./"  ; where output will be located

plot_dims = (/1,3/)  ; refers to number of panels in a plot...typically based on number of experiments

;;;;;;;;;;;;;;;;;;;;;;;;;;

system("mkdir -p "+storage)

nexpts = dimsizes(experiments)
nvars = dimsizes(variables)
nfhrs = dimsizes(fhrs)

min_val = new(nvars,"float",-999999)
max_val = new(nvars,"float",-999999)

plots = new( (/nvars,nexpts/),"graphic",-999999)

if ( fileexists(map_info_file)) then
   m = addfile(map_info_file+".nc","r")
   lats = m->XLAT(0,:,:)
   lons = m->XLONG(0,:,:)
else
   print(map_info_file+" not there. try another file for map information.")
   exit
end if

mpres = True
mpres@gsnMaximize = True
wrf_mapres_c(m,mpres,0)
mpres@mpOutlineBoundarySets = "GeophysicalAndUSStates" ; "AllBoundaries" ; "GeophysicalAndUSStates"
mpres@gsnFrame = False
mpres@gsnDraw = False
mpres@gsnSpreadColors = False ; True
mpres@gsnSpreadColorStart = 22
mpres@gsnSpreadColorEnd = -2
mpres@lbLabelBarOn = False ; turn off individual colorbars
mpres@mpGridAndLimbOn   = True   ; lat/lon lines
mpres@mpGridSpacingF = 10.0       ; for equal lat/lon spacing
;res@mpGridLatSpacingF = 10.0
;res@mpGridLonSpacingF = 10.0
mpres@mpGridLineColor   = "Grey40"
mpres@mpLimbLineColor   = mpres@mpGridLineColor
mpres@mpGeophysicalLineColor = mpres@mpGridLineColor
mpres@mpNationalLineColor = mpres@mpGridLineColor
mpres@mpPerimLineColor = mpres@mpGridLineColor
mpres@mpUSStateLineColor = mpres@mpGridLineColor
mpres@mpGridLineDashPattern  = 1
mpres@mpGridLineThicknessF   = 0.5
mpres@mpLimbLineThicknessF   = 0.5
mpres@pmTickMarkDisplayMode = "Always"
;mpres@cnLevelSelectionMode = "ManualLevels"
mpres@cnLevelSelectionMode = "ExplicitLevels"
mpres@tmXTOn = False
mpres@tmYROn = False
mpres@cnLineLabelInterval = 2  ; Label every line
mpres@cnLineLabelFontHeightF = .018
;mpres@cnLineLabelDensityF = 1.5 ; Can add more labels
mpres@cnLineLabelPlacementMode = "constant"
mpres@cnFillMode = "RasterFill"
mpres@mpFillDrawOrder = "PreDraw"

vecres = True
vecres@gsnDraw  = False
vecres@gsnFrame = False
vecres@vcGlyphStyle = "WindBarb"
vecres@vcWindBarbColor = "Grey70" ; "Black"
vecres@vcWindBarbLineThicknessF = 1.5
vecres@vcMonoWindBarbColor = True
vecres@vcRefAnnoOn = False
vecres@vcRefMagnitudeF = 10.
vecres@vcRefLengthF = 0.030
vecres@vcMinDistanceF = 0.030  ; wind barb density
vecres@vcWindBarbTickLengthF = 0.4
vecres@vcWindBarbTickSpacingF = 0.3
vecres@gsnRightString = ""
vecres@gsnLeftString = ""
;;;;vecres@tfDoNDCOverlay = True

pnlres = True
pnlres@gsnMaximize = True
pnlres@gsnPanelBottom = 0.05  ; add some space on bottom of the image
pnlres@gsnPanelLabelBar = True                                    ; Add panel (common) label bar
pnlres@lbOrientation = "horizontal" ; "vertical"
;pnlres@lbLabelStride = 2
pnlres@lbLabelAutoStride = True
pnlres@pmLabelBarOrthogonalPosF = -0.035  ; move label bar over a bit
;pnlres@lbLabelPosition = "left"
pnlres@lbLabelFontHeightF = .015
pnlres@lbTitleOn = True
;pnlres@lbTitlePosition = "right"
pnlres@lbTitleFontHeightF = .020      ; Change title size
pnlres@txFontHeightF = .015
pnlres@gsnPanelXWhiteSpacePercent = 5
pnlres@gsnPanelYWhiteSpacePercent = 5
;pnlres@vpWidthF =  0.58
;pnlres@vpHeightF =  0.18

;;;;;;;;;;;;;;

ii = 0
init = start_init
do while init .le. final_init ; LOOP OVER INITS

   do f = 0, nfhrs-1  ; LOOP OVER FORECAST HOURS

      fhr = fhrs(f)
      st_fhr = sprinti("%0.3i",fhr)

      valid_time = systemfunc(EXEC+" "+init+" "+fhr)

      outname = storage+"gust_front_variables_"+init+"_f"+st_fhr
      wks=gsn_open_wks("pdf",outname)
     ;gsn_define_colormap(wks,"rainbow+white+gray")
      gsn_define_colormap(wks,"precip3_16lev")

      do v = 0,nvars -1 ; LOOP OVER FIELDS/VARIABLES

	 variable = variables(v)

	 if ( str_lower(variable).eq."fronttot") then
	    contour_levels := ispan(0,100,5) ; should make more dynamic
	 end if
	 if ( str_lower(variable).eq."vec_wind_diff_mag") then
	    contour_levels := ispan(0,30,2) ; should make more dynamic
	 end if
	 if ( str_lower(variable).eq."binary_gust_front") then
	    contour_levels := (/ 0.5 /)
	 end if
	 if ( str_lower(variable).eq."tot_precip") then
	    contour_levels := (/ 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0 /)
	 end if
	 if ( str_lower(variable).eq."gradtheta") then
	    contour_levels := ispan(0,100,5) 
	 end if

	 do k = 0,nexpts-1  ; LOOP OVER EXPERIMENTS

	    expt = experiments(k)
	    mpres@cnLevels := contour_levels

	    ;fname = expt+ PATH TO FILE
	    fname = expt+"gf_tool_"+expt_prefix(k)+"_"+init+"_f"+st_fhr+".nc"

	    if ( fileexists(fname)) then
	       print("processing "+fname)
	    else
	       print(fname+" is missing.")
	       continue
	    end if

	    ; add the file
	    a = addfile(fname,"r")

	    ; get data for the variable
	    if ( isfilevar(a,variable) ) then
	       var_dims := getfilevardimsizes(a,variable)
	       if ( dimsizes(var_dims) .eq. 3 ) then ; assume first dimension is time.
		  data = a->$variable$(0,:,:) ; assumes a 2-D variable
	       else
		  data = a->$variable$(:,:) ; assumes a 2-D variable
	       end if
	    else
	       print(variable+" is not present in "+fname)
	       continue
	    end if

	    if ( overlay_binary_gust_front ) then
	       if ( isfilevar(a,"binary_gust_front") ) then
		  var_dims := getfilevardimsizes(a,"binary_gust_front")
		  if ( dimsizes(var_dims) .eq. 3 ) then ; assume first dimension is time.
		     binary_gf = a->binary_gust_front(0,:,:) ; assumes a 2-D variable
		  else
		     binary_gf = a->binary_gust_front(:,:) ; assumes a 2-D variable
		  end if
	       else
		  print("variable binary_gust_front is not present in "+fname)
		  continue
	       end if
	    end if

	    ; get winds, if requested
	    if ( overlay_winds ) then
	       var_dims := getfilevardimsizes(a,"U10")  ; assume same dimsizes for U10 and V10
	        if ( dimsizes(var_dims) .eq. 3 ) then ; assume first dimension is time.
		   u_plane = a->U10(0,:,:)*1.94386 ; convert m/s into kts
		   v_plane = a->V10(0,:,:)*1.94386 
	        else
		   u_plane = a->U10(:,:)*1.94386 ; convert m/s into kts
		   v_plane = a->V10(:,:)*1.94386 
		end if
	    end if

	    if ( any( str_lower(variables_to_contour_only).eq.str_lower(variable)) ) then
	      mpres@cnFillOn = False
	      mpres@cnLinesOn = True
	      mpres@cnLineLabelsOn = False
	      mpres@cnInfoLabelOn = False
	      mpres@mpFillOn = False
	     ;mpres@mpLandFillColor = "Grey40" ; 192
	    else
	      mpres@cnFillOn = True
	      mpres@cnLinesOn = False
	      mpres@cnLineLabelsOn = False
	      mpres@cnInfoLabelOn = False
	      mpres@mpFillOn = False
	    end if

	    if ( user_specified_lat_lon_corners ) then
		mpres@mpLeftCornerLatF = coords_SW(0)
		mpres@mpLeftCornerLonF = coords_SW(1)
		mpres@mpRightCornerLatF = coords_NE(0)
		mpres@mpRightCornerLonF = coords_NE(1)
	    end if

           ;mpres@gsnRightString = "Max : "+max(data)
            mpres@gsnRightString = variable
	    mpres@gsnLeftString = expt_labels(k)

	    data@lat2d = lats
	    data@lon2d = lons

	    ; define ranges for plotting. keep constant across all inits, fhrs, and experiments
	    ; CURRENTLY NOT USED
	   ;if ( k.eq.0 .and. ii.eq.0 .and. f.eq.0) then
	   ;   min_val(v) = min(data)
	   ;   max_val(v) = max(data)
	   ;end if
	   ;contour_spacing = nice_mnmxintvl(min_val(v),max_val(v),60,False)
	   ;mpres@cnMinLevelValF = contour_spacing(0)
	   ;mpres@cnMaxLevelValF = contour_spacing(1)
	   ;mpres@cnLevelSpacingF = contour_spacing(2)

	    plots(v,k) = gsn_csm_contour_map(wks,data,mpres)
	    delete(data)

	    if ( overlay_winds ) then
	       u_plane@lat2d = lats
	       v_plane@lat2d = lats
	       u_plane@lon2d = lons
	       v_plane@lon2d = lons
	       vector := gsn_csm_vector(wks,u_plane,v_plane,vecres)
	       overlay(plots(v,k),vector)
	       delete(u_plane)
	       delete(v_plane)
	    end if
	    
	    if ( overlay_binary_gust_front ) then
	       mpres@cnFillOn = False
	       mpres@cnLinesOn = True
	       mpres@cnLineLabelsOn = False
	       mpres@cnInfoLabelOn = False
	       mpres@mpFillOn = False
	      ;mpres@mpLandFillColor = "Grey40" ; 192
	       mpres@cnLevels := (/ 0.5 /)
	       binary_gf@lat2d = lats
	       binary_gf@lon2d = lons
	       mpres@gsnRightString = ""
	       mpres@gsnLeftString = ""
	      ;mpres@tfDoNDCOverlay = True
	       gf := gsn_csm_contour(wks,binary_gf,mpres)
	       overlay(plots(v,k),gf)
	      ;delete(gf)
	      ;mpres@tfDoNDCOverlay = False ; reset to default
	    end if

	 end do ; END LOOP OVER EXPERIMENTS (loop over "k")

	 pnlres@txString = st_fhr+" hr fcst initialized "+init+", valid "+valid_time
	 pnlres@lbTitleString = legend_labels(v)
	 gsn_panel(wks,(/plots(v,:)/),plot_dims,pnlres)     ; Now make the panel plot, which advances frame

      end do ; END LOOP OVER FIELDS (loop over "v")

     ;destroy(wks) ; get rid of wks. make new one for next fhr/init

   end do ; END LOOP OVER FHRS (loop over "f")

   init = systemfunc(EXEC+" "+init+" "+inc_init)
   ii = ii + 1

end do ; END LOOP OVER INITS (loop over "ii")

;;;;;;;

system("date")
end
