/; This code is part of the dynamico set of dynamical adjustment codes
--- Program dynamico_model.ncl
--- This code performs dynamical adjustment as described in Deser et al. 2016 J.Climate.
--- This code reads data from one or multiple perturbed simulations (e.g a simulation constrained by 
--- historical external forcings) as well as data from a reference simulation 
--- (e.g a pre-industrial control simulation with constant forcings).
--- The data are usually an atmospheric circulation variable (such as sea level pressure 
--- slp) and a surface variable (such as temperature -sat- or precipitation -pr-).
--- It then searches the closest N_a analogues of every slp monthly means from the 
--- historical simulation among the slp monthly means from the reference. 
--- It then draws N_b analogues among the N_a ones and uses the moore-penrose 
--- pseudo-inverse to get the best linear combination of the N_b analogues 
--- that fits the historical slp monthly mean. The same coefficients are 
--- then applied to sat or pr or other variable.
--- The algo is done niter times for each month of the historical period, 
--- then the mean over the niter iterations is performed.
--- The mean is then saved as the dynamical contribution for the variable
--- of interest. All iterations are also saved to handle uncertainties.
---
--- This version of the code performs dynamical adjustment for CMIP6 historical simulations using picontrol simulations as
--- reference simulations. Two models are currently included but of course the user can add any other model.
---
--- ASSUMPTIONS: 1. All data files (psl, tas)  must have the same order for lon (0->360)
---                 as the code assumes it is the case (use cdo or nco tricks to do it beforehand if needed)
---                 For latitude, the code assumes that the latitudes go as (-90->90). If not the case, use ncl ::-1 trick
---              2. The code assumes that all data files have correct time coordinates (can be read by cd_calendar)
---              3. The code assumes that simulations have no missing values so no infilling is performed 
---              4. The user can choose to detrend the data from the reference simulation (to remove drift for instance)
---              5. Data filenames are read and written out with specific names reflecting our own data structure 
---                 and organization. It is very likely that the user will have to change the file and directory names. 
---                 (Look at Input and Auxilary Input sections of the code to see or change required names)
---
--- Author: Laurent Terray 08/10/2014
--- 09/10/2014: Original code
--- 10/01/2021: Polished version put on github
----------------------------------------------------------------------------------------------------------------------
;/
;--- Modules
external TEWELESF "~/codes/dynamical_adjustment/tewelesf.so"
external DISTANCEF "~/codes/dynamical_adjustment/distancef.so"
;---
begin
;---------------------------------------------------------------------------------------------------------------------
;--- Section 1: Main Input parameters
      DIST = "EUCLIDE"          ; type of distance for analog search: TEWELES or EUCLIDE
      N_a  = 60                ; number of closest analogs = N_a
      N_b  = 30                ; Select N_b analogs among the N_a
      niter = 10               ; number of iterations
      varp = "psl"              ; atmospheric circulation variable
      var ="tas"                ; variable name to be dynamically adjusted
      scale = .01               ; scale psl (back to hPa)
      scale_t = 1.              ; scale the variable to be adjusted (useful forr precipitation to get mm/month)
      nmem = 2                  ; number of perturbed simulations
      ini_m = 0                 ; initial simulation number in 0,1,2,3 ....,N
      fin_m = 1                 ; final simulation number in 0,...28,29
      MODEL = "CNRM-CM6-1"      ; Model name
      DTRD = "linear"           ; detrending performed on reference data (none, linear, quadratic)
      EXPnm = "historical"      ; type of perturbed simulations: historical, rcp's , ssp's, pacemaker's, others ... 
      REFnm = "piControl"       ; name of the reference simulation
      region = "EU"             ; regional domain
      diro = "/data/scratch/globc/terray/codes/dynadj/results/"   ; Output directory for results
      FMTnm = "NetCDF4Classic"  ; format of output files; all netcdf ncl formats can be used
;
;----- End of Main Input -------------------------------------------------------------------------------------------
;
;----- Additional Input
;----- Section 2: the user can add new geographical domains
;
;--- Regions for slp analogue search: if psl longitude domain crosses 0 then flip must equal "yes"  
      if (region .eq. "EU") then
        flip = "yes"
        minlat_psl                = 25
        maxlat_psl                = 90
        minlon_psl                = -35
        maxlon_psl                = 50
        minlat_tas                = 35
        maxlat_tas                = 75
        minlon_tas                = -15
        maxlon_tas                = 30
      end if
      if (region .eq. "NH") then
        flip = "no"
        minlat_psl                = 25
        maxlat_psl                = 90
        minlon_psl                = 0
        maxlon_psl                = 359.75
        minlat_tas                = 35
        maxlat_tas                = 75
        minlon_tas                = -15
        maxlon_tas                = 30
      end if
;--------------------------------------------------------------------------------------------------------
;----- Section 3: the user can add here more reference simulations
;--- Reference simulation: here two models having piControl simulations from CMIP6
      if (MODEL .eq. "UKESM1-0-LL") then
        tStart_ti          = 196001     ;--- Ref. simulation initial date for filename
        tEnd_ti            = 383912     ;--- Ref. simulation final date for filename
        tStartc            = 196001     ;--- Initial date for reference simulation (needed for cd_calendar if date like 0XXX01)
        tEndc              = 383912     ;--- Final date for reference simulation (as above)
        tStart_t           = 185001     ;--- Pert. simulation initial date for filename
        tEnd_t             = 201412     ;--- Pert. simulation final date for filename
        tStart             = 190001     ;--- Initial date for perturbed simulation
        tEnd               = 201412     ;--- Final date for perturbed simulation
      end if
      if (MODEL .eq. "CNRM-CM6-1") then
        tStart_ti          = 185001 
        tEnd_ti            = 284912 
        tStartc            = 185001
        tEndc              = 284912
        tStart_t           = 185001
        tEnd_t             = 201412
        tStart             = 196501
        tEnd               = 200012
      end if
;----- Load data from reference simulations
      tStart_ti_str = sprinti("%0.6i", tStart_ti) 
      tEnd_ti_str = sprinti("%0.6i", tEnd_ti)
      if (MODEL .eq. "UKESM1-0-LL") then 
        rprefix="/data/scratch/globc/terray/cmip6/picntrl/"
        f_rpsl  = addfile(rprefix+varp+"_Amon_"+MODEL+"_"+REFnm+"_r1i1p1f2_gn_"+tStart_ti_str+"-"+tEnd_ti_str+".nc","r")
        f_rtas  = addfile(rprefix+var+"_Amon_"+MODEL+"_"+REFnm+"_r1i1p1f2_gn_"+tStart_ti_str+"-"+tEnd_ti_str+".nc","r")
      end if
      if (MODEL .eq. "CNRM-CM6-1") then 
        rprefix="/data/scratch/globc/terray/cmip6/picntrl/"
        f_rpsl  = addfile(rprefix+varp+"_Amon_"+MODEL+"_"+REFnm+"_r1i1p1f2_gr_"+tStart_ti_str+"-"+tEnd_ti_str+".nc","r")
        f_rtas  = addfile(rprefix+var+"_Amon_"+MODEL+"_"+REFnm+"_r1i1p1f2_gr_"+tStart_ti_str+"-"+tEnd_ti_str+".nc","r")
      end if  
   ;--------------------------------------------------
   ;------ Loading variables on training period
   ;--------------------------------------------------
      print("   Loading circulation variable "+varp+" and the other variable "+var+" on piCTRL training period")
      timp = f_rpsl->time
      cldr = cd_calendar(timp, -1)
      print("   Initial date for reference simulation = "+cldr(0) )
      irStart_t = ind(cldr .eq. tStartc)
      irEnd_t   = ind(cldr .eq. tEndc)
      ntime_t = irEnd_t - irStart_t + 1
      print("   Total number of months = "+ntime_t)
      delete([/ timp, cldr /])
;--- Circulation variable and adjusted variable for training period from reference simulation
      tmp                    = f_rpsl->$varp$
;--- Check if longitudes and latitudes are correctly ranked
      lon_chk = f_rpsl->lon
      lat_chk = f_rpsl->lat
      nlat_chk = dimsizes(lat_chk)
      if (min(lon_chk) .lt. 0) then
        print("   Some longitudes are negative : "+min(lon_chk) )
        print("   Longitudes must be within 0 - 360; The code stops")
        exit
      end if
      if (lat_chk(0) .gt. lat_chk(nlat_chk-1)) then
        print("   The latitude array goes from North to South")
        print("   The code requires the opposite; the code stops")
        exit
      end if
      delete([/lon_chk, lat_chk, nlat_chk /])
      tmp_t = f_rtas->$var$
      tmp = (/ tmp * scale /)
      tmp_t = (/ tmp_t * scale_t /)
      print("   Min/Max of circulation variable "+varp+" = "+min(tmp)+" - "+max(tmp) )
      print("   Min/Max of the other variable "+var+" = "+min(tmp_t)+" - "+max(tmp_t) )
      if (flip .eq. "yes") then
        tmp = lonFlip(tmp)
        tmp_t = lonFlip(tmp_t)
      end if
      if (DTRD .eq. "linear") then
        print("   Linear detrending performed on reference simulation")
        nyrs = ntime_t/12
        do imo = 0,11
          tmp(imo::12,:,:) = (/ dtrend_msg_n(ispan(1,nyrs,1),tmp(imo::12,:,:), False, False, 0)  /)
          tmp_t(imo::12,:,:) = (/ dtrend_msg_n(ispan(1,nyrs,1),tmp_t(imo::12,:,:), False, False, 0)  /)
        end do
        print("   Min/Max after detrending of circulation variable "+varp+" = "+min(tmp)+" - "+max(tmp) )
        print("   Min/Max after detrending of the other variable "+var+" = "+min(tmp_t)+" - "+max(tmp_t) )
      elseif (DTRD .eq. "quadratic") then
        print("   Quadratic detrending performed on reference simulation")
        do imo = 0,11
          tmp(imo::12,:,:) = (/ dtrend_quadratic_msg_n(tmp(imo::12,:,:), False, False, 0)  /)
          tmp_t(imo::12,:,:) = (/ dtrend_quadratic_msg_n(tmp_t(imo::12,:,:), False, False, 0)  /)
        end do
        print("   Min/Max after detrending of circulation variable "+varp+" = "+min(tmp)+" - "+max(tmp) )
        print("   Min/Max after detrending of the other variable "+var+" = "+min(tmp_t)+" - "+max(tmp_t) )
      elseif (DTRD .eq. "none")
        print("   No detrending on reference simulation")
      else
        print("   detrending not implemented; Stop the code")
        exit
      end if
      rpsl = tmp(time|irStart_t:irEnd_t, {lat|minlat_psl:maxlat_psl},{lon|minlon_psl:maxlon_psl})
      rtas = tmp_t(time|irStart_t:irEnd_t, {lat|minlat_tas:maxlat_tas},{lon|minlon_tas:maxlon_tas})
      delete([/ tmp, tmp_t, f_rtas, f_rpsl/])
;--- Getting latitudes & longitudes for reshape and weights
      lat              = dble2flt(rpsl&lat)
      lon              = dble2flt(rpsl&lon)
      nlat             = dimsizes(lat)
      nlon             = dimsizes(lon)
      lat_tas          = dble2flt(rtas&lat)
      lon_tas          = dble2flt(rtas&lon)
      nlatt            = dimsizes(lat_tas)
      nlont            = dimsizes(lon_tas)
;---- Set new seed for the random number generator
      random_setallseed(N_a*N_b*niter, N_a*tEnd_ti)
;--- Members
      member = ispan(0,nmem-1,1)
      member!0 ="member"
      miss_val = 1.e+20
   ;---------------------------------------------------------------------
   ;----- Load now perturbed simulations and loop over number of members
   ;---------------------------------------------------------------------
      tStart_t_str = sprinti("%0.6i", tStart_t) 
      tEnd_t_str = sprinti("%0.6i", tEnd_t) 
      do i_m = ini_m,fin_m
        i_mn = i_m + 1
        print("   Member number "+i_mn)
        if (MODEL .eq. "UKESM1-0-LL") then
          eprefix="/data/scratch/globc/terray/cmip6/"+EXPnm+"/"
          f_epsl  = addfile(eprefix+varp+"_Amon_"+MODEL+"_"+EXPnm+"_r"+i_mn+"i1p1f2_gn_"+tStart_t_str+"-"+tEnd_t_str+".nc","r")
        end if
        if (MODEL .eq. "CNRM-CM6-1") then
          eprefix="/data/scratch/globc/terray/cmip6/"+EXPnm+"/"
          f_epsl  = addfile(eprefix+varp+"_Amon_"+MODEL+"_"+EXPnm+"_r"+i_mn+"i1p1f2_gr_"+tStart_t_str+"-"+tEnd_t_str+".nc","r")
        end if
        tmp  = f_epsl->$varp$
        if (flip .eq. "yes") then
          tmp = lonFlip(tmp)
        end if
        if (i_m .eq. ini_m) then
          timp                    = f_epsl->time
          cldr                   = cd_calendar(timp, -1)
          print("   Initial date for perturbed simulation = "+cldr(0) )
          ieStart_f              = ind(cldr .eq. tStart)
          ieEnd_f                = ind(cldr .eq. tEnd)
          ntime_f                = ieEnd_f - ieStart_f + 1
          print("   Start of the period to be dunamically adjusted = "+cldr(ieStart_f))
          print("   End of the period to be dynamically adjusted   = "+cldr(ieEnd_f))
          print("   Total number of months to be adjusted = "+ntime_f)
          delete([/ timp, cldr /])
        end if
        delete([/ f_epsl/])
;
        fpsl= tmp(time|ieStart_f:ieEnd_f, {lat|minlat_psl:maxlat_psl},{lon|minlon_psl:maxlon_psl} )
        epsl= tmp(time|ieStart_f:ieEnd_f, {lat|minlat_tas:maxlat_tas},{lon|minlon_tas:maxlon_tas} )
        epsl = (/ epsl * scale /)  
        fpsl = (/ fpsl * scale /) 
        if(.not. isatt(fpsl,"FillValue")) then
          fpsl@_FillValue = miss_val
        end if
        delete([/ tmp/])
;
        if (i_m .eq. ini_m) then
          epsl_GM = dim_avg_n_Wrap( epsl,0)
          fpsl_GM = dim_avg_n_Wrap( fpsl,0)
;--- Create array to hold dynamically adjusted variable
          etas_adj=conform_dims( (/nmem,niter,ntime_f,nlatt,nlont/), epsl, (/2,3,4/) )
          etas_adj!0 = "member"
          etas_adj&member = member
          itera=ispan(1,niter,1)
          itera!0 = "iteration"
          etas_adj!1 = "iteration"
          etas_adj&iteration = itera
          etas_adj!2 = "time"
          etas_adj&time = epsl&time
          etas_adj!3 = "lat"
          etas_adj&lat = lat_tas
          etas_adj!4 = "lon"
          etas_adj&lon = lon_tas
;--- Circulation varibale now
          epsl_rec = conform_dims( (/nmem,niter,ntime_f,nlat,nlon/), fpsl, (/2,3,4/) )
          epsl_rec!0 = "member"
          epsl_rec&member = member
          epsl_rec!1 = "iteration"
          epsl_rec&iteration = itera
          epsl_rec!2 = "time"
          epsl_rec&time = epsl&time
          epsl_rec!3 = "lat"
          epsl_rec&lat = lat
          epsl_rec!4 = "lon"
          epsl_rec&lon = lon
;--- Add attributes
          today = systemfunc("date")
          etas_adj@long_name = rtas@long_name+" Dynamical contribution for all iterations"
          epsl_rec@long_name = "Constructed "+rpsl@long_name+" analog for all iterations"
          epsl_rec@description = " Dynamico run performed on "+today(0)+" with algorithm parameters: N_a = "+N_a+" , N_b = "+N_b+" , niter = "+niter
          etas_adj@description = epsl_rec@description
          time = epsl&time
      end if
;
;---- Calculate distances
;
;---- First define needed arrays      
;
      ntt = ntime_t/12
      ntf = ntime_f/12
      print("   Number of years in Reference period = "+ntt)
      print("   Number of years in Perturbed period = "+ntf)
      dist_m = new((/ntt,ntf,12/), float,miss_val)
;
       print("   CPU time: " + get_cpu_time())
       print("   Distance calculation Start ")
      if (DIST .eq. "EUCLIDE") then 
;----- Euclidian distance calculation
        print("   Calculating euclidian distance")
; reorder dimensions
        rpsl_2d = onedtond(ndtooned(rpsl),(/ntime_t,nlat*nlon/))
        fpsl_2d = onedtond(ndtooned(fpsl),(/ntime_f,nlat*nlon/)) 
        do im = 0, 11
          DISTANCEF::distancef(rpsl_2d(im::12,:), fpsl_2d(im::12,:), ntt, ntf, nlat*nlon, dist_m(:,:,im) )
        end do
        delete([/rpsl_2d,fpsl_2d/])
      end if
      if (DIST .eq. "TEWELES") then
 ;----- Teweles and Wobus score calculation
        print("   Calculating Teweles-Wobus score")
        do im = 0, 11
          TEWELESF::tewelesf(rpsl(im::12,:,:), fpsl(im::12,:,:), lat, lon, ntt, ntf, nlat,nlon,dist_m(:,:,im) )
         end do
       end if
      print("   CPU time: " + get_cpu_time())
      print("   Distances done ")
;------------------------------------------------------------------------------
; rank Distances 
      i_dst = dim_pqsort_n(dist_m,2,0) ; sort dist increasing order
      delete(dist_m)
      do im = 0, 11
        i_dst(:,:,im) = (/ i_dst(:,:,im)*12 + im /); reconstruct orig. address
      end do
;
      dpsl_z = conform_dims( (/N_a*ntime_f,nlat,nlon/),fpsl_GM,(/1,2/) ) 
      dtas_z = conform_dims( (/N_a*ntime_f,nlatt,nlont/),epsl_GM,(/1,2/) )
      dst_1d = ndtooned(i_dst(0:N_a - 1,:,:))
      dpsl_z = (/ rpsl(dst_1d,:,:) /)
      dtas_z = (/ rtas(dst_1d,:,:) /)
      dpsl = onedtond(ndtooned(dpsl_z),(/N_a,ntime_f,nlat,nlon /) )
      dtas = onedtond(ndtooned(dtas_z),(/N_a,ntime_f,nlatt,nlont /) )
      delete([/dst_1d,i_dst,dpsl_z,dtas_z/])
;-------------------------------------------------------------------------------
;--- Analog step
;--- Reconstruct monthly analogs using Moore-Penrose inverse
;--- 1. Loop on all perturbed simulation members 
;--- 2. Loop on all months
;--- 3. Loop on iteration of analog steps
;
      fpsl_m = onedtond(ndtooned(fpsl),(/ntime_f,nlat*nlon/))
      delete(fpsl)
      do ii = 0, ntime_f - 1      ; Start of time loop
        print("   time step = "+ii)
;--- First loop on sets of analogs :   dims of dpsl N_a,ntime_f,nlat,nlon
        rec_tas_2d = new((/niter,nlatt,nlont/), float, miss_val)
        rec_psl_2d = new((/niter,nlat,nlon/), float, miss_val)
        do ia = 0,niter-1         ; Start of iteration loop
;--- 1- select randomly N_b analogs among the N_a ones
          nba = generate_unique_indices(N_a)
          iadr = nba(0:N_b-1)
          A_mat = onedtond(ndtooned(dpsl(iadr,ii,:,:)),(/N_b,nlat*nlon/))
          dimss = dimsizes(A_mat)
          nca = dimss(0)
          nra = dimss(1)
          delete(dimss)
          s = new((/nca,nca/),float)
          s = (/ 0. /)
          u   = new ( (/nca,nra/) , typeof(A_mat) )  ; pre-allocate space 
          v   = new ( (/nca,nca/) , typeof(A_mat) )
          sgv  = svd_lapack (A_mat, "S" , "S", 1,  u, v)  ; if 0 v is here v_T
          delete(A_mat)
          do im = 0, nca - 1
            s(im,im) = (/ 1./sgv(im) /)
          end do
;--- Get pseudo-inverse
          x = fpsl_m(ii,:) # transpose(u) # s # v
          delete([/ s, u, v, sgv, nba/])
          rec_tas = x # onedtond(ndtooned(dtas(iadr,ii,:,:)),(/N_b,nlatt*nlont/))
          rec_psl = x # onedtond(ndtooned(dpsl(iadr,ii,:,:)),(/N_b,nlat*nlon/))
          delete([/ x, iadr/])
          rec_tas_2d(ia,:,:) = onedtond(rec_tas,(/nlatt,nlont/))
          rec_psl_2d(ia,:,:) = onedtond(rec_psl,(/nlat,nlon/))
          delete([/ rec_tas, rec_psl/])
        end do               ; end of loop on number of iteration on set of analogs
        etas_adj(i_m,:,ii,:,:) = (/ rec_tas_2d /)
        epsl_rec(i_m,:,ii,:,:) = (/ rec_psl_2d /)
        delete([/ rec_tas_2d, rec_psl_2d /])
      end do                 ; End of time loop
      delete(fpsl_m) 
;-----
      delete([/dpsl, dtas /])
      print("   CPU time: " + get_cpu_time())
      print("   Analog reconstruction done for member = "+(i_m+1))
;
;--- End of loop over members
;
    end do
;--- Average over niter iterations
    etas_adj_avg = dim_avg_n_Wrap(etas_adj,1)
    epsl_rec_avg = dim_avg_n_Wrap(epsl_rec,1)
    print("   Write Output files in folder "+diro)
;
;--- Write output data
;--- Opening output files
    setfileoption("nc", "Format",  FMTnm)
    filo = var+"_"+MODEL+"_"+EXPnm+"_"+DIST+"_"+N_a+"_"+N_b+"_"+niter+"_region_"+region+"_member_range_"+ini_m+"_"+fin_m+".nc"
    filp = varp+"_"+MODEL+"_"+EXPnm+"_"+DIST+"_"+N_a+"_"+N_b+"_"+niter+"_region_"+region+"_member_range_"+ini_m+"_"+fin_m+".nc"
    system("/bin/rm -f " + diro + filo)
    system("/bin/rm -f " + diro + filp)
    fout  = addfile (diro + filo, "c")
    foup  = addfile (diro + filp, "c")
    setfileoption(fout,"DefineMode",True)
    setfileoption(foup,"DefineMode",True)
;--- file for adjusted variable first
    dimNames = (/"member", "iteration","lat", "lon", "time"/)  
    dimSizes = (/ nmem , niter, nlatt,  nlont, ntime_f/) 
    dimUnlim = (/ False, False, False, False, True/)   
    filedimdef(fout,dimNames,dimSizes,dimUnlim)
    filevardef(fout, "member"  ,typeof(member),getvardims(member))
    filevardef(fout, "iteration",typeof(itera),getvardims(itera))
    filevardef(fout, "lat"  ,typeof(lat_tas),getvardims(lat_tas))
    filevardef(fout, "lon"  ,typeof(lon_tas),getvardims(lon_tas))
    filevardef(fout, "time"  ,typeof(time),getvardims(time))
    filevardef(fout, "tas_all"    ,typeof(etas_adj)  ,getvardims(etas_adj))
    filevardef(fout, "tas_adj"    ,typeof(etas_adj_avg)  ,getvardims(etas_adj_avg))
    filevarattdef(fout,"lat"  ,lat_tas)
    filevarattdef(fout,"lon"  ,lon_tas)
    filevarattdef(fout,"time"  ,time)
    filevarattdef(fout,"tas_adj"  ,etas_adj_avg)
    filevarattdef(fout,"tas_all"  ,etas_adj)
    setfileoption(fout,"DefineMode",False)
    fout->member   = (/member/)
    fout->iteration = (/itera/)
    fout->lat    = (/lat_tas/)
    fout->lon    = (/lon_tas/)
    fout->time   = (/time/)
    fout->tas_adj     = (/etas_adj_avg/)
    fout->tas_all     = (/etas_adj/)
    delete(fout)
;--- Circulation variable now
    dimSizes = (/ nmem , niter, nlat,  nlon, ntime_f/) 
    dimUnlim = (/ False, False, False, False, True/)   
    filedimdef(foup,dimNames,dimSizes,dimUnlim)
    filevardef(foup, "member"  ,typeof(member),getvardims(member))
    filevardef(foup, "iteration",typeof(itera),getvardims(itera))
    filevardef(foup, "lat"  ,typeof(lat),getvardims(lat))
    filevardef(foup, "lon"  ,typeof(lon),getvardims(lon))
    filevardef(foup, "time"  ,typeof(time),getvardims(time))
    filevardef(foup, "psl_all"    ,typeof(epsl_rec)  ,getvardims(epsl_rec))
    filevardef(foup, "psl_rec"    ,typeof(epsl_rec_avg)  ,getvardims(epsl_rec_avg))
    filevarattdef(foup,"lat"  ,lat)
    filevarattdef(foup,"lon"  ,lon)
    filevarattdef(foup,"time"  ,time)
    filevarattdef(foup,"psl_rec"  ,epsl_rec_avg)
    filevarattdef(foup,"psl_all"  ,epsl_rec)
    setfileoption(foup,"DefineMode",False)
    foup->member   = (/member/)
    foup->iteration = (/itera/)
    foup->lat    = (/lat/)
    foup->lon    = (/lon/)
    foup->time   = (/time/)
    foup->psl_rec     = (/epsl_rec_avg/)
    foup->psl_all     = (/epsl_rec/)
    delete(foup)
;********************************************************************************************
end
