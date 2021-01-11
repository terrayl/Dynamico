/;
--- Program dynamico.ncl
--- This code performs dynamical adjustment as described in Deser et al. 2016 J.Climate.
--- The code first reads the observed data for an atmospheric circulation variable, 
--- (here, sea level pressure, slp) and for a variable that is going to be dynamically adjusted,
--- (here, surface air temperature, sat). It also defines two geographical domain for slp and sat.
--- It searches N_a closest analogues of slp monthly means by excluding the current 
--- month among the slp monthly means from the whole period. 
--- It then draws N_b analogues among the N_a ones and uses the Moore-Penrose
--- pseudo-inverse to get the best linear combination of the N_b analogues that fits the historical
--- slp monthly data. The same coefficients are then applied to sat data.
--- In case external forcing is present, sat might have to be detrended beforehand and several parametric methods 
--- are available (linear, quadratic, cubic, spline). Another option is to use data files coming
--- from applying non parametric methods (eemd, model estimates from large ensembles, cmip means).
--- The algo is done niter times for each month of the historical period, 
--- then the mean over the niter iterations is performed.
--- The mean is then saved as the dynamical contribution for the variable
--- of interest. All iterations are also saved to handle uncertainties.
--- 
--- ASSUMPTIONS: 1. All data files (psl, tas)  must have the same order for lon (0->360)
---                 as the code assumes it is the case (use cdo or nco tricks to do it beforehand if needed)
---                 For latitude, the code uses ncl ::-1 trick to have the latitudes as (-90->90)
---              2. Assume that all files have correct time coordinates (can be read by cd_calendar)
---              3. The infilling is based on the poisson_grid_fill routine with default parameters (check the ncl page for details)
---              4. The psl array need to be infilled if missing values are present (use FILP = "yes")
---              5. Data filenames are read and written out with specific names reflecting our own data structure and organization.
---                 It is very likely that the user will have to change the file and directory names. 
---                 (Look at Input and Auxilary Input sections of the code to see or change required names)
---
--- Author: Laurent Terray 17/09/2017 CECI, CNRS/Cerfacs
-------------------------------------------------------------------------------------------------------------
;/
;--- Modules
external TEWELESF "~/codes/dynamical_adjustment/tewelesf.so"
external DISTANCEF "~/codes/dynamical_adjustment/distancef.so"
;---
begin
;--- Main Input parameters
;--- Algorithm
      DIST = "EUCLIDE"                     ; distance used to get analogues EUCLIDE or TEWELES
      N_a  = 70                            ; Number of closest analogs
      N_b  = 50                            ; Select N_b analogs among N_a
      niter = 30                           ; Number of algorithm iterations
      region = "EU"                        ; geographical domain
      FILL = "no"                          ; Infilling or not adjusted variable
      FILP = "no"                          ; Infilling or not circulation variable
      DEBUG = True                         ; Printing checks
;--- Input/output
      ROOT = "/data/scratch/globc/terray/" ; root directory
      DIRO = ROOT+"obs/"                   ; directory for observed datasets 
      DIRR = ROOT+"codes/dynadj/results/"  ; directory to save output file
      FMTnm = "NetCDF4Classic"             ; format of output files
      ANAP = "20CR_V3"                     ; Observed circulation variable dataset
      ANAT = "HadCRUT5"                    ; Observed to-be adjusted variable dataset
      var ="tas"                           ; adjusted variable name for dynamico input and output file
      varp = "psl"                         ; circulation variable name for dynamico input and output file
      tStart = 192001                      ; Initial date for the Observations (both circulation and to be adjusted variables)
      tEnd = 200012                        ; Final date for the Observations
      period = "185001-201812"             ; dataset full period, only used for adjusted variable file name
      period_p = "183601-201512"           ; dataset full period, only used for circulation variable file name
;--- Detrending adjusted variable: the detrending is performed separately for each month 
;--- It is also performed over the period given by tStart and tEnd
;--- if "none" no detrend; "linear", "quad", "cubic" remove trend of given order, 
;--- "splin"  remove cubic spline; "filty": high-pass filter; 
;--- "file_XXXX" remove value based on multimodel ensemble mean from file: file_cmip5, file_cmip6, file_cesm (meaning CESM-LE)
;--- "eemd" remove trend using eemd file with nr IMFs included (one or two)
      REG = "linear"                  ; detrending method for the to-be adjusted variable
      nr = 2                          ; keep one or two IMFs to represent the trend
      cutoff_hp = 120.                ; Freq. cutoff in months for filtm
      cutoff_hpy = 10.                ; Freq. cutoff in months for filty (filter each month)
      knots = 4                       ; number of connecting knots for the spline
      nwt = 11                        ; number of weights for the HF filter (must be odd)
;--- Detrending slp: "none", "linear", "quad", "filty", "file_XXXX", "eemd"
      REGP = "none"                   ; detrending method for the circulation variable
;--- Other input parameters
      scale_m = 0.01                  ; scale circulation variable from model file (e.g to get hPa from Pa units)
      dir_eemd = ROOT+"/codes/eemd/"  ; directory where all eemd files are stored
;--- End of Main Input ----------------------------------------------------------------------------------------
;
;--- Auxilary Input: can be modified if needed
      if (DEBUG) then
        print("   --- Dynamico run on region: "+region)
        print("   --- Algo parameters: "+N_a+", "+N_b+" ,"+niter)
        print("   --- Using "+varp+" from "+ANAP+" and "+var+" from "+ANAT)
        print("   --- Detrending "+var+" based on: "+REG+" and "+varp+": "+REGP)
        if (REG .eq. "eemd" .or. REGP .eq. "eemd") then
          print("   --- EEMD detrending with "+nr+" IMFs removed")
        end if
      end if
;--- File names to be defined by the user if needed (if REG or REGP equals "file")
      if (REG .eq. "file_cmip5" .or. REGP .eq. "file_cmip5") then
        FICS = ROOT+"cmip5/tas_psl_1m_hist_rcp85_cmip5EM_188001-210012.nc"
        FICP = ROOT+"cmip5/tas_psl_1m_hist_rcp85_cmip5EM_188001-210012.nc"
      end if
      if (REG .eq. "file_cesm" .or. REGP .eq. "file_cesm") then
        FICS = ROOT+"ncar_lgens/hist/tas_Amon_CESM1_historical_r_allEM_192001-210012.nc"
        FICP = ROOT+"ncar_lgens/hist/psl_Amon_CESM1_historical_r_allEM_192001-210012.nc"
      end if     
;--- Regions to be defined/modified by the user
      if (region .eq. "EU") then
        flip = "yes"
        minlat_psl = 25
        maxlat_psl = 90
        minlon_psl = -50
        maxlon_psl = 60
        minlat_tas = 30
        maxlat_tas = 75
        minlon_tas = -30
        maxlon_tas = 40
      end if
      if (region .eq. "EA") then
        flip = "yes"
        minlat_psl = 20
        maxlat_psl = 90
        minlon_psl = -90
        maxlon_psl = 180
        minlat_tas = 20
        maxlat_tas = 90
        minlon_tas = -50
        maxlon_tas = 180
      end if
      if (region .eq. "NA") then
        flip = "no"
        minlat_psl = 20
        maxlat_psl = 90
        minlon_psl = 180
        maxlon_psl = 350
        minlat_tas = 20
        maxlat_tas = 90
        minlon_tas = 180
        maxlon_tas = 310
      end if
      if (region .eq. "NH") then
        flip = "no"
        minlat_psl = 20
        maxlat_psl = 90
        minlon_psl = 0
        maxlon_psl = 359.95
        minlat_tas = 20
        maxlat_tas = 90
        minlon_tas = 0
        maxlon_tas = 359.95
      end if   
;--- Load circulation and adjusted variable data files
      fnm_rpsl = DIRO+ANAP+"/"+varp+"_Amon_"+ANAP+"_"+period_p+".nc"
      f_rpsl  = addfile(fnm_rpsl,"r")
      fnm_rtas = DIRO+ANAT+"/"+var+"_Amon_"+ANAT+"_"+period+".nc"
      f_rtas  = addfile(fnm_rtas,"r")
      if (DEBUG) then
        print(" Input "+varp+" file is: "+fnm_rpsl)
        print(" Input "+var+" file is: "+fnm_rtas)
      end if
;--- Random sequence for selecting analogs
      random_setallseed(N_a*N_b*niter, N_a)
;--- Define parameters for Poisson infilling routine.
      guess     = 1                ; use zonal means
      is_cyclic = True             ; if cyclic [global]
      nscan     = 500              ; number of iterations
      eps       = 1.e-2            ; variable dependent
      relc      = 0.6              ; relaxation coefficient
      opti       = 0               ; not used
;--- NO MODIFICATION NEEDED BELOW THIS POINT -----------------------------------------------------------
;------------------------------------------------------------
;--- 1. Loading variables on full period and infill if needed
;------------------------------------------------------------
      timp = f_rpsl->time
      timt = f_rtas->time
      cldr = cd_calendar(timp, -1)
      cldt = cd_calendar(timt, -1)
      irStart_t = ind(cldt .eq. tStart)
      irEnd_t   = ind(cldt .eq. tEnd)
      irStart   = ind(cldr .eq. tStart)
      irEnd     = ind(cldr .eq. tEnd)
      ntime_f   = irEnd - irStart + 1
      ntime_t   = ntime_f - 12
      nyrs = ntime_f/12
      if (DEBUG) then
        print("   Number of months in training period = "+ntime_t)
        print("   Number of months and years in full period = "+ntime_f+" - "+nyrs)
      end if
      delete([/ cldr, cldt/])
;--- Focus on the required period
      timt := timt(irStart_t:irEnd_t)
      timp := timp(irStart:irEnd)
;--- Circulation variable
      if (ANAP .eq. "20CR_V3") then
        var_psl = "psl"
        tmp = f_rpsl->$var_psl$(irStart:irEnd,:,:)
        lat_o = f_rpsl->lat
        lon_o = f_rpsl->lon
        scale = 0.01
      elseif (ANAP .eq. "20CR_V2C") then
        var_psl = "psl"
        tmp = f_rpsl->$var_psl$(irStart:irEnd,:,:)
        lat_o = f_rpsl->lat
        lon_o = f_rpsl->lon
        scale = 0.01
      elseif (ANAP .eq. "HadSLP2") then
        var_psl = "slp"
        tmp = short2flt(f_rpsl->$var_psl$(irStart:irEnd,::-1,:))
        lat_o = f_rpsl->lat(::-1)
        lon_o = f_rpsl->lon
        scale = 1.
      end if
      if (.not.isatt(tmp,"_FillValue")) then
        tmp@_FillValue = tmp@missing_value
      end if
      tmp = (/ tmp * scale /)
      if (DEBUG) then
        print("   Min/Max "+varp+" raw and scaled:")
        printMinMax(tmp, True)
      end if
      if (FILP .eq. "yes") then
        print("   Infilling "+varp)
        poisson_grid_fill(tmp, is_cyclic, guess, nscan, eps, relc, opti)
        print("   Min/Max "+varp+" after Poisson fill")
        printMinMax(tmp, True)
      end if
 ;--- Adjusted variable: here only temperature datasets
      if (ANAT .eq. "GISST") then
        var_tas = "air"
        tmp_t  = f_rtas->$var_tas$(irStart_t:irEnd_t,::-1,:)
        lon_i = f_rtas->lon
        lat_i = f_rtas->lat(::-1)
        tStart_clm = 195101
        tEnd_clm   = 198012
      elseif (ANAT .eq. "BEST") then
        var_tas = "temperature"
        tmp_t  = f_rtas->$var_tas$(irStart_t:irEnd_t,:,:)
        lon_i = f_rtas->lon
        lat_i = f_rtas->lat
        tStart_clm = 195101
        tEnd_clm   = 198012
      elseif (ANAT .eq. "HadCRUT5") then
        var_tas = "tas"
        tmp_t  = dble2flt(f_rtas->$var_tas$(irStart_t:irEnd_t,:,:))
        lon_i = f_rtas->lon
        lat_i = f_rtas->lat
        tStart_clm = 196101
        tEnd_clm   = 199012
      elseif (ANAT .eq. "COWT") then
        var_tas = "temperature_anomaly"
        tmp_t  = f_rtas->$var_tas$(irStart_t:irEnd_t,:,:)
        lon_i = f_rtas->lon
        lat_i = f_rtas->lat
        tStart_clm = 196101
        tEnd_clm   = 199012
      elseif (ANAT .eq. "NOAA") then
        var_tas = "air"
        tmp_t  = f_rtas->$var_tas$(irStart_t:irEnd_t,:,:)
        lon_i = f_rtas->lon
        lat_i = f_rtas->lat
        tStart_clm = 196101
        tEnd_clm   = 199012
      else
        print("   Unknown dataset, exit now !!!") 
        exit      
      end if
      if (.not.isatt(tmp_t,"_FillValue")) then
        tmp_t@_FillValue = tmp_t@missing_value
      end if
      if (DEBUG) then
        print("   Min/Max "+var+" raw values:")
        printMinMax(tmp_t, True)
      end if
      if (FILL .eq. "yes") then
        print("   Infilling "+var)
        poisson_grid_fill(tmp_t, is_cyclic, guess, nscan, eps, relc, opti)
        print("   Min/Max "+var+" after Poisson fill")
        printMinMax(tmp_t, True)
      end if
;--------------------------------------------------
;--- 2. Detrend slp data
;--------------------------------------------------
      if (REGP .eq. "none") then
        print("   No "+varp+" detrending is performed")
      elseif (REGP .eq. "linear") then
        print("   linear detrending of "+varp+" is performed")
        tmp_dtrd = tmp
        do imo = 0,11
          tmp_dtrd(imo::12,:,:) = (/ dtrend_msg_n(ispan(1,nyrs,1),tmp(imo::12,:,:), False, False, 0)  /)
        end do
        tmp := tmp_dtrd
        delete(tmp_dtrd)
      elseif (REGP .eq. "quad") then
        print("   quadratic detrending of "+varp+" is performed")
        tmp_dtrd = tmp
        do imo = 0,11
          tmp_dtrd(imo::12,:,:) = (/ dtrend_quadratic_msg_n(tmp(imo::12,:,:), False, False, 0)  /)
        end do
        tmp := tmp_dtrd
        delete(tmp_dtrd)
      elseif (REGP .eq. "filty") then
        print("   high-pass filtering of "+varp+" is performed")
        opt   = 0
        fca = 1./cutoff_hpy
        tmp_tbf = tmp
        tmp_tbf@long_name = "BW Highpass: "+cutoff_hpy+" month with "+nwt+" weights"
        ihp = 1       
        nsigma = 1.
        wgt_hp = filwgts_lanczos (nwt, ihp, fca, -999., nsigma)
        do imo = 0,11
          tmp_tbf(imo::12,:,:) = (/ wgt_runave_n_Wrap (tmp(imo::12,:,:),wgt_hp,opt,0) /)
        end do
        tmp := tmp_tbf
        delete([/ tmp_tbf, opt, wgt_hp /])
      elseif (isStrSubset(REGP,"file_")) then
        fic_p = addfile(FICP,"r")
        lon_m = fic_p->lon
        lat_m = fic_p->lat
        tim_m = fic_p->time
        gwt  = NormCosWgtGlobe(lat_m)
        cldm = cd_calendar(tim_m, -1)
        ieStart_m = ind(cldm .eq. tStart)
        ieEnd_m   = ind(cldm .eq. tEnd)
        tmp_p = fic_p->psl
        tmp_p = (/ tmp_p * scale_m /) 
        pslm_clm     = clmMonTLL(tmp_p)
        psl_a = calcMonAnomTLL(tmp_p,pslm_clm)
        pslmi= area_hi2lores_Wrap (lon_m,lat_m,psl_a, True, gwt, lon_o,lat_o , False)
        if (DEBUG) then
          print("   forced model "+varp+" response Min/Max before interpolation:")
          printMinMax(psl_a, True)
          print("   forced model "+varp+" response Min/Max after interpolation:")
          printMinMax(pslmi, True)
        end if
;--- subtract forced psl from Observations
        tmp = (/ tmp - pslmi(ieStart_m:ieEnd_m,:,:) /)
        delete([/ pslmi, tmp_p, pslm_clm, psl_a, cldm, gwt, lat_m, lon_m, tim_m, fic_p /])
      else 
        print("   Unknown option to detrend "+varp+" : Exit NOW !")
        exit
      end if
      if (DEBUG) then
        print("   Min/Max "+varp+" after detrending")
        printMinMax(tmp, True)
      end if
;--------------------------------------------------
;--- 3. Detrend adjusted variable data
;--------------------------------------------------
      if (REG .eq. "none") then
        print("   No "+var+" detrending is performed")
      elseif (REG .eq. "linear") then
        print("   linear detrending of "+var+" is performed")
        tmp_dtrd = tmp_t
        do imo = 0,11
          tmp_dtrd(imo::12,:,:) = (/ dtrend_msg_n(ispan(1,nyrs,1),tmp_t(imo::12,:,:), False, False, 0)  /)
        end do
        tmp_t := tmp_dtrd
        delete(tmp_dtrd)
      elseif (REG .eq. "quad") then
        print("   quadratic detrending of "+var+" is performed")
        tmp_dtrd = tmp_t
        do imo = 0,11
          tmp_dtrd(imo::12,:,:) = (/ dtrend_quadratic_msg_n(tmp_t(imo::12,:,:), False, False, 0)  /)
        end do
        tmp_t := tmp_dtrd
        delete(tmp_dtrd)
      elseif (REG .eq. "cubic") then
        dt = ispan(1,nyrs,1)
        do imo = 0,11
          coef = lspoly_n(dt, tmp_t(lat|:,lon|:,time|imo::12), 1, ndeg, 2)
          poly3 = tmp_t
          do it = 0, nyrs -1
            poly3(it,:,:) = (/coef(:,:,0) + coef(:,:,1) * dt(it) + coef(:,:,2) * dt(it)^2 + coef(:,:,3) * dt(it)^3 /)
          end do 
          cubic = tmp_t(imo::12,:,:)
          cubic = (/ tmp_t(imo::12,:,:) - poly3 /)
          tmp_t(imo::12,:,:) = (/ cubic /)
          delete ([/ cubic, poly3, coef /])
        end do
      elseif (REG .eq. "splin") then
        print("   spline detrending of "+var+" is performed")
        xi = fspan(1, nyrs,nyrs)
        yo = tmp_t(lat|:,lon|:,time|0::12)
        do imo = 0,11
          yz = tmp_t(lat|:,lon|:,time|imo::12)
          yo = (/ csa1s(xi,yz,knots,xi) /)
          tmp_clm = dim_avg_n_Wrap(tmp_t(imo::12,:,:), 0)
          tmp_t(imo::12,:,:) = (/ tmp_t(imo::12,:,:) - yo(time|:,lat|:,lon|:) /)
          tmp_t(imo::12,:,:) = (/ tmp_t(imo::12,:,:) + conform_dims(dimsizes(tmp_t(imo::12,:,:)), tmp_clm,(/1,2/) ) /)
          delete([/ yz, tmp_clm /])
        end do
        delete([/ yo, xi /])
      elseif (REG .eq. "filty") then
        print("   high-pass filtering of "+var+" is performed")
        opt   = 0
        fca = 1./cutoff_hpy
        tmp_tbf = tmp_t
        tmp_tbf@long_name = "BW Highpass: "+cutoff_hpy+" month"
        nwt = 33  
        ihp = 1       
        nsigma = 1.
        wgt_hp = filwgts_lanczos (nwt, ihp, fca, -999., nsigma)
        do imo = 0,11
          tmp_tbf(imo::12,:,:) = (/ wgt_runave_n_Wrap (tmp_t(imo::12,:,:),wgt_hp,opt,0) /)
        end do
        tmp_t := tmp_tbf
        delete([/ tmp_tbf, opt, wgt_hp /])
      elseif (REG .eq. "eemd") then          ; eemd trend is X(t) - X(0)
        print("   eemd detrending of "+var+" is performed")
        trdz = tmp_t
        fice = dir_eemd+var+"_"+ANAT+"_monthly_eemd_"+period+".nc"
        fic_e = addfile(fice,"r")
        zbuf = fic_e->tas_imf
        timd = fic_e->time
        clde = cd_calendar(timd,-1)
        ieStart_m = ind(clde .eq. tStart)
        ieEnd_m   = ind(clde .eq. tEnd)
        tas_eemd = zbuf(time|ieStart_m:ieEnd_m,lat|:,lon|:,imf|:)
        dimss = dimsizes(tas_eemd)
        nb_imf = dimss(3)
        print("   Number of IMFs is : "+nb_imf)
        if (nr .eq. 1) then 
          trdz(irStart_t:irEnd_t,:,:) = (/ tas_eemd(:,:,:,nb_imf-1) /)
        end if
        if (nr .eq. 2) then  
          trdz(irStart_t:irEnd_t,:,:) = (/ dim_sum_n_Wrap(tas_eemd(:,:,:,nb_imf-2:),3) /)
        end if
        do imo = 0,11
          trdzi = trdz(imo,:,:)
          trdz(imo::12,:,:) = (/ trdz(imo::12,:,:) - conform_dims(dimsizes(trdz(imo::12,:,:)), trdzi, (/1,2/) ) /)
        end do
        if (DEBUG) then
          print("   EEMD "+var+" Trend Min/Max :")
          printMinMax(trdz, True)
        end if
        tmp_t = (/ tmp_t - trdz /)
        delete([/ trdz, tas_eemd, trdzi, zbuf, dimss /])
      elseif (isStrSubset(REG,"file_")) then
        fic_p = addfile(FICS,"r")
        lon_m = fic_p->lon
        lat_m = fic_p->lat
        tim_m = fic_p->time
        gwt  = NormCosWgtGlobe(lat_m)
        cldm = cd_calendar(tim_m, -1)
        ieStart_m = ind(cldm .eq. tStart)
        ieEnd_m   = ind(cldm .eq. tEnd)
        d_clm = ind(cldm .eq. tStart_clm)
        e_clm = ind(cldm .eq. tEnd_clm)
        tas_p = fic_p->tas
        tasm_clm     = clmMonTLL(tas_p(d_clm:e_clm,:,:))
        tas_a = calcMonAnomTLL(tas_p,tasm_clm)
        tasmi= area_hi2lores_Wrap (lon_m,lat_m,tas_a, True, gwt, lon_i,lat_i , False)
        if (DEBUG) then
          print("   forced model "+var+" response Min/Max before interpolation:")
          printMinMax(tas_a, True)
          print("   forced model "+var+" response Min/Max after interpolation:")
          printMinMax(tasmi, True)
        end if
;--- subtract forced adjusted variable from Observations
        tmp_t = (/ tmp_t - tasmi(ieStart_m:ieEnd_m,:,:) /)
        delete([/ tasmi, tas_p, tasm_clm, tas_a, cldm, gwt, lat_m, lon_m, tim_m, fic_p /])
      else 
        print("   Unknown option to detrend "+var+" : Exit NOW !")
        exit
      end if 
      if (DEBUG) then
        print("   Min/Max "+var+" after detrending")
        printMinMax(tmp_t, True)
      end if
;--- flip tmp and tmp_t if needed
      if (flip .eq. "yes") then
        tmp = lonFlip(tmp)
        tmp_t = lonFlip(tmp_t)
      end if
;--------------------------------------------------------
;--- 4. Pre-processing before time loop and analog search
;--------------------------------------------------------
;--- prepare array for analog algorithm: select regions for both variables
      rpsl_all = tmp(time|:, {lat|minlat_psl:maxlat_psl},{lon|minlon_psl:maxlon_psl})
      rtas_all = tmp_t(time|:, {lat|minlat_tas:maxlat_tas},{lon|minlon_tas:maxlon_tas})
      rpsl = tmp(time|12:, {lat|minlat_psl:maxlat_psl},{lon|minlon_psl:maxlon_psl})
      rtas = tmp_t(time|12:, {lat|minlat_tas:maxlat_tas},{lon|minlon_tas:maxlon_tas})
;--- Getting lat & lon for reshape and weights
      lat              = rpsl_all&lat
      lon              = rpsl_all&lon
      nlat             = dimsizes(lat)
      nlon             = dimsizes(lon)
      lat_tas              = rtas_all&lat
      lon_tas              = rtas_all&lon
      nlatt             = dimsizes(lat_tas)
      nlont             = dimsizes(lon_tas)
      delete([/ tmp, tmp_t, f_rtas, f_rpsl/]) 
;--- Create array to hold dynamically adjusted and reconstructed circulation variable
      etas_adj=conform_dims( (/niter,ntime_f,nlatt,nlont/),rtas_all, (/1,2,3/) )
      etas_adj!0 = "iteration"
      itera=ispan(1,niter,1)
      itera!0 = "iteration"
      etas_adj&iteration = itera
      etas_adj!1 = "time"
      etas_adj&time = timt
      etas_adj!2 = "lat"
      etas_adj&lat = lat_tas
      etas_adj!3 = "lon"
      etas_adj&lon = lon_tas
;--- 
      epsl_rec=conform_dims( (/niter,ntime_f,nlat,nlon/),rpsl_all, (/1,2,3/) )
      epsl_rec!0 = "iteration"
      itera!0 = "iteration"
      epsl_rec&iteration = itera
      epsl_rec!1 = "time"
      epsl_rec&time = timp
      epsl_rec!2 = "lat"
      epsl_rec&lat = lat
      epsl_rec!3 = "lon"
      epsl_rec&lon = lon
      ntt = nyrs -1
/;-------------------------------------------
;--- 5. The leave one out method for analog search
;---    fpsl is the target array, rpsl the training array
;-------------------------------------------
;/
      print("   Analogue algorithm: time loop begins")
      do ii = 0, ntime_f - 1
        iyear = ii/12
        iyearn = tStart/100 + iyear
        isea = ii - (iyear*12)
        print("   time step = "+ii+" - Year: "+iyearn+" - Month: "+(isea+1) )
;--- get the array to be reconstructed
        fpsl = rpsl_all(time|ii, lat|:, lon|:)
        etas = rtas_all(time|ii, lat|:, lon|: )
;--- get the training sample and deal with boundaries
        if (iyear .eq. 0) then      ; target is in first year
          rpsl = (/ rpsl_all(time|12:, lat|:, lon|:) /)
          rtas = (/ rtas_all(time|12:, lat|:, lon|:) /)
        end if
        if (iyear .eq. ntt) then    ; target is in last year
          nend = ntime_f - 13
          rpsl = (/ rpsl_all(time|:nend, lat|:, lon|:) /)
          rtas = (/ rtas_all(time|:nend, lat|:, lon|:) /)
        end if
        if (iyear .gt. 0 .and. iyear .lt. ntt) then  ; normal month
          imb = iyear*12 - 1
          ima = iyear*12 + 12
          rpsl(:imb,:,:) =  (/ rpsl_all(time|:imb, lat|:, lon|:) /)
          rpsl(imb+1:,:,:) =  (/ rpsl_all(time|ima:, lat|:, lon|:) /)
          rtas(:imb,:,:) =  (/ rtas_all(time|:imb, lat|:, lon|:) /)
          rtas(imb+1:,:,:) =  (/ rtas_all(time|ima:, lat|:, lon|:) /)
        end if
;--- Calculate distances
        dist_m = new((/ntt,1/), float)
        if (DIST .eq. "EUCLIDE") then 
          rpsl_2d = onedtond(ndtooned(rpsl),(/ntime_t,nlat*nlon/))
          fpsl_1d = onedtond(ndtooned(fpsl),(/1,nlat*nlon/)) 
          DISTANCEF::distancef(rpsl_2d(isea::12,:), fpsl_1d, ntt, 1, nlat*nlon, dist_m )
          delete([/rpsl_2d,fpsl_1d/])
        end if
        if (DIST .eq. "TEWELES") then
          fpsl_3d = onedtond(ndtooned(fpsl),(/1,nlat,nlon/))
          TEWELESF::tewelesf(rpsl(isea::12,:,:), fpsl_3d, lat, lon, ntt, 1, nlat,nlon,dist_m )
        end if
;--- rank Distances: sort dist increasing order 
        i_dst = dim_pqsort_n(dist_m,2,0)
        delete(dist_m)
        i_dst = (/ i_dst*12 + isea /); reconstruct orig. address
;--- Select the N_a closest analogues
        dst_1d = ndtooned(i_dst(0:N_a - 1,:))
        dpsl = (/ rpsl(dst_1d,:,:) /)
        dtas = (/ rtas(dst_1d,:,:) /)
        delete([/dst_1d,i_dst/])
/;-----------------------------------------------------------
 Constructed analog step using Moore-Penrose pseudo-inverse
 and perform loop on iteration of analog steps
;/
        fpsl_m = onedtond(ndtooned(fpsl),(/1,nlat*nlon/))
        rec_psl_2d = new((/niter,nlat,nlon/), float, rpsl@_FillValue) 
        rec_tas_2d = new((/niter,nlatt,nlont/), float, rtas@_FillValue)
        do ia = 0,niter-1         ; Start of iteration loop
;--- Select randomly N_b analogs among the N_a ones
          nba = generate_unique_indices(N_a)
          iadr = nba(0:N_b-1)
          A_mat = onedtond(ndtooned(dpsl(iadr,:,:)),(/N_b,nlat*nlon/))
          dimss = dimsizes(A_mat)
          nca = dimss(0)
          nra = dimss(1)
          s = new((/nca,nca/),float)
          s = (/ 0. /)
          u   = new ( (/nca,nra/) , typeof(A_mat) )  ; pre-allocate space 
          v   = new ( (/nca,nca/) , typeof(A_mat) )
          sgv  = svd_lapack (A_mat, "S" , "S", 1,  u, v)  ; if 0 v is here v_T
          delete([/ A_mat, dimss /])
          do im = 0, nca - 1
            s(im,im) = (/ 1./sgv(im) /)
          end do
;--- Perform pseudo-inverse to get linear coefficients and apply to analogs
          x = fpsl_m(0,:) # transpose(u) # s # v       ; x has dimension N_b
          delete([/ s, u, v, sgv, nba/])
          rec_tas = x # onedtond(ndtooned(dtas(iadr,:,:)),(/N_b,nlatt*nlont/))
          rec_psl = x # onedtond(ndtooned(dpsl(iadr,:,:)),(/N_b,nlat*nlon/))
          delete([/ x, iadr/])
          rec_tas_2d(ia,:,:) = onedtond(rec_tas,(/nlatt,nlont/))
          rec_psl_2d(ia,:,:) = onedtond(rec_psl,(/nlat,nlon/))
          delete([/ rec_tas, rec_psl/])
        end do               ; end of loop on number of iterations
        etas_adj(:,ii,:,:) = (/rec_tas_2d /)
        epsl_rec(:,ii,:,:) = (/rec_psl_2d /)
        delete([/ rec_tas_2d, rec_psl_2d, fpsl_m, fpsl/])
      end do                 ; end of time loop
;--- 
      delete([/dpsl, dtas /])
      print("   CPU time: " + get_cpu_time())
      print("   Analog done ")
      etas_adj_avg = dim_avg_n_Wrap(etas_adj,0)
      epsl_rec_avg = dim_avg_n_Wrap(epsl_rec,0)
      print("   Write Output ")
;--- Add attributes to output arrays
      if (var .eq. "tas") then
        etas_adj_avg@standard_name = "Average of the dynamical component of air temperature "
        etas_adj_avg@units = "K" 
        etas_adj@standard_name = "Dynamical component of air temperature for all iterations"
        etas_adj@units = "K"  
      end if    
      if (varp .eq. "psl") then
        epsl_rec_avg@standard_name = "Average of the constructed sea level pressure analogues"
        epsl_rec_avg@units = "hPa" 
        epsl_rec@standard_name = "Constructed sea level pressure analogues for all iterations"
        epsl_rec@units = "hPa" 
      end if
/;
-------------------------
---- 6.  Write results
-------------------------
;/
      if (REG .eq. "eemd" .and. nr .eq. 1) then
        REGNM = REG+"_ST"
      elseif (REG .eq. "eemd" .and. nr .eq. 2) then
        REGNM = REG+"_ST+MDV"
      else
        REGNM = REG
      end if
      setfileoption("nc", "Format",  FMTnm)
      filo = var+"_"+DIST+"_"+N_a+"_"+N_b+"_"+niter+"_dynadj_"+ANAP+"_"+REGP+"_"+ANAT+"_"+REGNM+"_"+region+"_"+tStart+"-"+tEnd+".nc"
      filp = varp+"_"+DIST+"_"+N_a+"_"+N_b+"_"+niter+"_recons_"+ANAP+"_"+REGP+"_"+ANAT+"_"+REGNM+"_"+region+"_"+tStart+"-"+tEnd+".nc"
      system("/bin/rm -f " + DIRR + filo)
      system("/bin/rm -f " + DIRR + filp)
      fout  = addfile (DIRR + filo, "c")
      foup  = addfile (DIRR + filp, "c")
      setfileoption(fout,"DefineMode",True)
      setfileoption(foup,"DefineMode",True)
      dimNames = (/"iteration", "lat", "lon", "time"/)
;--- Adjusted variable first  
      dimSizes = (/ niter ,  nlatt,  nlont, ntime_f/) 
      dimUnlim = (/ False, False, False, True/)   
      filedimdef(fout,dimNames,dimSizes,dimUnlim)
      filevardef(fout, "iteration"  ,typeof(itera),getvardims(itera))
      filevardef(fout, "lat"  ,typeof(lat_tas),getvardims(lat_tas))
      filevardef(fout, "lon"  ,typeof(lon_tas),getvardims(lon_tas))
      filevardef(fout, "time"  ,typeof(timt),getvardims(timt))
      filevardef(fout, "tas_a"    ,typeof(etas_adj)  ,getvardims(etas_adj))
      filevardef(fout, "tas_avg"    ,typeof(etas_adj_avg)  ,getvardims(etas_adj_avg))
      filevarattdef(fout,"lat"  ,lat_tas)
      filevarattdef(fout,"lon"  ,lon_tas)
      filevarattdef(fout,"time"  ,timt)
      filevarattdef(fout,"tas_a"  ,etas_adj)
      filevarattdef(fout,"tas_avg"  ,etas_adj_avg)
      setfileoption(fout,"DefineMode",False)
      fout->iteration    = (/itera/)
      fout->lat    = (/lat_tas/)
      fout->lon    = (/lon_tas/)
      fout->time   = (/timt/)
      fout->tas_a     = (/etas_adj/)
      fout->tas_avg     = (/etas_adj_avg/)
      delete(fout)
;--- Circulation variable now
      dimSizes = (/ niter ,  nlat,  nlon, ntime_f/)
      dimUnlim = (/ False, False, False, True/)   
      filedimdef(foup,dimNames,dimSizes,dimUnlim)
      filevardef(foup, "iteration"  ,typeof(itera),getvardims(itera))
      filevardef(foup, "lat"  ,typeof(lat),getvardims(lat))
      filevardef(foup, "lon"  ,typeof(lon),getvardims(lon))
      filevardef(foup, "time"  ,typeof(timp),getvardims(timp))
      filevardef(foup, "psl_a"    ,typeof(epsl_rec)  ,getvardims(epsl_rec))
      filevardef(foup, "psl_avg"    ,typeof(epsl_rec_avg)  ,getvardims(epsl_rec_avg))
      filevarattdef(foup,"lat"  ,lat)
      filevarattdef(foup,"lon"  ,lon)
      filevarattdef(foup,"time"  ,timp)
      filevarattdef(foup,"psl_a"  ,epsl_rec)
      filevarattdef(foup,"psl_avg"  ,epsl_rec_avg)
      setfileoption(foup,"DefineMode",False)
      foup->iteration    = (/itera/)
      foup->lat    = (/lat/)
      foup->lon    = (/lon/)
      foup->time   = (/timp/)
      foup->psl_a   = (/epsl_rec/)
      foup->psl_avg = (/epsl_rec_avg/)
      delete(foup)
;**************************************************************************
end
