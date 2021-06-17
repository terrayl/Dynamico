/;
 --- Program dynamico_daily.ncl
 ---
 --- Code repository: github.com/terrayl/Dynamico
 --- Reference: Terray, L.,2021: A dynamical adjustment perspective on extreme event attribution.
 ---                             Weather and Climate Dynamics, submitted, June 2021
 ---
 --- This code performs dynamical adjustment as described in Deser et al. 2016 J.Climate, but does it using daily data.
 --- The code first reads the observation data for an atmospheric circulation variable, for example daily sea level pressure (slp), and 
 --- another surface variable, named the adjusted variable, for example surface air temperature (sat) or precipitation (pr). 
 --- It also defines two different geographical domains for the circulation and adjusted variables.
 ---
 --- It searches N_a closest analogues of the circulation variable daily means in the same period (a N-day period around the given day)
 --- in all the other years, thereby excluding the current year in the analog search. 
 ---
 --- It then draws N_b analogues among the N_a ones and uses the moore-penrose
 --- pseudo-inverse to get the best linear combination of the N_b analogues that fits the historical
 --- circulation variables daily data. The same coefficients are then applied to the adjusted variable data.
 ---
 --- In case external forcing is present, the adjusted variable might have to be detrended beforehand and several parametric methods 
 --- are available (linear, quadratic, cubic, loess).
 --- The algo is done niter times for each day of the period under scrutiny, then the mean over the niter iterations is performed.
 --- The means are then saved as the best reconstructed circulation analogue and dynamical contribution for the adjusted variable
 --- of interest. and all iterations are also saved to handle uncertainties.
 --- 
 --- ASSUMPTIONS: 1. All data files must have the same order for longitudes (0->360)
 ---                 as the code assumes it is the case (use cdo or nco tricks to do it beforehand if needed)
 ---                 For latitude, the code uses ncl ::-1 trick to have the latitudes as (-90->90)
 ---              2. Assume that all files have correct time coordinates (can be read by cd_calendar)
 ---              3. The infilling is based on the poisson_grid_fill routine with default parameters (check the ncl page for details)
 ---              4. The circulation variable array need to be extrapolated if missing values are present (use FILP = "yes")
 ---              5. Data filenames are read and written out with specific names reflecting our own data structure and organization.
 ---                 It is very likely that the user will have to change the file and directory names.
 ---                 (Look at Input and Auxilary Input of the code to see or change required names)
 ---              6. Users wanting to use the loess filter need to install the NCSTAT package
 ---                 https://terray.locean-ipsl.upmc.fr/ncstat/index.html 
 ---                 NCSTAT loess package is used to detrend (each month is detrended separately)
 ---                 Input parameter nl gives the filter length (in years) that will be multiplied by the number of days per month
 ---                 to get the low-pass filter timescale in days
 ---                 Additional smoothing can also be perfomed (see parameter smt in NCSTAT procedure)
 ---
 --- Author: Laurent Terray 5/04/2020 CECI, CNRS/Cerfacs
 --- The code has been adapted from the monthly mean code. Main difference is that most input parameters are given as arguments of the code, thereby allowing scripting.
 ---
 --- Calling line: ncl -Q var=\"tasmax\" varp=\"psl\" ANAP=\"NCEP\" ANAT=\"NCEP\" daIni=20100715 daFin=20100814 DIST=\"TEWELES\" N_d=15 N_a=400 N_b=200 niter=100 nl=45 region=\"WR\" FILL=\"no\" FILP=\"no\" REG=\"loess\" REGP=\"none\" DEBUG=True FMTnm=\"NetCDF4Classic\" dynamico_daily.ncl
 --- Main Input parameters
 --- Calendar:
 ---      daIni = 20100715  ; initial date to estimate dynamical adjustment
 ---      daFin = 20100814  ; final date to estimate dynamical adjustment
 --- Algorithm:
 ---      DIST = "TEWELES"  ; distance used to get analogues EUCLIDE or TEWELES
 ---      N_d  = 15         ; range for analog days: +/-  N_d days
 ---      N_a  = 400        ; Number of closest analogs
 ---      N_b  = 200        ; Select N_b analogs among N_a
 ---      niter = 100       ; Number of algorithm iterations
 ---      nl = 45           ; Loess filter length (in years)
 ---      region = "WR"     ; geographical domain
 ---      FILL = "no"       ; Infilling or not sat datasets
 ---      FILP = "no"       ; Infilling or not slp datasets
 ---      REG = "loess"     ; detrending method for adjusted variable
 ---      REGP = "none"     ; detrending method for circulation variable
 ---      DEBUG = True      ; Printing checks
 --- Datasets:
 ---      var ="tasmax"        ; adjusted variable name for dynamico input and output file
 ---      varp = "psl"      ; circulation variable name for dynamico input and output file
 ---      ANAP = "NCEP"     ; Observed circulation variable dataset
 ---      ANAT = "NCEP"     ; Observed adjusted variable dataset
 --- Format output file:
 ---      FMTnm = "NetCDF4Classic"
-------------------------------------------------------------------------------------------------------------
;/
;--- Modules
external TEWELESF "~/codes/dynamical_adjustment/tewelesf.so"
external DISTANCEF "~/codes/dynamical_adjustment/distancef.so"
;---
begin
;--- Auxilary input:
;--------------------------------------------------
;--- 1.  Input/output data and directories (need to be adapted by the user tp its own data structure)
;--------------------------------------------------
;---
      ROOT = "/data/scratch/globc/terray/" ; root directory
      DOBS = ROOT+"obs/"                   ; directory for observed datasets
      WORK = ROOT+"codes/dynadj/work/"     ; temporary work directory
      DIRO = ROOT+"codes/dynadj/results/"  ; output files directory
      setfileoption("nc", "Format", FMTnm) ; setting up format of output and work files
;--- names for output files
      fict_out = var+"_distance_"+DIST+"_days_"+N_d+"_"+N_a+"_"+N_b+"_"+niter+"_dynadj_"+ANAP+"_"+REGP+"_"+ANAT+"_"+REG+"_"+region+"_"+daIni+"-"+daFin+".nc"
      ficp_out = varp+"_"+var+"_distance_"+DIST+"_days_"+N_d+"_"+N_a+"_"+N_b+"_"+niter+"_recons_"+ANAP+"_"+REGP+"_"+ANAT+"_"+REG+"_"+region+"_"+daIni+"-"+daFin+".nc"
;--- Period for atmospheric circulation variable 
      if (ANAP .eq. "ERAI") then
        pStart = 19790101                ;--- Dataset initial date
        pEnd   = 20181231                ;--- Dataset final date
      elseif (ANAP .eq. "ERA5") then
        pStart = 19790101
        pEnd   = 20191231
      elseif (ANAP .eq. "20CR_V2C") then
        pStart = 18510101 
        pEnd   = 20181231
      elseif (ANAP .eq. "20CR_V3") then
        pStart = 18360101
        pEnd   = 20181231
      elseif (ANAP .eq. "NCEP") then 
        pStart = 19480101
        pEnd   = 20181231
      else
        print(" Dataset "+ANAP+" not included yet; Stop the code!")
        exit
      end if
      period_p = pStart+"_"+pEnd
      scale = 100.               ;--- scaling(division) for circulation variable i.e to go from Pa to hPa 
;--- Adjusted variable now
      if (ANAT .eq. "ERAI") then
        tStart = 19790101             ;--- Dataset initial date
        tEnd   = 20181231             ;--- Dataset final date
      elseif (ANAT .eq. "ERA5") then
        tStart = 19790101
        tEnd   = 20191231
      elseif (ANAT .eq. "NCEP") then
        tStart = 19480101
        tEnd   = 20181231
      elseif (ANAT .eq. "NCEP") then
        tStart = 19790101                      
        tEnd   = 20181231                       
      elseif (ANAT .eq. "BEST") then
        tStart = 19000101                     
        tEnd   = 20181231
      elseif (ANAT .eq. "20CR_V2C") then
        period = "18510101_20181231"
        tStart = 18510101                      
        tEnd   = 20181231 
      elseif (ANAT .eq. "20CR_V3") then
        tStart = 18360101                     
        tEnd   = 20181231   
      elseif (ANAT .eq. "EOBS") then
        tStart = 19500101                      
        tEnd   = 20181231
      elseif (ANAT .eq. "HadGHCND") then
        tStart = 19500101  
        tEnd   = 20141231
      else
        print(" Dataset "+ANAT+" not included yet; Stop !")
        exit
      end if
      period = tStart+"_"+tEnd
;--- Get period under scrutiny
      iStart = max( (/pStart, tStart/) )
      iEnd   = min( (/pEnd, tEnd/) )
      print(" Common period for dynamical adjustment is: "+iStart+" - "+iEnd)
;--------------------------------------------------
;--- 2.  Region and geographical domain definitions
;--------------------------------------------------
;--- Regions to be defined/modified by the user
      if (region .eq. "EU") then
        flip = "yes"
        minlat_p = 25
        maxlat_p = 90
        minlon_p = -60
        maxlon_p = 100
        minlat_t = 25
        maxlat_t = 90
        minlon_t = -60
        maxlon_t = 50
      end if
      if (region .eq. "WE") then
        flip = "yes"
        minlat_p = 25
        maxlat_p = 80
        minlon_p = -30
        maxlon_p = 45
        minlat_t = 30
        maxlat_t = 75
        minlon_t = -20
        maxlon_t = 35
      end if      
      if (region .eq. "EA") then
        flip = "yes"
        minlat_p = 25
        maxlat_p = 90
        minlon_p = -90
        maxlon_p = 150
        minlat_t = 25
        maxlat_t = 90
        minlon_t = -60
        maxlon_t = 100
      end if
      if (region .eq. "WA") then
        flip = "yes"
        minlat_p = 25
        maxlat_p = 90
        minlon_p = -20
        maxlon_p = 80
        minlat_t = 25
        maxlat_t = 90
        minlon_t = -10
        maxlon_t = 70
      end if
      if (region .eq. "WR") then
        flip = "no"
        minlat_p = 25
        maxlat_p = 80
        minlon_p = 10
        maxlon_p = 65
        minlat_t = 25
        maxlat_t = 80
        minlon_t = 15
        maxlon_t = 60
      end if      
      if (region .eq. "EO") then
        flip = "yes"
        minlat_p = 25
        maxlat_p = 90
        minlon_p = -70
        maxlon_p = 120
        minlat_t = 25
        maxlat_t = 90
        minlon_t = -40
        maxlon_t = 70
      end if
      if (region .eq. "NA") then
        flip = "no"
        minlat_p = 20
        maxlat_p = 90
        minlon_p = 180
        maxlon_p = 350
        minlat_t = 20
        maxlat_t = 90
        minlon_t = 180
        maxlon_t = 310
      end if
      if (region .eq. "NH") then
        flip = "no"
        minlat_p = 20
        maxlat_p = 90
        minlon_p = 0
        maxlon_p = 359.95
        minlat_t = 20
        maxlat_t = 90
        minlon_t = 0
        maxlon_t = 359.95
      end if
;---
;--- End of Auxilary Input ----------------------------------------------------------------------------------------
/;
 ---
 --- *********************************** Beginning of the Main code ************************************************
 ---
;/
      if (DEBUG) then
        print("   --- Dynamico run on region: "+region)
        print("   --- Using daily "+varp+" from "+ANAP+" and "+var+" from "+ANAT)
        print("   --- Detrending "+var+" based on: "+REG+" and "+varp+": "+REGP)
      end if
;---   
;--- Load data files
;---
      fnm_rpsl = DOBS+ANAP+"/"+varp+"_1d_"+period_p+"_"+ANAP+".nc"
      f_rpsl  = addfile(fnm_rpsl,"r")
      fnm_rtas = DOBS+ANAT+"/"+var+"_1d_"+period+"_"+ANAT+".nc"
      f_rtas  = addfile(fnm_rtas,"r")
      if (DEBUG) then
        print(" Input "+varp+" file is: "+fnm_rpsl)
        print(" Input "+var+" file is: "+fnm_rtas)
      end if
;--- Random sequence for selecting analogs
      random_setallseed(N_a*N_b*niter, N_a)
;--- Define parameters for NCL Poisson infilling routine.
      guess     = 1                ; use zonal means
      is_cyclic = True             ; if cyclic [global]
      nscan     = 10               ; number of iterations
      eps       = 5.e-2            ; variable dependent
      relc      = 0.6              ; relaxation coefficient
      opti       = 0               ; not used
;--- NO MODIFICATION NEEDED BELOW THIS POINT -----------------------------------------------------------
/;
---------------------------------------------------
 --- 1. Loading variables on full period and infill
---------------------------------------------------
;/
      timp = f_rpsl->time
      timt = f_rtas->time
      cldr = cd_calendar(timp, -2)
      cldt = cd_calendar(timt, -2)
      irStart_t = ind(cldt .eq. iStart)
      irEnd_t   = ind(cldt .eq. iEnd)
      irStart   = ind(cldr .eq. iStart)
      irEnd     = ind(cldr .eq. iEnd)
      timp := timp(irStart:irEnd)
      timt := timt(irStart_t:irEnd_t)
      ntime_f   = irEnd - irStart + 1
      year_ini = iStart/10000
      year_fin = iEnd/10000
      nyrs = year_fin - year_ini + 1
      if (nl .gt. nyrs) then
        print(" Loess filter length "+nl+" greater than number of years "+nyrs+" ; Stop code Now !!!")
        exit
      end if
      ntime_t   = (nyrs-1) * (2*N_d +1)
      if ( ntime_t .lt. N_a) then
        print(" Training test too small for required number of analogs; Stop code Now !!!")
        exit
      end if
      if (DEBUG) then
        print(" Number of days in training period = "+ntime_t)
        print(" Number of days in full period     = "+ntime_f)
        print(" Number of years in full period    = "+nyrs)
      end if
      delete([/ cldr, cldt/])
;---
;--- Reading circulation variable
;---
      if (ANAP .eq. "20CR_V2C" .or. ANAP .eq. "20CR_V3") then
        tmp = f_rpsl->$varp$
        lat_o = f_rpsl->lat
      elseif (ANAP .eq. "ERAI" .or. ANAP .eq. "ERA5" .or. ANAP .eq. "NCEP") then
        tmp = f_rpsl->$varp$(:,::-1,:)
        lat_o = f_rpsl->lat(::-1)
      end if
      time_p = f_rpsl->time
      lon_o = f_rpsl->lon     
      if (.not.isatt(tmp,"_FillValue")) then
        tmp@_FillValue = tmp@missing_value
      end if
      tmp = (/ tmp/scale /)
      if (DEBUG) then
        print(" Raw "+varp+" Min/Max")
        print(" Min/Max "+varp+"(hPa) = "+min(tmp)+" - "+max(tmp) )
      end if
      if (FILP .eq. "yes") then
        poisson_grid_fill(tmp, is_cyclic, guess, nscan, eps, relc, opti)
        print(" Min/Max after Poisson fill")
        print(" Min/Max "+varp+"(hPa) = "+min(tmp)+" - "+max(tmp) )
      end if
;---
;--- Reading adjusted variable
;---
      if (ANAT .eq. "20CR_V2C" .or. ANAT .eq. "20CR_V3" .or. ANAT .eq. "BEST" .or. ANAT .eq. "EOBS") then
        tmp_t  = f_rtas->$var$
        lat_i = f_rtas->lat
      elseif (ANAT .eq. "ERAI" .or. ANAT .eq. "ERA5" .or. ANAT .eq. "HadGHCND" .or. ANAT .eq. "NCEP") then
        tmp_t  = f_rtas->$var$(:,::-1,:)
        lat_i = f_rtas->lat(::-1)
      else
        print(" Unknown dataset for adjusted variable, exit now !!!") 
        exit      
      end if
      time_t = f_rtas->time
      lon_i  = f_rtas->lon      
      if (.not.isatt(tmp_t,"_FillValue")) then
        tmp_t@_FillValue = tmp_t@missing_value
      end if
      if (DEBUG) then
        print(" Raw "+var+" Min/Max")
        print(" Min/Max "+var+"(K) = "+min(tmp_t)+" - "+max(tmp_t) )
      end if
      if (FILL .eq. "yes") then
        poisson_grid_fill(tmp_t, is_cyclic, guess, nscan, eps, relc, opti)
        print(" Min/Max after Poisson fill")
        print(" Min/Max "+var+"(K) = "+min(tmp_t)+" - "+max(tmp_t) )
      end if
/;
--------------------------------------------------
 --- 2. Detrend circulation variable data
--------------------------------------------------
;/
      modays = (/31,28,31,30,31,30,31,31,30,31,30,31 /)       ;--- Nb days for each month
      if (REGP .eq. "none") then
        print(" No detrending is performed for variable "+varp)
      elseif (REGP .eq. "loess") then
        print(" Loess detrending is performed for variable "+varp)
;--- Get specific string for file naming
        datestr = systemfunc("date")
        datestr = str_sub_str(datestr, " ", "_")
        datestr = str_sub_str(datestr, ":", "_")
;--- Detrend each month separately (to account for seasonality of trends
;--- Select month
        TIME   = cd_calendar(time_p, 0)
        month   = floattointeger( TIME(:,1) )
        finp_nm = WORK+varp+"_inp"+datestr+".nc"
        fout_nm = WORK+varp+"_out"+datestr+".nc"
        v_nm = varp+"_"+"residual"
        do im = 0, 11
          ntl = nl * modays(im)
          system("rm -rf " + finp_nm)
          system("rm -rf " + fout_nm)
          id_d = ind(month .eq. (im+1) )
          tmp_mo = tmp(id_d,:,:)
          dim_p = dimsizes(tmp_mo)
          ntim = dim_p(0)
          finp = addfile(finp_nm,"c")
          filedimdef(finp,"time",-1,True)
          finp->$varp$ = tmp_mo
          cmd = "comp_trend_3d -f="+finp_nm+" -v="+varp+" -nt="+ntl+" -smt="+modays(im)+" -a=residual -robust -o="+fout_nm
          system(cmd)
          fout = addfile(fout_nm,"r")
          resi = fout->$v_nm$
          tmp(id_d,:,:) = (/ resi /)
          delete([/id_d, tmp_mo, dim_p, ntim, finp, fout, resi/])
        end do
        system("rm -rf " + finp_nm)
        system("rm -rf " + fout_nm)
        delete([/TIME, month, finp_nm, fout_nm, v_nm, time_p/])
      elseif (REGP .eq. "linear") then
        print(" Linear detrending is performed for variable: "+varp)
;--- Select month
        TIME   = cd_calendar(time_p, 0)
        month   = floattointeger( TIME(:,1) )
        do im = 0,11
          id_d = ind(month .eq. (im+1) )
          tmp_mo = tmp(id_d,:,:)
          dim_p = dimsizes(tmp_mo)
          ntim = dim_p(0)
          tmp(id_d,:,:) = (/ dtrend_msg_n(ispan(1,ntim,1),tmp_mo, False, False, 0)  /)
          delete([/id_d, tmp_mo, dim_p, ntim/])
        end do
        delete([/TIME, month/])
      elseif (REGP .eq. "quad") then
        print(" Quadratic detrending is performed for variable: "+varp)
        TIME   = cd_calendar(time_p, 0)
        month   = floattointeger( TIME(:,1) )
        do im = 0,11
          id_d = ind(month .eq. (im+1) )
          tmp_mo = tmp(id_d,:,:)          
          tmp(id_d,:,:) = (/ dtrend_quadratic_msg_n(tmp_mo, False, False, 0)  /)
          delete([/id_d, tmp_mo/])
        end do
        delete([/TIME, month/])
      elseif (REG .eq. "cubic") then
        print(" Cubic detrending is performed for variable: "+varp)
        TIME   = cd_calendar(time_t, 0)
        month   = floattointeger( TIME(:,1) )
        do im = 0,11
          id_d = ind(month .eq. (im+1) )
          tmp_mo = tmp(id_d,:,:)
          dim_p = dimsizes(tmp_mo)
          ntim = dim_p(0)
          dt = ispan(1, ntim, 1)
          coef = lspoly_n(dt, tmp_mo(lat|:,lon|:,time|:), 1, 4, 2)  ;--- 4 is here the number of coefficients (4 for a 3rd order polynomial)
          poly3 = tmp
          do it = 0, ntim -1
            poly3(it,:,:) = (/coef(:,:,0) + coef(:,:,1) * dt(it) + coef(:,:,2) * dt(it)^2 + coef(:,:,3) * dt(it)^3 /)
          end do 
          tmp(id_d,:,:) = (/ tmp_mo - poly3 /)
          delete ([/ poly3, coef, id_d, tmp_mo, dim_p, ntim /])
        end do
        delete([/TIME, month/])
      else 
        print("Unknown option to detrend variable "+varp+" : Exit NOW !")
        exit
      end if 
      if (DEBUG) then
        print(" Min/max "+varp+" values after detrending")
        print(" Min/Max = "+min(tmp)+" - "+max(tmp))
      end if
/;
--------------------------------------------------
 --- 3. Detrend adjusted variable data
--------------------------------------------------
;/
      if (REG .eq. "none") then
        print(" No detrending is performed for variable "+var)
      elseif (REG .eq. "loess") then
        print(" Loess detrending is performed for variable "+var)
;--- Get specific string for file naming
        datestr = systemfunc("date")
        datestr = str_sub_str(datestr, " ", "_")
        datestr = str_sub_str(datestr, ":", "_")
;--- Detrend each month separately (to account for seasonality of trends
;--- Select month
        TIME   = cd_calendar(time_t, 0)
        month   = floattointeger( TIME(:,1) )
        finp_nm = WORK+var+"_inp"+datestr+".nc"
        fout_nm = WORK+var+"_out"+datestr+".nc"
        v_nm = var+"_"+"residual"
        system("rm -rf " + finp_nm)
        system("rm -rf " + fout_nm)
        do im = 0, 11
          ntl = nl * modays(im)
          id_d = ind(month .eq. (im+1) )
          tmp_mo = tmp_t(id_d,:,:)
;--- Get grand mean for given month
          tmp_gm = dim_avg_n_Wrap(tmp_mo, 0)
          print(" Min/Max "+var+" before loess filter for month:"+(im+1)+" = "+min(tmp_mo)+" - "+max(tmp_mo) ) 
          dim_t = dimsizes(tmp_mo)
          ntim = dim_t(0)
          finp = addfile(finp_nm,"c")
          filedimdef(finp,"time",-1,True)
          finp->$var$ = tmp_mo
          cmd = "comp_trend_3d -f="+finp_nm+" -v="+var+" -nt="+ntl+" -smt="+modays(im)+" -a=residual -robust -o="+fout_nm
          system(cmd)
          fout = addfile(fout_nm,"r")
          resi = fout->$v_nm$
;--- Add back grand mean
          resi = (/ resi + conform_dims(dimsizes(resi), tmp_gm, (/1, 2/) ) /)
          if (DEBUG) then
            print(" Min/Max "+var+" after loess filter for month:"+(im+1)+" = "+min(resi)+" - "+max(resi) )
          end if
          tmp_t(id_d,:,:) = (/ resi /)
          system("rm -rf " + finp_nm)
          system("rm -rf " + fout_nm)
          delete([/id_d, tmp_mo, dim_t, ntim, finp, fout, resi, cmd/])
        end do
        delete([/TIME, month, finp_nm, fout_nm, v_nm, time_t/])
      elseif (REG .eq. "linear") then
        print(" Linear detrending is performed for variable: "+var)
;--- Select month
        TIME   = cd_calendar(time_t, 0)
        month   = floattointeger( TIME(:,1) )
        do im = 0,11
          id_d = ind(month .eq. (im+1) )
          tmp_mo = tmp_t(id_d,:,:)
          dim_t = dimsizes(tmp_mo)
          ntim = dim_t(0)
          tmp_t(id_d,:,:) = (/ dtrend_msg_n(ispan(1,ntim,1),tmp_mo, False, False, 0)  /)
          delete([/id_d, tmp_mo, dim_t, ntim/])
        end do
        delete([/TIME, month/])
      elseif (REG .eq. "quad") then
        print(" Quadratic detrending is performed for variable: "+var)
        TIME   = cd_calendar(time_t, 0)
        month   = floattointeger( TIME(:,1) )        
        do im = 0,11
          id_d = ind(month .eq. (im+1) )
          tmp_mo = tmp_t(id_d,:,:)         
          tmp_dtrd(id_d,:,:) = (/ dtrend_quadratic_msg_n(tmp_mo, False, False, 0)  /)
          delete([/id_d, tmp_mo/])
        end do
        delete([/TIME, month/])
      elseif (REG .eq. "cubic") then
        print(" Cubic detrending is performed for variable: "+var)
        TIME   = cd_calendar(time_t, 0)
        month   = floattointeger( TIME(:,1) )
        do im = 0,11
          id_d = ind(month .eq. (im+1) )
          tmp_mo = tmp_t(id_d,:,:)
          dim_t = dimsizes(tmp_mo)
          ntim = dim_t(0)
          dt = ispan(1, ntim, 1)
          coef = lspoly_n(dt, tmp_mo(lat|:,lon|:,time|:), 1, 4, 2) ;--- 4 is here the number of coefficients (4 for a 3rd order polynomial)
          poly3 = tmp_t
          do it = 0, ntim -1
            poly3(it,:,:) = (/coef(:,:,0) + coef(:,:,1) * dt(it) + coef(:,:,2) * dt(it)^2 + coef(:,:,3) * dt(it)^3 /)
          end do 
          tmp_t(id_d,:,:) = (/ tmp_mo - poly3 /)
          delete ([/ poly3, coef, id_d, tmp_mo, dim_t, ntim /])
        end do
        delete([/TIME, month/])
      else 
        print("Unknown option to detrend variable "+var+" : Exit NOW !")
        exit
      end if 
      if (DEBUG) then
        print(" Min/max "+var+" values after detrending")
        print(" Min/Max = "+min(tmp_t)+" - "+max(tmp_t))
      end if
;--- flip tmp and tmp_t if needed
      if (flip .eq. "yes") then
        tmp = lonFlip(tmp)
        if (ANAT .ne. "EOBS") then
          tmp_t = lonFlip(tmp_t)
        end if
      end if
/;
---------------------------------------------------------
 --- 4. Pre-processing before time loop and analog search
 ---    Prepare arrays for analog algorithm
---------------------------------------------------------
;/
      rpsl_all = tmp(time|irStart:irEnd, {lat|minlat_p:maxlat_p},{lon|minlon_p:maxlon_p})
      rtas_all = tmp_t(time|irStart_t:irEnd_t, {lat|minlat_t:maxlat_t},{lon|minlon_t:maxlon_t})
;--- buffer working arrays
      rpsl = tmp(time|:ntime_t - 1, {lat|minlat_p:maxlat_p},{lon|minlon_p:maxlon_p})
      rtas = tmp_t(time|:ntime_t - 1, {lat|minlat_t:maxlat_t},{lon|minlon_t:maxlon_t})
;--- Getting lat & lon for reshape and weights
      lat              = rpsl_all&lat
      lon              = rpsl_all&lon
      nlat             = dimsizes(lat)
      nlon             = dimsizes(lon)
      lat_t              = rtas_all&lat
      lon_t              = rtas_all&lon
      nlatt             = dimsizes(lat_t)
      nlont             = dimsizes(lon_t)
      delete([/ tmp, tmp_t, f_rtas, f_rpsl/]) 
;--- Create array to hold dynamically adjusted variable and reconstructed circulation variable
;--- First get dates for dynamical adjustment
      clda_p = cd_calendar(timp,-2)
      da_ini = ind(clda_p .eq. daIni)
      da_fin = ind(clda_p .eq. daFin)
      ntime_da = da_fin - da_ini + 1
      clda_t = cd_calendar(timt,-2)
      dat_ini = ind(clda_t .eq. daIni)
      dat_fin = ind(clda_t .eq. daFin)
      delete([/ clda_p, clda_t /])
      if (DEBUG) then
        print(" Dynamical adjustment from "+daIni+" to "+daFin+" - Total number : "+ntime_da+" days")
      end if
;--- Define arrays, first for adjusted variable
      timt_da = timt(dat_ini:dat_fin)
      etas_adj= new( (/niter,ntime_da,nlatt,nlont/), float, rtas_all@_FillValue)
      etas_adj!0 = "iteration"
      itera=ispan(1,niter,1)
      itera!0 = "iteration"
      etas_adj&iteration = itera
      etas_adj!1 = "time"
      etas_adj&time = timt_da
      etas_adj!2 = "lat"
      etas_adj&lat = lat_t
      etas_adj!3 = "lon"
      etas_adj&lon = lon_t
;--- Circulation variable now
      timp_da = timp(da_ini:da_fin)
      epsl_rec= new( (/niter,ntime_da,nlat,nlon/), float, rpsl_all@_FillValue)
      epsl_rec!0 = "iteration"
      itera!0 = "iteration"
      epsl_rec&iteration = itera
      epsl_rec!1 = "time"
      epsl_rec&time = timp_da
      epsl_rec!2 = "lat"
      epsl_rec&lat = lat
      epsl_rec!3 = "lon"
      epsl_rec&lon = lon
      ntt = nyrs -1
;--- Get calendar for circulation variable
      TIME   = cd_calendar(timp, 0)
      year    = floattointeger( TIME(:,0) )
      month   = floattointeger( TIME(:,1) )
      day     = floattointeger( TIME(:,2) ) 
      ddd     = day_of_year(year, month, day)
      year@calendar = timp@calendar
/;
---------------------------------------------------------------------------------
 --- 5. The leave one out method for analog search
 ---    fpsl is the target array, rpsl the training array
 ---    Handle the calendar subtleties: leap years, ending and beginning of years
---------------------------------------------------------------------------------
;/
      do ida = 0 , ntime_da - 1
;--- add index of initial date to get correct days
        ii = ida + da_ini
;--- Test for leap years and get number of days for the current year
        n_leap = days_in_month(year(ii),2)
        n_tot = 337 + n_leap
        d_day = ddd(ii)
        print(" time step = "+ii+" - Year: "+year(ii)+" - Month: "+month(ii)+" - Day: "+day(ii)+" - D-Day = "+d_day )
;--- get the array to be reconstructed
        fpsl = rpsl_all(time|ii, lat|:, lon|:)
        etas = rtas_all(time|ii, lat|:, lon|: )
;--- get the training sample and deal with beginning and end of years
;--- beginning of the first year
        if (d_day .le. N_d .and. year(ii) .eq. year_ini) then
          idx_t = ind( year .gt. year_ini .and. month .eq. 1 .and. day .le. (2* N_d + 1) ) 
;--- beginning of a normal year
        elseif(d_day .le. N_d .and. year(ii) .gt. year_ini) then
          idx_t =  ind( ( year .ne. year(ii) .and. month .eq. 1 .and. day .le. day(ii)+N_d) .or. ( year .ne. year(ii)-1 .and. month .eq. 12 .and. day .ge. (31 - (N_d - day(ii) ) ) ) )
;--- end of the last year
        elseif(day(ii) .gt. 31-N_d .and. month(ii) .eq. 12 .and. year(ii) .eq. year_fin) then
          idx_t = ind( year .lt. year_fin .and. month .eq. 12 .and. day .gt. (31 - (2*N_d + 1)) )
;--- end of a normal year
        elseif(day(ii) .gt. 31-N_d .and. month(ii) .eq. 12 .and. year(ii) .lt. year_fin .and. n_leap .eq. 28) then
          idx_t =  ind( ( year .ne. year(ii) .and. month .eq. 12 .and. day .ge. (day(ii)-N_d)) .or. ( year .ne. year(ii)+1 .and. month .eq. 1 .and. day .le. (N_d - (31 - day(ii)) ) ) )
        elseif(day(ii) .ge. 31-N_d .and. month(ii) .eq. 12 .and. year(ii) .lt. year_fin .and. n_leap .eq. 29) then
          idx_t =  ind( ( year .ne. year(ii) .and. month .eq. 12 .and. day .ge. (day(ii)-N_d)) .or. ( year .ne. year(ii)+1 .and. month .eq. 1 .and. day .le. (N_d - (31 - day(ii)) ) ) )
;--- Normal case
        else
          idx_t = ind( year .ne. year(ii) .and. (ddd .ge. d_day-N_d .and. ddd .le. d_day+N_d) )
        end if
        rpsl = (/ rpsl_all(time|idx_t, lat|:, lon|:) /)
        rtas = (/ rtas_all(time|idx_t, lat|:, lon|:) /)
        delete(idx_t)
;--- Calculate distances
        dist_m = new((/ntime_t,1/), float)
        if (DIST .eq. "EUCLIDE") then 
          rpsl_2d = onedtond(ndtooned(rpsl),(/ntime_t,nlat*nlon/))
          fpsl_1d = onedtond(ndtooned(fpsl),(/1,nlat*nlon/)) 
          DISTANCEF::distancef(rpsl_2d, fpsl_1d, ntime_t, 1, nlat*nlon, dist_m )
          delete([/rpsl_2d,fpsl_1d/])
        end if
        if (DIST .eq. "TEWELES") then
          fpsl_3d = onedtond(ndtooned(fpsl),(/1,nlat,nlon/))
          TEWELESF::tewelesf(rpsl, fpsl_3d, lat, lon, ntime_t, 1, nlat,nlon,dist_m )
        end if
;--- rank Distances: sort dist increasing order 
        i_dst = dim_pqsort_n(dist_m,2,0)
        delete(dist_m)
;--- Select the N_a closest analogues
        dst_1d = ndtooned(i_dst(0:N_a - 1,:))
        dpsl = (/ rpsl(dst_1d,:,:) /)
        dtas = (/ rtas(dst_1d,:,:) /)
        delete([/dst_1d,i_dst/])
/;-----------------------------------------------------------
 --- Constructed analog step using Moore-Penrose pseudo-inverse
 --- and perform loop on iteration of analog steps
;/
        fpsl_m = onedtond(ndtooned(fpsl),(/1,nlat*nlon/))
        rec_psl_2d = new((/niter,nlat,nlon/), float, rpsl@_FillValue) 
        rec_tas_2d = new((/niter,nlatt,nlont/), float, rtas@_FillValue)
        do ia = 0,niter-1         ; Start of iteration loop
;---
;--- Select randomly N_b analogs among the N_a ones
;---
          nba = generate_unique_indices(N_a)
          iadr = nba(0:N_b-1)
          A_mat = onedtond(ndtooned(dpsl(iadr,:,:)),(/N_b,nlat*nlon/))
          dimss = dimsizes(A_mat)
          nca = dimss(0)
          nra = dimss(1)
          s = new((/nca,nca/),float)
          s = (/ 0. /)
          u   = new ( (/nca,nra/) , typeof(A_mat) )  ;--- pre-allocate space 
          v   = new ( (/nca,nca/) , typeof(A_mat) )
          sgv  = svd_lapack (A_mat, "S" , "S", 1,  u, v)  ;--- performing singular value decomposition 
          delete([/ A_mat, dimss /])                      ;--- 1 returns V, 0 returns v_T
          do im = 0, nca - 1
            s(im,im) = (/ 1./sgv(im) /)
          end do
;---
;--- Perform pseudo-inverse to get linear coefficients and apply to analogs
;---
          x = fpsl_m(0,:) # transpose(u) # s # v       ; x has dimension N_b
          delete([/ s, u, v, sgv, nba/])
          rec_tas = x # onedtond(ndtooned(dtas(iadr,:,:)),(/N_b,nlatt*nlont/))
          rec_psl = x # onedtond(ndtooned(dpsl(iadr,:,:)),(/N_b,nlat*nlon/))
          delete([/ x, iadr/])
          rec_tas_2d(ia,:,:) = onedtond(rec_tas,(/nlatt,nlont/))
          rec_psl_2d(ia,:,:) = onedtond(rec_psl,(/nlat,nlon/))
          delete([/ rec_tas, rec_psl/])
        end do               ; end of loop on number of iterations
        etas_adj(:,ida,:,:) = (/rec_tas_2d /)
        epsl_rec(:,ida,:,:) = (/rec_psl_2d /)
        delete([/ rec_tas_2d, rec_psl_2d, fpsl_m, fpsl/])
      end do                 ; end of time loop
;--- 
      delete([/dpsl, dtas /])
      print("CPU time: " + get_cpu_time())
      print(" Analog done ")
      etas_adj_avg = dim_avg_n_Wrap(etas_adj,0)
      epsl_rec_avg = dim_avg_n_Wrap(epsl_rec,0)
      print(" Write Output ")
/;
-------------------------
 --- 6.  Write results
-------------------------
;/
      da_var = var+"_avg"
      ia_var = var+"_a"
      system("/bin/rm -f " + DIRO + fict_out)
      system("/bin/rm -f " + DIRO + ficp_out)
      fout  = addfile (DIRO + fict_out, "c")
      foup  = addfile (DIRO + ficp_out, "c")
      setfileoption(fout,"DefineMode",True)
      setfileoption(foup,"DefineMode",True)
      dimNames = (/"iteration", "lat", "lon", "time"/)
;--- Adjusted variable 
      dimSizes = (/ niter ,  nlatt,  nlont, ntime_da/) 
      dimUnlim = (/ False, False, False, True/)   
      filedimdef(fout,dimNames,dimSizes,dimUnlim)
      filevardef(fout, "iteration"  ,typeof(itera),getvardims(itera))
      filevardef(fout, "lat"  ,typeof(lat_t),getvardims(lat_t))
      filevardef(fout, "lon"  ,typeof(lon_t),getvardims(lon_t))
      filevardef(fout, "time"  ,typeof(timt_da),getvardims(timt_da))
      filevardef(fout, ia_var  ,typeof(etas_adj)  ,getvardims(etas_adj))
      filevardef(fout, da_var  ,typeof(etas_adj_avg)  ,getvardims(etas_adj_avg))
      filevarattdef(fout,"lat"  ,lat_t)
      filevarattdef(fout,"lon"  ,lon_t)
      filevarattdef(fout,"time"  ,timt_da)
      filevarattdef(fout,ia_var ,etas_adj)
      filevarattdef(fout,da_var ,etas_adj_avg)
      setfileoption(fout,"DefineMode",False)
      fout->iteration    = (/itera/)
      fout->lat    = (/lat_t/)
      fout->lon    = (/lon_t/)
      fout->time   = (/timt_da/)
      fout->$ia_var$     = (/etas_adj/)
      fout->$da_var$     = (/etas_adj_avg/)
;--- Circulation variable
      da_var = varp+"_avg"
      ia_var = varp+"_a"
      dimSizes = (/ niter ,  nlat,  nlon, ntime_da/)
      dimUnlim = (/ False, False, False, True/)   
      filedimdef(foup,dimNames,dimSizes,dimUnlim)
      filevardef(foup, "iteration"  ,typeof(itera),getvardims(itera))
      filevardef(foup, "lat"  ,typeof(lat),getvardims(lat))
      filevardef(foup, "lon"  ,typeof(lon),getvardims(lon))
      filevardef(foup, "time"  ,typeof(timp_da),getvardims(timp_da))
      filevardef(foup, ia_var    ,typeof(epsl_rec)  ,getvardims(epsl_rec))
      filevardef(foup, da_var    ,typeof(epsl_rec_avg)  ,getvardims(epsl_rec_avg))
      filevarattdef(foup,"lat"  ,lat)
      filevarattdef(foup,"lon"  ,lon)
      filevarattdef(foup,"time"  ,timp_da)
      filevarattdef(foup,ia_var  ,epsl_rec)
      filevarattdef(foup,da_var  ,epsl_rec_avg)
      setfileoption(foup,"DefineMode",False)
      foup->iteration    = (/itera/)
      foup->lat    = (/lat/)
      foup->lon    = (/lon/)
      foup->time   = (/timp_da/)
      foup->$ia_var$   = (/epsl_rec/)
      foup->$da_var$ = (/epsl_rec_avg/)
      delete([/ fout, foup /])            ;--- Closing files
;**************************************************************************
end
