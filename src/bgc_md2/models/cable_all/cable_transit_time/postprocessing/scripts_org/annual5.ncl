load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
IsCreate = False
startyr=(/1903,1901,1901,1901,1901,1901/)
;endyr=(/1950,1950,1950,1950,1950,1950/)
endyr=(/2011,2011,2011,2011,2100,2011/)
nyr=endyr-startyr+1
nland=5656
SimName=(/"C/S1","C/S2","C/S3","CN/S1","CN/S2","CN/S3"/)
;SimActive=(/False,False,False,True,True,True/)
SimActive=(/True,True,True,True,True,True/)
SimActive=(/False,False,False,False,True,False/)
nsimactive=num(SimActive)
nsim=dimsizes(SimName)
nplantcasa=3
nlittercasa=3
nsoilcasa=3
nday=365
FilePath="/datastore/lu02b/CABLE-traceability/"+SimName
datainfo=asciiread("cnppool1979_steady_single.csv",(/nland,46/),"float")

area=datainfo(:,7)*1.e+9 ;in m2
area@units="m2"
latitude=datainfo(:,5)
latitude@units="degrees north"
longitude=datainfo(:,6)
longitude@units="degrees east"

do isim=0,nsim-1
   if(SimActive(isim))then
      do iyr=startyr(isim),endyr(isim)
         if(iyr .eq. startyr(isim) .and. startyr(isim) .lt. 2000)then
            system("dmget "+FilePath(isim)+"/output/out_ncar_1???_ndep.nc&")
         else
            if(iyr .eq. 2000 .or. iyr .eq. startyr(isim) .and. startyr(isim) .gt. 2000)then
               system("dmput -r -K 0:40959 "+FilePath(isim)+"/output/out_ncar_1???_ndep.nc&")
               system("dmget "+FilePath(isim)+"/output/out_ncar_2???_ndep.nc&")
            end if
         end if   
         print((/FilePath(isim)+"/output/out_ncar_"+iyr+"_ndep.nc"/))
         fin=addfile(FilePath(isim)+"/output/out_ncar_"+iyr+"_ndep.nc","r")
         NPPdata=dim_sum_n_Wrap(fin->NPP,0)
         NPPdata=(/NPPdata*12.e-6*3600*24/)
         NPPdata@units="gC/m2/year"
         GPPdata=dim_sum_n_Wrap(fin->GPP,0)
         GPPdata=(/GPPdata*12.e-6*3600*24/)
         GPPdata@units="gC/m2/year"
         HeteroRespdata=dim_sum_n_Wrap(fin->HeteroResp,0)
         HeteroRespdata=(/HeteroRespdata*12.e-6*3600*24/)
         HeteroRespdata@units="gC/m2/year"
         LAIdata=dim_avg_n_Wrap(fin->LAI(1:,:,:),0)
         fracCallocdata=dim_avg_n_Wrap(fin->fracCalloc(1:,:,:,:),0)
         cplantdata=dim_avg_n_Wrap(fin->Cplant(1:,:,:,:),0)
         clitterdata=dim_avg_n_Wrap(fin->Clitter(1:,:,:,:),0)
         csoildata=dim_avg_n_Wrap(fin->Csoil(1:,:,:,:),0)
         kplantdata=dim_avg_n_Wrap(fin->kplant(1:,:,:,:),0)
         klitterdata=dim_avg_n_Wrap(fin->klitter(1:,:,:,:),0)
         ksoildata=dim_avg_n_Wrap(fin->ksoil(1:,:,:,:),0)
         xktempdata=dim_avg_n_Wrap(fin->xktemp(1:,:,:),0)
         xkwaterdata=dim_avg_n_Wrap(fin->xkwater(1:,:,:),0)
         xkNlimitingdata=dim_avg_n_Wrap(fin->xkNlimiting(1:,:,:),0)
         xkleafcolddata=dim_avg_n_Wrap(fin->xkleafcold(1:,:,:),0)
         xkleafdrydata=dim_avg_n_Wrap(fin->xkleafdry(1:,:,:),0)
         nplantdata=dim_avg_n_Wrap(fin->Nplant(1:,:,:,:),0)
         nlitterdata=dim_avg_n_Wrap(fin->Nlitter(1:,:,:,:),0)
         nsoildata=dim_avg_n_Wrap(fin->Nsoil(1:,:,:,:),0)
         fromLeaftoLdata=dim_avg_n_Wrap(fin->fromLeaftoL(1:,:,:,:),0)
         fromWoodtoLdata=dim_avg_n_Wrap(fin->fromWoodtoL(1:,:,:,:),0)
         fromRoottoLdata=dim_avg_n_Wrap(fin->fromRoottoL(1:,:,:,:),0)
         fromMettoSdata=dim_avg_n_Wrap(fin->fromMettoS(1:,:,:,:),0)
         fromStrtoSdata=dim_avg_n_Wrap(fin->fromStrtoS(1:,:,:,:),0)
         fromCWDtoSdata=dim_avg_n_Wrap(fin->fromCWDtoS(1:,:,:,:),0)
         fromSOMtoSOMdata=dim_avg_n_Wrap(fin->fromSOMtoSOM(1:,:,:,:),0)
         
         if(iyr .eq. startyr(isim))then
            iveg = fin->iveg
            patchfrac = fin->patchfrac
         end if

         if(.not. IsCreate)then
            system("cp -p out_template_annual.nc "+FilePath(isim)+"/annual/out_ncar_"+iyr+"_annual.nc")
            fout = addfile(FilePath(isim)+"/annual/out_ncar_"+iyr+"_annual.nc","w")
            fout@titile = SimName(isim)+" annual output for traceability analysis"
            fout@creation_date = systemfunc ("date")
            fout@source_file   = FilePath(isim)+"/output/out_ncar_"+iyr+"_ndep.nc"
         else
            system("if [ -f "+FilePath(isim)+"/annual/out_ncar_"+iyr+"_annual.nc ];then rm "+FilePath(isim)+"/annual/out_ncar_"+iyr+"_annual.nc;fi")
            fout = addfile(FilePath(isim)+"/annual/out_ncar_"+iyr+"_annual.nc","c")

;===================================================================
    ; explicitly declare file definition mode. Improve efficiency.
    ;===================================================================
        setfileoption(fout,"DefineMode",True)

    ;===================================================================
    ; create global attributes of the file
    ;===================================================================
        fAtt               = True            ; assign file attributes
	fAtt@title         = "C/S1/ annual output for traceability analysis"
	fAtt@source_file   = FilePath(isim)+"/output/out_ncar_"+iyr+"_ndep.nc"
	fAtt@Conventions   = "None"   
	fAtt@creation_date = systemfunc ("date")        
	fileattdef( fout, fAtt )            ; copy file attributes    

    ;===================================================================
    ; predefine the coordinate variables and their dimensionality
    ; Note: to get an UNLIMITED record dimension, we set the dimensionality
    ; to -1 (or the actual size) and set the dimension name to True.
    ;===================================================================
        dimNames = (/"plant_casa","litter_casa","soil_casa","patch","land"/)  
	dimSizes = (/3,3,3,10,5656/) 
	dimUnlim = (/False,False,False,False,False/)   
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

    ;===================================================================
    ; predefine the the dimensionality of the variables to be written out
    ;===================================================================
    ; Here we are using NCL functions to facilitate defining 
    ; each variable's dimension name(s) and type. 
    ; The following could be replaced with explicit, user defined dimension 
    ; names different from those associated with the variable in memory. 
    ; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
    ; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/)) 
    ;===================================================================
       print((/getvardims(iveg)/))
       filevardef(fout, "iveg"             ,typeof(iveg)             ,(/"patch","land"/)) 
       filevardef(fout, "patchfrac"        ,typeof(patchfrac)        ,(/"patch","land"/))
       filevardef(fout, "latitude"         ,typeof(latitude)         ,(/"land"/))
       filevardef(fout, "longitude"        ,typeof(longitude)        ,(/"land"/))
       filevardef(fout, "area"             ,typeof(area)             ,(/"land"/))
       filevardef(fout, "NPP"              ,typeof(NPPdata)          ,(/"patch","land"/)) 
       filevardef(fout, "GPP"              ,typeof(GPPdata)          ,(/"patch","land"/))                           
       filevardef(fout, "HeteroResp"       ,typeof(HeteroRespdata)   ,(/"patch","land"/))                          
       filevardef(fout, "LAI"              ,typeof(LAIdata)          ,(/"patch","land"/))                              
       filevardef(fout, "fracCalloc"       ,typeof(fracCallocdata)   ,(/"plant_casa","patch","land"/))                              
       filevardef(fout, "Cplant"           ,typeof(cplantdata)       ,(/"plant_casa","patch","land"/))                          
       filevardef(fout, "Clitter"          ,typeof(clitterdata)      ,(/"litter_casa","patch","land"/))                  
       filevardef(fout, "Csoil"            ,typeof(csoildata)        ,(/"soil_casa","patch","land"/))                
       filevardef(fout, "kplant"           ,typeof(kplantdata)       ,(/"plant_casa","patch","land"/))                 
       filevardef(fout, "klitter"          ,typeof(klitterdata)      ,(/"litter_casa","patch","land"/))                 
       filevardef(fout, "ksoil"            ,typeof(ksoildata)        ,(/"soil_casa","patch","land"/))                 
       filevardef(fout, "xktemp"           ,typeof(xktempdata)       ,(/"patch","land"/))                          
       filevardef(fout, "xkwater"          ,typeof(xkwaterdata)      ,(/"patch","land"/))                          
       filevardef(fout, "xkNlimiting"      ,typeof(xkNlimitingdata)  ,(/"patch","land"/))                          
       filevardef(fout, "xkleafcold"       ,typeof(xkleafcolddata)   ,(/"patch","land"/))                          
       filevardef(fout, "xkleafdry"        ,typeof(xkleafdrydata)    ,(/"patch","land"/))                          
       filevardef(fout, "Nplant"           ,typeof(nplantdata)       ,(/"plant_casa","patch","land"/))                 
       filevardef(fout, "Nlitter"          ,typeof(nlitterdata)      ,(/"litter_casa","patch","land"/))                 
       filevardef(fout, "Nsoil"            ,typeof(nsoildata)        ,(/"soil_casa","patch","land"/))                 
       filevardef(fout, "fromLeaftoL"      ,typeof(fromLeaftoLdata)  ,(/"litter_casa","patch","land"/))                     
       filevardef(fout, "fromWoodtoL"      ,typeof(fromWoodtoLdata)  ,(/"litter_casa","patch","land"/))                     
       filevardef(fout, "fromRoottoL"      ,typeof(fromRoottoLdata)  ,(/"litter_casa","patch","land"/))                     
       filevardef(fout, "fromMettoS"       ,typeof(fromMettoSdata)   ,(/"soil_casa","patch","land"/))                      
       filevardef(fout, "fromStrtoS"       ,typeof(fromStrtoSdata)   ,(/"soil_casa","patch","land"/))                      
       filevardef(fout, "fromCWDtoS"       ,typeof(fromCWDtoSdata)   ,(/"soil_casa","patch","land"/))                      
       filevardef(fout, "fromSOMtoSOM"     ,typeof(fromSOMtoSOMdata) ,(/"soil_casa","patch","land"/))                        
       
    ;===================================================================
    ; Copy attributes associated with each variable to the file
    ; All attributes associated with each variable will be copied.
    ;====================================================================
       filevarattdef(fout, "iveg"             ,iveg) 
       filevarattdef(fout, "patchfrac"        ,patchfrac)   
       filevarattdef(fout, "latitude"         ,latitude)   
       filevarattdef(fout, "longitude"        ,longitude)
       filevarattdef(fout, "area"             ,area)       
       filevarattdef(fout, "NPP"              ,NPPdata)   
       filevarattdef(fout, "GPP"              ,GPPdata)  
       filevarattdef(fout, "HeteroResp"       ,HeteroRespdata)
       filevarattdef(fout, "LAI"              ,LAIdata)      
       filevarattdef(fout, "fracCalloc"       ,fracCallocdata)      
       filevarattdef(fout, "Cplant"           ,cplantdata)   
       filevarattdef(fout, "Clitter"          ,clitterdata) 
       filevarattdef(fout, "Csoil"            ,csoildata)  
       filevarattdef(fout, "kplant"           ,kplantdata)
       filevarattdef(fout, "klitter"          ,klitterdata)
       filevarattdef(fout, "ksoil"            ,ksoildata) 
       filevarattdef(fout, "xktemp"           ,xktempdata)
       filevarattdef(fout, "xkwater"          ,xkwaterdata)     
       filevarattdef(fout, "xkNlimiting"      ,xkNlimitingdata)     
       filevarattdef(fout, "xkleafcold"       ,xkleafcolddata) 
       filevarattdef(fout, "xkleafdry"        ,xkleafdrydata) 
       filevarattdef(fout, "Nplant"           ,nplantdata)   
       filevarattdef(fout, "Nlitter"          ,nlitterdata) 
       filevarattdef(fout, "Nsoil"            ,nsoildata)  
       filevarattdef(fout, "fromLeaftoL"      ,fromLeaftoLdata) 
       filevarattdef(fout, "fromWoodtoL"      ,fromWoodtoLdata)
       filevarattdef(fout, "fromRoottoL"      ,fromRoottoLdata) 
       filevarattdef(fout, "fromMettoS"       ,fromMettoSdata) 
       filevarattdef(fout, "fromStrtoS"       ,fromStrtoSdata)
       filevarattdef(fout, "fromCWDtoS"       ,fromCWDtoSdata)
       filevarattdef(fout, "fromSOMtoSOM"     ,fromSOMtoSOMdata)
       

    ;===================================================================
    ; explicitly exit file definition mode. **NOT REQUIRED**
    ;===================================================================
        setfileoption(fout,"DefineMode",False)

       end if

    ;===================================================================
    ; output only the data values since the dimensionality and such have
    ; been predefined. The "(/", "/)" syntax tells NCL to only output the
    ; data values to the predefined locations on the file.
    ;====================================================================
       fout->iveg         = (/iveg/)     
       fout->patchfrac    = (/patchfrac/)
       fout->latitude     = (/latitude/)
       fout->longitude    = (/longitude/) 
       fout->area         = (/area/)
       fout->NPP          = (/NPPdata/)
       fout->GPP          = (/GPPdata/)
       fout->HeteroResp   = (/HeteroRespdata/) 
       fout->LAI          = (/LAIdata/)
       fout->fracCalloc   = (/fracCallocdata/)
       fout->Cplant       = (/cplantdata/)
       fout->Clitter      = (/clitterdata/)
       fout->Csoil        = (/csoildata/)
       fout->kplant       = (/kplantdata/)
       fout->klitter      = (/klitterdata/)
       fout->ksoil        = (/ksoildata/)
       fout->xktemp       = (/xktempdata/)
       fout->xkwater      = (/xkwaterdata/)
       fout->xkNlimiting  = (/xkNlimitingdata/)
       fout->xkleafcold   = (/xkleafcolddata/)
       fout->xkleafdry    = (/xkleafdrydata/)
       fout->Nplant       = (/nplantdata/)
       fout->Nlitter      = (/nlitterdata/)
       fout->Nsoil        = (/nsoildata/)
       fout->fromLeaftoL  = (/fromLeaftoLdata/)
       fout->fromWoodtoL  = (/fromWoodtoLdata/)
       fout->fromRoottoL  = (/fromRoottoLdata/)
       fout->fromMettoS   = (/fromMettoSdata/)
       fout->fromStrtoS   = (/fromStrtoSdata/)
       fout->fromCWDtoS   = (/fromCWDtoSdata/)
       fout->fromSOMtoSOM = (/fromSOMtoSOMdata/)

      end do
   end if
end do
end 
