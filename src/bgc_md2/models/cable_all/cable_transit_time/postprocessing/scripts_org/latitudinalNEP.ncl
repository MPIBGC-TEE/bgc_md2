load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
print("start ncl")
startyr=(/1901,1901,1901,1901,1901,1901/)
endyr=(/2100,2100,2100,2100,2100,2100/)
nyr=endyr-startyr+1
Active=(/True,True,True,True,True,True/)
npatch=10
nland=5656
npool=9
ntime=365
npool=30
FilePath="/datastore/lu02b/CABLE-traceability/"
;SimName=(/"C/S1/","C/S2/","C/S3/","CN/S1/","CN/S2/","CN/S3/"/)
SimName=(/"CN/S1/","CN/S2/","CN/S3/"/)
nsim=dimsizes(SimName)

fgrid=addfile("gridinfo_NCAR_1.9x2.5_landfrac_revised.nc","r")
finfo=addfile("/datastore/lu02b/CABLE-traceability/C/S1/output/restart_ncar_2101_ndep.nc","r")
latitude=finfo->latitude
nlat=dimsizes(fgrid->latitude)
longitude=where(finfo->longitude .lt. 0,finfo->longitude+360,finfo->longitude)
indlat=ind_nearest_coord(latitude,fgrid->latitude,0)
indlon=ind_nearest_coord(longitude,fgrid->longitude,0)
landfrac=new(nland,float)
area=new(nland,float)

do iland=0,nland-1
   landfrac(iland)=fgrid->land_fraction(indlat(iland),indlon(iland))
   area(iland)=fgrid->area(indlat(iland),indlon(iland))*landfrac(iland)
end do

NPPlat=new((/nsim,nlat,nyr(0)/),float)
HRlat=new((/nsim,nlat,nyr(0)/),float)

isDataReady=False
if(isDataReady)then
   print("hmmm")
   NPPlat=asciiread("NPPlat.txt",(/nsim,nlat,nyr(0)/),"float")
   HRlat=asciiread("HRlat.txt",(/nsim,nlat,nyr(0)/),"float")
else
   do isim=0,nsim-1
      print((/SimName(isim)/))
      do iyr=startyr(isim),endyr(isim)
         print((/iyr/))
         fin=addfile(FilePath+SimName(isim)+"/annual/out_ncar_"+iyr+"_annual.nc","r")
         if(iyr .eq. startyr(isim))then
            iveg=fin->iveg
            patchfrac=fin->patchfrac
            patchfrac=where(iveg .eq. 17, patchfrac@_FillValue,patchfrac)
         end if
         NPP=fin->NPP
         HR=fin->HeteroResp
  
         isnan=isnan_ieee(NPP) .or. isnan_ieee(HR)
         NPP=where(isnan,NPP@_FillValue,NPP)
         HR=where(isnan,HR@_FillValue,HR)

         do ilat=0,nlat-1
            indlandlat:=ind(indlat .eq. ilat)
            if(all(.not. ismissing(indlandlat)))then
               NPPlat(isim,ilat,iyr-startyr(isim))=sum(NPP(:,indlandlat)* conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat)) \
                                                  /sum(conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat))
               HRlat(isim,ilat,iyr-startyr(isim)) =sum(HR (:,indlandlat)* conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat)) \
                                                  /sum(conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat))
            end if
         end do
         print(avg(NPPlat(isim,:,iyr-startyr(isim))))
      end do
   end do
   asciiwrite("NPPlat.txt",NPPlat)
   asciiwrite("HRlat.txt",HRlat)
   asciiwrite("NEPlat.txt",NPPlat-HRlat)
end if

wks = gsn_open_wks("pdf","NEP_latitudinal")
gsn_define_colormap(wks,"rainbow+white+gray")    ; Change color map
  res                          = True

  res@tiMainFontHeightF        = 0.025
  res@tiXAxisString            = "NPP-HR (gC/m2)"
  res@tiYAxisString            = "Latitude"
  res@xyLineThicknesses        = (/4,4,4/)
  res@tmLabelAutoStride                = True

  latitude_new=fgrid->latitude

  latitude_new(0:17)=latitude_new@_FillValue
  latitude_new(91:95)=latitude_new@_FillValue
  res@trXMinF                   = -120
  res@trXMaxF                   = 120
  res@tiMainString = "Latitudinal variation of NEP averaged over 200-year simulation"
  res@xyLineColors := (/"black","black","black","black"/)
  res@xyDashPatterns = (/0,2,16/)
;  print(NPPlat(2,72,:))
;  print(HRlat(2,72,:))
  plot_xy=new(6,graphic)
  plot_xy(0)=gsn_csm_xy(wks,(/dim_avg(NPPlat(:,:,0:199))-dim_avg(HRlat(:,:,0:199))/),latitude_new,res)
  gsn_panel(wks,plot_xy,(/2,3/),False)
end
