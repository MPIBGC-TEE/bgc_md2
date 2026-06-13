load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
startyr=2008
nyr=2011-startyr+1
endyr=2011+nyr-1

PFTName=(/"ENF","EBF","DNF","DBF","SHRUB","C3 grass","C4 grass","Tundra","Crop","not used","not used","not used","not used","Barren","not used","not used","Ice"/)
;fin=addfile("annual_allnew.nc","r")
;fin=addfile("annual_all_ratio_adjusted.nc","r")
fin=addfile("gridinfo_NCAR_1.9x2.5_landfrac_revised.nc","r")
;Tdata=fin->TBOT
;Pdata=fin->PBOT
;FSdata=fin->FSDS
;FLdata=fin->FLDS
;Qdata=fin->QBOT
;Rdata=fin->RAIN
;SNdata=fin->SNOW
;Wdata=fin->WIND
patchfrac=fin->patchfrac
;Tdata&time=fspan(1901,2100,200)
;Pdata&time=fspan(1901,2100,200)
;FSdata&time=fspan(1901,2100,200)
;FLdata&time=fspan(1901,2100,200)
;Qdata&time=fspan(1901,2100,200)
;Rdata&time=fspan(1901,2100,200)
;SNdata&time=fspan(1901,2100,200)
;Wdata&time=fspan(1901,2100,200)
tmp=patchfrac(0,:,:)
tmp!0="latitude"
tmp!1="longitude"
tmp&latitude:=fin->latitude
tmp&longitude:=fin->longitude

wks = gsn_open_wks("pdf","patchfrac")
gsn_define_colormap(wks,"Rainbow+white+gray")    ; Change color map

      res                      = True

      res@cnFillMode           = "RasterFill"
      res@cnFillOn             = True      ; Turn on contour fill
      res@cnLinesOn            = False     ; Turn off contour lines
      res@lbLabelAutoStride    = True      ; Control labelbar labels
      res@lbBoxLinesOn         = False     ; Turn off lbar box lines
      res@gsnSpreadColors      = True  ; Span full color map
      res@gsnAddCyclic         = True
      res@cnLevelSelectionMode = "ExplicitLevels"
      res@gsnSpreadColorStart  = 48     ; Start at color index 14
      res@gsnSpreadColorEnd    = -3    ; Stop at 3rd color from end

;do iyr=startyr,endyr
;   res@tiMainString = "year="+iyr+" Temperature"
;   res@cnLevels=fspan(240,310,51)
;   plot=gsn_csm_contour_map(wks,Tdata({iyr},:,:),res)     
;end do
;delete(res@cnLevels)

;do iyr=startyr,endyr
;   res@tiMainString = "year="+iyr+" Pressure"
;   res@cnLevels=fspan(60000,104000,51)
 ;  plot=gsn_csm_contour_map(wks,Pdata({iyr},:,:),res)     
;end do

;do iyr=startyr,endyr
;   res@tiMainString = "year="+iyr+" Shortwave radiation"
 ;  res@cnLevels=fspan(80,320,51)
 ;;  plot=gsn_csm_contour_map(wks,FSdata({iyr},:,:),res)     
;end do

;do iyr=startyr,endyr
;   res@tiMainString = "year="+iyr+" Longwave radiation"
;   res@cnLevels=fspan(100,475,51)
;   plot=gsn_csm_contour_map(wks,FLdata({iyr},:,:),res)     
;end do

;do iyr=startyr,endyr
;   res@tiMainString = "year="+iyr+" Sepcific humidity"
;   res@cnLevels=fspan(0,0.02,51)
;   plot=gsn_csm_contour_map(wks,Qdata({iyr},:,:),res)     
;end do

;do iyr=startyr,endyr
;   res@tiMainString = "year="+iyr+" Rainfall"
;   res@cnLevels=fspan(0,0.0001,51)
;   plot=gsn_csm_contour_map(wks,Rdata({iyr},:,:),res)     
;end do

;do iyr=startyr,endyr
;   res@tiMainString = "year="+iyr+" Snowfall"
;   res@cnLevels=fspan(0,1e-5,51)
;   plot=gsn_csm_contour_map(wks,SNdata({iyr},:,:),res)     
;end do

;do iyr=startyr,endyr
;;   res@tiMainString = "year="+iyr+" Wind"
 ;  res@cnLevels=fspan(0,10,51)
;   plot=gsn_csm_contour_map(wks,Wdata({iyr},:,:),res)     
;end do

;res@cnLevels=fspan(-5,5,51)
;res@tiMainString = "year=2012-2011 Temperature"
;tmp=Tdata({2016},:,:)
;tmp=(/dim_avg_n(Tdata({2012:2012},:,:),0)-dim_avg_n(Tdata({2011:2011},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2011-2010 Temperature"
;tmp=(/dim_avg_n(Tdata({2011:2011},:,:),0)-dim_avg_n(Tdata({2010:2010},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2010-2009 Temperature"
;tmp=(/dim_avg_n(Tdata({2010:2010},:,:),0)-dim_avg_n(Tdata({2009:2009},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)

;res@cnLevels=fspan(-50,50,51)
;res@tiMainString = "year=(2013-2012) Shortwave radiation(W)"
;tmp=FSdata({2016},:,:)
;tmp=(/dim_avg_n(FSdata({2013:2013},:,:),0)-dim_avg_n(FSdata({2012:2012},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=(2012-2011) Shortwave radiation(W)"
;tmp=(/dim_avg_n(FSdata({2012:2012},:,:),0)-dim_avg_n(FSdata({2011:2011},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=(2011-2010) Shortwave radiation(W)"
;tmp=(/dim_avg_n(FSdata({2011:2011},:,:),0)-dim_avg_n(FSdata({2010:2010},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=(2010-2009) Shortwave radiation(W)"
;tmp=(/dim_avg_n(FSdata({2010:2010},:,:),0)-dim_avg_n(FSdata({2009:2009},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;
;res@cnLevels=fspan(-1000,1000,51)
;tmp=Rdata({2016},:,:)
;res@tiMainString = "year=2013-2012 Rainfall(mm)"
;tmp=(/dim_avg_n(Rdata({2013:2013},:,:),0)-dim_avg_n(Rdata({2012:2012},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2012-2011 Rainfall(mm)"
;tmp=(/dim_avg_n(Rdata({2012:2012},:,:),0)-dim_avg_n(Rdata({2011:2011},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2011-2010 Rainfall(mm)"
;tmp=(/dim_avg_n(Rdata({2011:2011},:,:),0)-dim_avg_n(Rdata({2010:2010},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2010-2009 Rainfall(mm)"
;tmp=(/dim_avg_n(Rdata({2010:2010},:,:),0)-dim_avg_n(Rdata({2009:2009},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)

;res@cnLevels=fspan(-200,200,51)
;tmp=SNdata({2016},:,:)
;res@tiMainString = "year=2013-2012 Snowfall(mm)"
;tmp=(/dim_avg_n(SNdata({2013:2013},:,:),0)-dim_avg_n(SNdata({2012:2012},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2012-2011 Snowfall(mm)"
;tmp=(/dim_avg_n(SNdata({2012:2012},:,:),0)-dim_avg_n(SNdata({2011:2011},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2011-2010 Snowfall(mm)"
;tmp=(/dim_avg_n(SNdata({2011:2011},:,:),0)-dim_avg_n(SNdata({2010:2010},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2010-2009 Snowfall(mm)"
;tmp=(/dim_avg_n(SNdata({2010:2010},:,:),0)-dim_avg_n(SNdata({2009:2009},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@cnLevels=fspan(-3,5,51)
;res@tiMainString = "year=2020-2010 Wind"
;tmp=Wdata({2016},:,:)
;tmp=(/dim_avg_n(Wdata({2016:2020},:,:),0)-dim_avg_n(Wdata({2006:2010},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2030-2020 Wind"
;tmp=(/dim_avg_n(Wdata({2026:2030},:,:),0)-dim_avg_n(Wdata({2016:2020},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
;res@tiMainString = "year=2010-2000 Wind"
;tmp=(/dim_avg_n(Wdata({2006:2010},:,:),0)-dim_avg_n(Wdata({1996:2000},:,:),0)/)
;plot=gsn_csm_contour_map(wks,tmp,res)
npft=dimsizes(patchfrac(:,0,0))
do ipft=0,npft-1
   res@tiMainString = "Fraction of functional type:"+PFTName(ipft)
   tmp=(/patchfrac(ipft,:,:)/)
   plot=gsn_csm_contour_map(wks,tmp,res)
end do

res@tiMainString = "Index of PFT distribution"
tmp=(/dim_maxind(patchfrac,0)/)
plot=gsn_csm_contour_map(wks,tmp,res)


end
   




