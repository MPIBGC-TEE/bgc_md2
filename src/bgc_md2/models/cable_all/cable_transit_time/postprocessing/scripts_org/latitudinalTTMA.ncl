load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
print("start ncl")
startyr=(/1901,1901,1980,1901,1901,1901/)
endyr=(/2100,2100,2100,2100,2100,2100/)
nyr=endyr-startyr+1
Active=(/True,True,True,True,True,True/)
npatch=10
nland=5656
npool=9
ntime=365
npool=30
nC14=157
FilePath="/datastore/lu02b/CABLE-traceability/"
SimName=(/"C/S1/","C/S2/","C/S3/","CN/S1/","CN/S2/","CN/S3/"/)
nsim=dimsizes(SimName)

fgrid=addfile("gridinfo_NCAR_1.9x2.5_landfrac_revised.nc","r")
finfo=addfile("/datastore/lu02b/CABLE-traceability/C/S1/output/restart_ncar_2101_ndep.nc","r")
latitude=finfo->latitude
nlat=dimsizes(fgrid->latitude)
longitude=where(finfo->longitude .lt. 0,finfo->longitude+360,finfo->longitude)
dataob=asciiread("Latidinal_turnover.csv",(/136,2/),"float")
nob=dimsizes(dataob(:,0))
datamodel=asciiread("Latidinal_turnover_model.csv",(/119,2/),"float")
nmo=dimsizes(datamodel(:,0))
dataob_lb=asciiread("Latidinal_turnover_lb.csv",(/163,2/),"float")
nob_lb=dimsizes(dataob_lb(:,0))
dataob_rb=asciiread("Latidinal_turnover_rb.csv",(/132,2/),"float")
nob_rb=dimsizes(dataob_rb(:,0))
do i=0,nob-2
   if(dataob(i+1,1) .le. dataob(i,1))then
      dataob(i+1,1)=dataob(i,1)+0.00001
   end if
end do
do i=0,nmo-2
   if(datamodel(i+1,1) .le. datamodel(i,1))then
      datamodel(i+1,1)=datamodel(i,1)+0.00001
   end if
end do
do i=0,nob_lb-2
   if(dataob_lb(i+1,1) .le. dataob_lb(i,1))then
      dataob_lb(i+1,1)=dataob_lb(i,1)+0.00001
   end if
end do
do i=0,nob_rb-2
   if(dataob_rb(i+1,1) .le. dataob_rb(i,1))then
      dataob_rb(i+1,1)=dataob_rb(i,1)+0.00001
   end if
end do


print(dataob)

indlat=ind_nearest_coord(latitude,fgrid->latitude,0)
indlon=ind_nearest_coord(longitude,fgrid->longitude,0)
landfrac=new(nland,float)
area=new(nland,float)

;print(latitude)
do iland=0,nland-1
   landfrac(iland)=fgrid->land_fraction(indlat(iland),indlon(iland))
   area(iland)=fgrid->area(indlat(iland),indlon(iland))*landfrac(iland)
end do

Cplantagelat=new((/nsim,nlat,nyr(0)/),float)
Clitteragelat=new((/nsim,nlat,nyr(0)/),float)
Csoilagelat=new((/nsim,nlat,nyr(0)/),float)
ClitterTTlat=new((/nsim,nlat,nyr(0)/),float)
CsoilTTlat=new((/nsim,nlat,nyr(0)/),float)
Cturnoveragelat=new((/nsim,nlat,nyr(0)/),float)

isDataReady=True
if(isDataReady)then
   print("hmmm")
   Cplantagelat=asciiread("Cplantagelat.txt",(/nsim,nlat,nyr(0)/),"float")
   Clitteragelat=asciiread("Clitteragelat.txt",(/nsim,nlat,nyr(0)/),"float")
   Csoilagelat=asciiread("Csoilagelat.txt",(/nsim,nlat,nyr(0)/),"float")
   ClitterTTlat=asciiread("ClitterTTlat.txt",(/nsim,nlat,nyr(0)/),"float")
   CsoilTTlat=asciiread("CsoilTTlat.txt",(/nsim,nlat,nyr(0)/),"float")
   Cturnoveragelat=asciiread("Cturnoveragelat.txt",(/nsim,nlat,nyr(0)/),"float")
else
   do isim=0,nsim-1
      print((/SimName(isim)/))
      do iyr=startyr(isim),endyr(isim)
         print((/iyr/))
         finage=addfile(FilePath+SimName(isim)+"/annual_TT/out_ncar_"+iyr+"_annual_age.nc","r")
         fincpool=addfile(FilePath+SimName(isim)+"/annual/out_ncar_"+iyr+"_annual.nc","r")
         finhr=addfile(FilePath+SimName(isim)+"/annual_TT/out_ncar_"+iyr+"_annualTT.nc","r")
         if(iyr .eq. startyr(isim))then
            iveg=fincpool->iveg
            patchfrac=fincpool->patchfrac
            patchfrac=where(iveg .eq. 17, patchfrac@_FillValue,patchfrac)
         end if
         Cplantage=finage->Cplant
         Clitterage=finage->Clitter
         Csoilage=finage->Csoil
         Cplant=fincpool->Cplant
         Clitter=fincpool->Clitter
         Csoil=fincpool->Csoil
         hr=fincpool->HeteroResp
  
         isnan_age=(dim_num_n(isnan_ieee(Cplantage),0) .gt. 0) .or. (dim_num_n(abs(Cplantage) .gt. 1.e+6,0) .gt. 0) \
              .or. (dim_num_n(isnan_ieee(Clitterage),0) .gt. 0) .or. (dim_num_n(abs(Clitterage) .gt. 1.e+6,0) .gt. 0) \
              .or. (dim_num_n(isnan_ieee(Csoilage),0) .gt. 0) .or. (dim_num_n(abs(Csoilage) .gt. 1.e+6,0) .gt. 0)
         isnan_cpool=(dim_num_n(isnan_ieee(Cplant),0) .gt. 0) .or. (dim_num_n(isnan_ieee(Clitter),0) .gt. 0) .or. (dim_num_n(isnan_ieee(Csoil),0) .gt. 0)
         isnan=conform(Cplantage,isnan_age .or. isnan_cpool,(/1,2/))
         Cplantage=where(isnan,Cplantage@_FillValue,Cplantage)
         Clitterage=where(isnan,Clitterage@_FillValue,Clitterage)
         Csoilage=where(isnan,Csoilage@_FillValue,Csoilage)

         Cplanthr=finhr->Cplant
         Clitterhr=finhr->Clitter
         Csoilhr=finhr->Csoil

         Cplanthr=where(isnan, Cplanthr@_FillValue,Cplanthr)
         Clitterhr=where(isnan, Clitterhr@_FillValue,Clitterhr)
         Csoilhr=where(isnan,Csoilhr@_FillValue,Csoilhr)

         Cplant(1,:,:)=fincpool->Cplant(2,:,:)
         Cplant(2,:,:)=fincpool->Cplant(1,:,:)
         Cplant=where(isnan,Cplant@_FillValue,Cplant)
         Clitter=where(isnan,Clitter@_FillValue,Clitter)
         Csoil=where(isnan,Csoil@_FillValue,Csoil)
         Cplant(2,:,:)=where(iveg .eq. 6 .or. iveg .eq. 7 .or. iveg .eq. 9,0,Cplant(2,:,:))
         Cplant(0,:,:)=where(iveg .eq. 14,0,Cplant(0,:,:))
         Cplant(1,:,:)=where(iveg .eq. 14,0,Cplant(1,:,:))
         Cplant(2,:,:)=where(iveg .eq. 14,0,Cplant(2,:,:))
         Clitter(0,:,:)=where(iveg .eq. 14,0,Clitter(0,:,:))
         Clitter(1,:,:)=where(iveg .eq. 14,0,Clitter(1,:,:))
         Clitter(2,:,:)=where(iveg .eq. 14,0,Clitter(2,:,:))
         Csoil(0,:,:)=where(iveg .eq. 14,0,Csoil(0,:,:))
         Csoil(1,:,:)=where(iveg .eq. 14,0,Csoil(1,:,:))
         Csoil(2,:,:)=where(iveg .eq. 14,0,Csoil(2,:,:))
         Clitter(2,:,:)=where(iveg .eq. 6 .or. iveg .eq. 7 .or. iveg .eq. 9,0,Clitter(2,:,:))

         Cplant=where(Cplant .eq. 0, 0.000001,Cplant)
         Clitter=where(Clitter .eq. 0, 0.000001,Clitter)
         Csoil=where(Csoil .eq. 0,0.000001,Csoil)

         do ilat=0,nlat-1
            indlandlat:=ind(indlat .eq. ilat)
            if(all(.not. ismissing(indlandlat)))then
               ClitterTTlat(isim,ilat,iyr-startyr(isim))=sum(Clitterage(:,:,indlandlat)*Clitterhr(:,:,indlandlat)* \
                          conform(Clitterhr(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/))) \
                         /sum((Cplanthr(:,:,indlandlat)+Clitterhr(:,:,indlandlat)+Csoilhr(:,:,indlandlat))* \
                          conform(Cplanthr(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/)))
               CsoilTTlat(isim,ilat,iyr-startyr(isim))=sum(Csoilage(:,:,indlandlat)*Csoilhr(:,:,indlandlat)* \
                          conform(Csoilhr(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/))) \
                         /sum((Cplanthr(:,:,indlandlat)+Clitterhr(:,:,indlandlat)+Csoilhr(:,:,indlandlat))* \
                          conform(Cplant(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/)))
               Cplantagelat(isim,ilat,iyr-startyr(isim))=sum(Cplantage(:,:,indlandlat)*Cplant(:,:,indlandlat)* \
                          conform(Cplant(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/))) \
                         /sum((Cplant(:,:,indlandlat)+Clitter(:,:,indlandlat)+Csoil(:,:,indlandlat))* \
                          conform(Cplant(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/)))
               Clitteragelat(isim,ilat,iyr-startyr(isim))=sum(Clitterage(:,:,indlandlat)*Clitter(:,:,indlandlat)* \
                          conform(Clitter(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/))) \
                         /sum((Cplant(:,:,indlandlat)+Clitter(:,:,indlandlat)+Csoil(:,:,indlandlat))* \
                          conform(Cplant(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/)))
               Csoilagelat(isim,ilat,iyr-startyr(isim))=sum(Csoilage(:,:,indlandlat)*Csoil(:,:,indlandlat)* \
                          conform(Csoil(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/))) \
                         /sum((Cplant(:,:,indlandlat)+Clitter(:,:,indlandlat)+Csoil(:,:,indlandlat))* \
                          conform(Cplant(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/)))
               Cturnoveragelat(isim,ilat,iyr-startyr(isim))=sum((Csoil(:,:,indlandlat)+Clitter(:,:,indlandlat)+Cplant(:,:,indlandlat))* \
                          conform(Csoil(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/))) \
                         /sum((Clitterhr(:,:,indlandlat)+Csoilhr(:,:,indlandlat))* \
                          conform(Cplant(:,:,indlandlat),conform(patchfrac(:,indlandlat),area(indlandlat),1)*patchfrac(:,indlandlat),(/1,2/)))
            end if
         end do
      end do
   end do
   asciiwrite("Ctotalagelat.txt",Cplantagelat+Clitteragelat+Csoilagelat)
   asciiwrite("Cplantagelat.txt",Cplantagelat)
   asciiwrite("Clitteragelat.txt",Clitteragelat)
   asciiwrite("Csoilagelat.txt",Csoilagelat)
   asciiwrite("CtotalTTlat.txt",ClitterTTlat+CsoilTTlat)
   asciiwrite("ClitterTTlat.txt",ClitterTTlat)
   asciiwrite("CsoilTTlat.txt",CsoilTTlat)
   asciiwrite("Cturnoveragelat.txt",Cturnoveragelat)
end if
   print(Cplantagelat(0,:,0)+Clitteragelat(0,:,0)+Csoilagelat(0,:,0))
   print("globalTT")
   print(ClitterTTlat(0,:,0)+CsoilTTlat(0,:,0))
   print(fgrid->latitude)

  Cturnoveragelat(:,0:17,:)=Cturnoveragelat@_FillValue
  Cturnoveragelat(:,88:,:)=Cturnoveragelat@_FillValue

data14C=asciiread("C14model/HeC14.txt",(/nC14,8/),"float")
lat14C=data14C(:,2)
lon14C=data14C(:,3)
C14age=data14C(:,7)

wks = gsn_open_wks("pdf","TTMA_latitudinal")
gsn_define_colormap(wks,"rainbow+white+gray")    ; Change color map
  res                          = True

  res@tiMainFontHeightF                = 0.025
  res@tiXAxisString            = "Carbon mean age (year)"
  res@tiYAxisString            = "Latitude"
  res@xyLineColors             = (/"black","grey26","black"/)
  res@xyLineThicknesses        = (/3,3,3/)
;  res@trXMinF                  = 1
;  res@trXMaxF                  = nyr
;  res@trYMinF                   = 0
  res@tmLabelAutoStride                = True
  res@xyXStyle = "Log"
;  res@tmXBMode = "Explicit"
;  res@tmXBValues = (/10,100,1000/)
;  res@tmXBLables = (/"10","100","1000"/)

  latitude_new=fgrid->latitude

  print(datamodel(:,1))
  print(latitude_new)
  dataob_int=linint1(dataob(:,1),dataob(:,0),False,latitude_new,0)
  datamodel_int=linint1(datamodel(:,1),datamodel(:,0),False,latitude_new,0)
  dataob_lb_int=linint1(dataob_lb(:,1),dataob_lb(:,0),False,latitude_new,0)
  dataob_rb_int=linint1(dataob_rb(:,1),dataob_rb(:,0),False,latitude_new,0)
  latitude_new(0:17)=latitude_new@_FillValue
  latitude_new(91:95)=latitude_new@_FillValue
  print(datamodel_int)
  res@trXMinF                   = 100
  res@trXMaxF                   = 10000
print("here0")
  res@tiMainString = "Latitudinal variation in C Mean age in 1900s"
print("here1")
  res@xyLineColors := (/"black","black","black","black"/)
  plot_xy=gsn_csm_xy(wks,(/dim_avg(Cplantagelat(3,:,81:104)+Clitteragelat(3,:,81:104)+Csoilagelat(3,:,81:104))/),latitude_new,res)
  plot3=plot_xy
print("here2")
  asciiwrite("meanagelat1982-2005.txt",dim_avg(Cplantagelat(3,:,81:104)+Clitteragelat(3,:,81:104)+Csoilagelat(3,:,81:104)))
  res@tiMainString = "Latitudinal variation in C Mean age in 2090s"
  res@xyLineColors = (/"black","black","black","black"/)
  res@xyDashPatterns = (/0,2,16/)
  plot_xy=gsn_csm_xy(wks,(/dim_avg(Cplantagelat(3,:,190:199)+Clitteragelat(3,:,190:199)+Csoilagelat(3,:,190:199)) \
                          ,dim_avg(Cplantagelat(4,:,190:199)+Clitteragelat(4,:,190:199)+Csoilagelat(4,:,190:199)) \ 
                          ,dim_avg(Cplantagelat(5,:,190:199)+Clitteragelat(5,:,190:199)+Csoilagelat(5,:,190:199))/),latitude_new,res)
  plot4=plot_xy
  res@tiMainString = "Latitudinal variation in C Mean age in 2090s"
  plot_xy=gsn_csm_xy(wks,dim_avg(Cplantagelat(3,:,190:199)+Clitteragelat(3,:,190:199)+Csoilagelat(3,:,190:199)),latitude_new,res)
  res@tiMainString = "Latitudinal variation in C Mean age change from 1900s to 2090s"
  res@xyXStyle = "Linear"
  res@trXMinF                   = -240
  res@trXMaxF                   = 0
  plot_xy=gsn_csm_xy(wks,(dim_avg(Cplantagelat(3,:,190:199)+Clitteragelat(3,:,190:199)+Csoilagelat(3,:,190:199)) \
                         -dim_avg(Cplantagelat(3,:,0:9)+Clitteragelat(3,:,0:9)+Csoilagelat(3,:,0:9))),latitude_new,res)

  res@trXMinF                  := -0.5
  res@trXMaxF                  = 1
  res@tiXAxisString = "Mean age relative change"
  plot_xy=gsn_csm_xy(wks,((/dim_avg(Cplantagelat(3,:,190:199)+Clitteragelat(3,:,190:199)+Csoilagelat(3,:,190:199)), \
                            dim_avg(Cplantagelat(4,:,190:199)+Clitteragelat(4,:,190:199)+Csoilagelat(4,:,190:199)), \
                            dim_avg(Cplantagelat(5,:,190:199)+Clitteragelat(5,:,190:199)+Csoilagelat(5,:,190:199))/) \
                         -conform_dims((/3,nlat/),dim_avg(Cplantagelat(3,:,0:9)+Clitteragelat(3,:,0:9)+Csoilagelat(3,:,0:9)),1)) \
                         /conform_dims((/3,nlat/),dim_avg(Cplantagelat(3,:,0:9)+Clitteragelat(3,:,0:9)+Csoilagelat(3,:,0:9)),1),latitude_new,res)
  res@tiXAxisString = "Mean age absolute change"
  res@trXMinF                   = -500
  res@trXMaxF                   = 500
;  plot_xy=gsn_csm_xy(wks,((/dim_avg(Cplantagelat(3,:,190:199)+Clitteragelat(3,:,190:199)+Csoilagelat(3,:,190:199)), \
;                            dim_avg(Cplantagelat(4,:,190:199)+Clitteragelat(4,:,190:199)+Csoilagelat(4,:,190:199)), \
;                            dim_avg(Cplantagelat(5,:,190:199)+Clitteragelat(5,:,190:199)+Csoilagelat(5,:,190:199))/) \
;                         -conform_dims((/3,nlat/),dim_avg(Cplantagelat(3,:,0:9)+Clitteragelat(3,:,0:9)+Csoilagelat(3,:,0:9)),1)),latitude_new,res)
  plot4=plot_xy

  plot_xy=gsn_csm_xy(wks,(dim_avg(Cplantagelat(3,:,190:199)+Clitteragelat(3,:,190:199)+Csoilagelat(3,:,190:199)) \
                         /dim_avg(Cplantagelat(3,:,0:9)+Clitteragelat(3,:,0:9)+Csoilagelat(3,:,0:9))),latitude_new,res)


  res@xyXStyle = "Log"
  res@trXMinF                   = 4
  res@trXMaxF                   = 1000
  res@tiMainString = "Latitudinal variation in C Transit time"
  res@tiXAxisString            = "Transit time / Turnover time (year)"

  plot_xy=gsn_csm_xy(wks,(/dim_avg(ClitterTTlat(3,:,81:104)+CsoilTTlat(3,:,81:104)),10^dataob_int,10^datamodel_int/),latitude_new,res)
  plot_xy=gsn_csm_xy(wks,(/dim_avg(ClitterTTlat(3,:,81:104)+CsoilTTlat(3,:,81:104)),10^dataob_int/),latitude_new,res)
  ;plot_xy=gsn_csm_xy(wks,(/dim_avg(ClitterTTlat(3,:,0:9)+CsoilTTlat(3,:,0:9)),dataob_int,datamodel_int/),latitude_new,res)

  gsres                   = True                        ; poly res
  gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres@gsFillColor       = "grey75"                 ; color chosen
  print("here3")
  print(10^dataob_lb(:,0))
  print("here4")
  print(10^dataob_rb(::-1,0))
  
  tmp=new((/nob_lb+nob_rb,2/),float)
  tmp(0:nob_lb-1,:)=dataob_lb(:,:)
  tmp(nob_lb:,:)=dataob_rb(::-1,:)
  dummy = gsn_add_polygon (wks,plot_xy,10^tmp(:,0),tmp(:,1),gsres)
;  dummy = gsn_add_polygon (wks,plot_xy,tmp(:,0),tmp(:,1),gsres)
  print(tmp(:,0))
  print(tmp(:,1))
  print(dataob_int)
  draw(plot_xy)                                            ; draw
  frame(wks)                                            ; advance frame


;  plot_xy_new1=gsn_csm_xy(wks,dataob(:,0),dataob(:,1),res)
;  res@xyLineColors = (/"red","red","red"/)
;  plot_xy_new2=gsn_csm_xy(wks,datamodel(:,0),datamodel(:,1),res)
;  overlay(plot_xy,plot_xy_new2)
;  draw(plot_xy)
;  frame(wks)
;  res@xyLineColors             = (/"blue","red","green"/)
  plot1=plot_xy
  res@tiXAxisString            = "Transit time (year)"
  res@xyLineColors := (/"black","black","black","black"/)
  res@tiMainString = "Latitudinal variation in C Transit time in 2010s"
  plot1=gsn_csm_xy(wks,dim_avg(ClitterTTlat(3,:,0:9)+CsoilTTlat(3,:,0:9)),latitude_new,res)
  plot_xy=gsn_csm_xy(wks,dim_avg(ClitterTTlat(3,:,110:119)+CsoilTTlat(3,:,110:119)),latitude_new,res)
  res@xyLineColors = (/"black","black","black","black"/)
  res@xyDashPatterns = (/0,2,16/)
  res@tiMainString = "Latitudinal variation in C Transit time in 2090s"
  plot_xy=gsn_csm_xy(wks,(/dim_avg(ClitterTTlat(3,:,190:199)+CsoilTTlat(3,:,190:199))\
                          ,dim_avg(ClitterTTlat(4,:,190:199)+CsoilTTlat(4,:,190:199))\
                          ,dim_avg(ClitterTTlat(5,:,190:199)+CsoilTTlat(5,:,190:199))/),latitude_new,res)
  plot2=plot_xy
  res@xyXStyle = "Linear"
  res@trXMinF                  = -12
  res@trXMaxF                  = 6
  res@tiMainString = "Latitudinal variation in C Transit time change from 1900s to 2090s"
  plot_xy=gsn_csm_xy(wks,dim_avg(ClitterTTlat(3,:,190:199)+CsoilTTlat(3,:,190:199))\
                        -dim_avg(ClitterTTlat(3,:,0:9)+CsoilTTlat(3,:,0:9)),latitude_new,res)

  res@trXMinF                  = -0.5
  res@trXMaxF                  = 1
  res@tiMainString = "Latitudinal variation of transit time and mean age relative change"
  res@tiXAxisString = "Transit time relative change"
  plot_xy=gsn_csm_xy(wks,((/dim_avg(ClitterTTlat(3,:,190:199)+CsoilTTlat(3,:,190:199)),\
                            dim_avg(ClitterTTlat(4,:,190:199)+CsoilTTlat(4,:,190:199)),\
                            dim_avg(ClitterTTlat(5,:,190:199)+CsoilTTlat(5,:,190:199))/)\
                         -conform_dims((/3,nlat/),dim_avg(ClitterTTlat(3,:,0:9)+CsoilTTlat(3,:,0:9)),1)) \
                         /conform_dims((/3,nlat/),dim_avg(ClitterTTlat(3,:,0:9)+CsoilTTlat(3,:,0:9)),1) ,latitude_new,res)
  res@tiXAxisString = "Transit time absolute change"
  res@trXMinF                  = -30
  res@trXMaxF                  = 18
  plot_xy=gsn_csm_xy(wks,((/dim_avg(ClitterTTlat(3,:,190:199)+CsoilTTlat(3,:,190:199)),\
                            dim_avg(ClitterTTlat(4,:,190:199)+CsoilTTlat(4,:,190:199)),\
                            dim_avg(ClitterTTlat(5,:,190:199)+CsoilTTlat(5,:,190:199))/)\
                         -conform_dims((/3,nlat/),dim_avg(ClitterTTlat(3,:,0:9)+CsoilTTlat(3,:,0:9)),1)),latitude_new,res)
 plot2=plot_xy

  res@xyXStyle = "Log"
  res@trXMinF                   = 4
  res@trXMaxF                   = 1000
  res@tiMainString = "Latitudinal variation in C turnover time"
  res@tiXAxisString            = "Turnover time (year)"
  print(Cturnoveragelat(3,:,0:9))
  plot_xy=gsn_csm_xy(wks,dim_avg(Cturnoveragelat(3,:,0:9)),latitude_new,res)

  res@xyXStyle = "Linear"
  res@trXMinF                  = -0.5
  res@trXMaxF                  = 1
  res@tiMainString = "Latitudinal variation of turnover time  change"
  res@tiXAxisString = "Turnover time relative change"
  plot_xy=gsn_csm_xy(wks,((/dim_avg(Cturnoveragelat(3,:,190:199)),\
                            dim_avg(Cturnoveragelat(4,:,190:199)),\
                            dim_avg(Cturnoveragelat(5,:,190:199))/)\
                         -conform_dims((/3,nlat/),dim_avg(Cturnoveragelat(3,:,0:9)),1)) \
                         /conform_dims((/3,nlat/),dim_avg(Cturnoveragelat(3,:,0:9)),1) ,latitude_new,res)
  res@tiXAxisString = "Turnover time absolute change"
  res@trXMinF                  = -30
  res@trXMaxF                  = 18
  plot_xy=gsn_csm_xy(wks,((/dim_avg(Cturnoveragelat(3,:,190:199)),\
                            dim_avg(Cturnoveragelat(4,:,190:199)),\
                            dim_avg(Cturnoveragelat(5,:,190:199))/)\
                         -conform_dims((/3,nlat/),dim_avg(Cturnoveragelat(3,:,0:9)),1)),latitude_new,res)
  plot2=plot_xy

  plot_xy=gsn_csm_xy(wks,(dim_avg(ClitterTTlat(3,:,190:199)+CsoilTTlat(3,:,190:199))\
                         /dim_avg(ClitterTTlat(3,:,0:9)+CsoilTTlat(3,:,0:9))),latitude_new,res)

gsn_panel(wks,(/plot1,plot2,plot3,plot4/),(/2,2/),False)


  
  res@trXMinF                   = 20
  res@trXMaxF                   = 30000
  res@tiMainString = "Latitudinal variation in soil C Mean age in 1900s"
  res@xyLineColors := (/"black","gray25","black","black"/)
  res@gsnFrame = False
  res@gsnDraw  = False
  res@tmXBLabelFontHeightF = 0.03
  res@tmYLLabelFontHeightF = 0.03
  res@tmXBLabelFont = 21
  res@tmYLLabelFont = 21
  res@xyXStyle = "Log"
  res@xyLineThicknesses        = (/6,6,6/)
  res@xyDashPatterns  = (/0,16,16/)
  y=latitude_new
  fittingage = 8.e-8 * y^6 - 9.e-6 * y^5 +0.0005 * y^4 -0.0052 * y^3 -0.9432 * y^2 + 34.243 * y + 2112.9
  print("here0")
  fittingage = where( y .gt. max(lat14C) .or. y .lt. min(lat14C),Csoilagelat@_FillValue,fittingage)
  print("here1")
  print(fittingage)
  plot_xy=gsn_csm_xy(wks,(/dim_avg(Csoilagelat(3,:,57:110)),fittingage/),latitude_new,res)
  print("here2")
  
  
 res1                     = True                   ; plot mods desired
 res1@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
 res1@xyMarkers           = 16                     ; choose type of marker
 res1@xyMarkerColor       = "gray65"                  ; Marker color
 res1@xyLineColor         = "gray25"                  ; Marker color
 res1@xyMonoLineColor     = True
 res1@xyMonoMarkerColor   = True
 res1@xyMarkerSizeF       = 0.010                  ; Marker size (default 0.01)
 res1@xyDashPattern       = 16
 res1@xyMonoDashPattern   = True
 res1@xyLineThicknessF    = 10               ; set second line to 2
 res1@xyMonoLineThickness = True
 res1@tmXBLabelFontHeightF = 0.03
 res1@tmYLLabelFontHeightF = 0.03
 res1@tmXBLabelFont = 21
 res1@tmYLLabelFont = 21
 res1@trXMinF = 20.0
 res1@trXMaxF = 30000
 res1@gsnFrame = False
 res1@gsnDraw  = False
 res1@xyXStyle = "Log"
 res1@trYMaxF = 90

 lat14C = where (C14age .le. 0, lat14C@_FillValue,lat14C)
 lat14C@units = "degrees_north"
 lat14C@long_name = "Latitude"
 C14age = where (C14age .le. 0, C14age@_FillValue,C14age)
 print(C14age)
 print(lat14C)
 plot1 = gsn_csm_xy(wks,C14age,lat14C,res1)
 
 print("here3")

 overlay(plot1,plot_xy)
 print("here4")
 draw(plot1)
 print("here5")
 frame(wks)
  

end
