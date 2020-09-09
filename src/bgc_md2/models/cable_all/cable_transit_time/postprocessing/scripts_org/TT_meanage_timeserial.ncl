load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
startyr=1901
endyr=2100
nyr=(endyr-startyr)+1
nland=5656
pftind=(/1,2,3,4,5,6,7,8,9,14/)
PFTName=(/"ENF","EBF","DNF","DBF","SHRUB","C3","C4","TUNDRA","CROP","BARREN"/)
npft=dimsizes(pftind)
;SimName=(/"C/S1","C/S2","C/S3","CN/S1","CN/S2","CN/S3"/)
SimName=(/"CN/S1","CN/S2","CN/S3"/)
nsim=dimsizes(SimName)
FilePath="/datastore/lu02b/CABLE-traceability/"

Cplantagepft=new((/nsim,npft,3,nyr/),float)
Clitteragepft=new((/nsim,npft,3,nyr/),float)
Csoilagepft=new((/nsim,npft,3,nyr/),float)
Ctotagepft=new((/nsim,npft,nyr/),float)
Cplantageglobal=new((/nsim,3,nyr/),float)
Clitterageglobal=new((/nsim,3,nyr/),float)
Csoilageglobal=new((/nsim,3,nyr/),float)
Ctotageglobal=new((/nsim,nyr/),float)

CplantTTpft=new((/nsim,npft,3,nyr/),float)
ClitterTTpft=new((/nsim,npft,3,nyr/),float)
CsoilTTpft=new((/nsim,npft,3,nyr/),float)
CtotTTpft=new((/nsim,npft,nyr/),float)
CplantTTglobal=new((/nsim,3,nyr/),float)
ClitterTTglobal=new((/nsim,3,nyr/),float)
CsoilTTglobal=new((/nsim,3,nyr/),float)
CtotTTglobal=new((/nsim,nyr/),float)

Cplanthrpft=new((/nsim,npft,3,nyr/),float)
Clitterhrpft=new((/nsim,npft,3,nyr/),float)
Csoilhrpft=new((/nsim,npft,3,nyr/),float)
Ctothrpft=new((/nsim,npft,nyr/),float)
Cplanthrglobal=new((/nsim,3,nyr/),float)
Clitterhrglobal=new((/nsim,3,nyr/),float)
Csoilhrglobal=new((/nsim,3,nyr/),float)
Ctothrglobal=new((/nsim,nyr/),float)

Cplantpft=new((/nsim,npft,3,nyr/),float)
Clitterpft=new((/nsim,npft,3,nyr/),float)
Csoilpft=new((/nsim,npft,3,nyr/),float)

Ctotcwtpft=new((/nsim,npft,nyr/),float)
Ctotcwtglobal=new((/nsim,nyr/),float)

FilePath1="/datastore/lu02b/CABLE-traceability/CN/S2/annual_TT/"
FilePath2="/datastore/lu02b/CABLE-traceability/CN/S2/annual/"

fgrid=addfile("gridinfo_NCAR_1.9x2.5_landfrac_revised.nc","r")
finfo=addfile("/datastore/lu02b/CABLE-traceability/C/S1/output/restart_ncar_2101_ndep.nc","r")
latitude=finfo->latitude
longitude=where(finfo->longitude .lt. 0,finfo->longitude+360,finfo->longitude)

indlat=ind_nearest_coord(latitude,fgrid->latitude,0)
indlon=ind_nearest_coord(longitude,fgrid->longitude,0)
landfrac=new(nland,float)
area=new(nland,float)

;print(latitude)
do iland=0,nland-1
   landfrac(iland)=fgrid->land_fraction(indlat(iland),indlon(iland))
   area(iland)=fgrid->area(indlat(iland),indlon(iland))*landfrac(iland)
end do

isDataReady=True
if(isDataReady)then
   Cplantagepft=asciiread("Cplantagepft.txt",(/nsim,npft,3,nyr/),"float")
   Clitteragepft=asciiread("Clitteragepft.txt",(/nsim,npft,3,nyr/),"float")
   Csoilagepft=asciiread("Csoilagepft.txt",(/nsim,npft,3,nyr/),"float")
   Cplantageglobal=asciiread("Cplantageglobal.txt",(/nsim,3,nyr/),"float")
   Clitterageglobal=asciiread("Clitterageglobal.txt",(/nsim,3,nyr/),"float")
   Csoilageglobal=asciiread("Csoilageglobal.txt",(/nsim,3,nyr/),"float")
   Ctotageglobal=asciiread("Ctotageglobal.txt",(/nsim,nyr/),"float")

   CplantTTpft=asciiread("CplantTTpft.txt",(/nsim,npft,3,nyr/),"float")
   ClitterTTpft=asciiread("ClitterTTpft.txt",(/nsim,npft,3,nyr/),"float")
   CsoilTTpft=asciiread("CsoilTTpft.txt",(/nsim,npft,3,nyr/),"float")
   Ctotpft=asciiread("CtotTTpft.txt",(/nsim,npft,nyr/),"float")
   CplantTTglobal=asciiread("CplantTTglobal.txt",(/nsim,3,nyr/),"float")
   ClitterTTglobal=asciiread("ClitterTTglobal.txt",(/nsim,3,nyr/),"float")
   CsoilTTglobal=asciiread("CsoilTTglobal.txt",(/nsim,3,nyr/),"float")
   CtotTTglobal=asciiread("CtotTTglobal.txt",(/nsim,nyr/),"float")

   Ctotcwtpft=asciiread("Ctotcwtpft.txt",(/nsim,npft,nyr/),"float")
   Ctotcwtglobal=asciiread("Ctotcwtglobal.txt",(/nsim,nyr/),"float")

   Cplanthrpft=asciiread("Cplanthrpft.txt",(/nsim,npft,3,nyr/),"float")
   Clitterhrpft=asciiread("Clitterhrpft.txt",(/nsim,npft,3,nyr/),"float")
   Csoilhrpft=asciiread("Csoilhrpft.txt",(/nsim,npft,3,nyr/),"float")
   Ctotpft=asciiread("Ctothrpft.txt",(/nsim,npft,nyr/),"float")
   Cplanthrglobal=asciiread("Cplanthrglobal.txt",(/nsim,3,nyr/),"float")
   Clitterhrglobal=asciiread("Clitterhrglobal.txt",(/nsim,3,nyr/),"float")
   Csoilhrglobal=asciiread("Csoilhrglobal.txt",(/nsim,3,nyr/),"float")
   Ctothrglobal=asciiread("Ctothrglobal.txt",(/nsim,nyr/),"float")

   NPPglobal = asciiread("NPPCNtotal.txt",(/3,nsim,nyr/),"float")
   HRglobal  = asciiread("HeteroRespCNtotal.txt",(/3,nsim,nyr/),"float")
   Cecoglobal=asciiread("CecosystemCN.txt",(/3,nsim,nyr/),"float")
   
else
   do isim=0,nsim-1
      print((/SimName(isim)/))
      do iyr=startyr,endyr
         print((/iyr/))
         finage=addfile(FilePath+SimName(isim)+"/annual_TT/out_ncar_"+iyr+"_annual_age.nc","r")
         fincpool=addfile(FilePath+SimName(isim)+"/annual/out_ncar_"+iyr+"_annual.nc","r")
         finhr=addfile(FilePath+SimName(isim)+"/annual_TT/out_ncar_"+iyr+"_annualTT.nc","r")
         if(iyr .eq. startyr)then
            iveg=fincpool->iveg
            patchfrac=fincpool->patchfrac
            patchfrac=where(iveg .eq. 17,patchfrac@_FillValue,patchfrac)
            ;print(sum(conform(patchfrac,area,1)*patchfrac))
         end if
         Cplantage=finage->Cplant
         Clitterage=finage->Clitter
         Csoilage=finage->Csoil
         Cplant=fincpool->Cplant
         Clitter=fincpool->Clitter
         Csoil=fincpool->Csoil
         kplant=fincpool->kplant
         klitter=fincpool->klitter
         ksoil=fincpool->ksoil
         LeaftoL=dim_sum_n(finage->fromLeaftoL,0)
         WoodtoL=dim_sum_n(finage->fromWoodtoL,0)
         RoottoL=dim_sum_n(finage->fromRoottoL,0)
         MettoS=dim_sum_n(finage->fromMettoS,0)
         StrtoS=dim_sum_n(finage->fromStrtoS,0)
         CWDtoS=dim_sum_n(finage->fromCWDtoS,0)
         FasttoS=finage->fromSOMtoSOM(0,:,:)+finage->fromSOMtoSOM(1,:,:)
         SlowtoS=finage->fromSOMtoSOM(2,:,:)
         PasstoS=SlowtoS
         PasstoS=0
     

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

         Cplantto=Cplanthr+(/LeaftoL,RoottoL,WoodtoL/)
         Cplantto=where(Cplantto .eq. 0,Cplantto@_FillValue,Cplantto)
         Clitterto=Clitterhr+(/MettoS,StrtoS,CWDtoS/)
         Clitterto=where(Clitterto .eq. 0,Clitterto@_FillValue,Clitterto)
         Csoilto=Csoilhr+(/FasttoS,SlowtoS,PasstoS/)
         Csoilto=where(Csoilto .eq. 0,Csoilto@_FillValue,Csoilto)

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

         do ipft=0,npft-1
            print((/PFTName(ipft)/))
            patchfractmp=where(iveg .eq. pftind(ipft),patchfrac,patchfrac@_FillValue)
            asciiwrite(PFTName(ipft)+"_area.txt",sum(conform(patchfractmp,area,1)*patchfractmp))
;            print(sum(conform(patchfractmp,area,1)*patchfractmp))
            Cplantagepft(isim,ipft,:,iyr-startyr)=dim_sum_n(Cplantage*Cplant*conform(Cplant,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/)) \
                                                        /dim_sum_n(Cplant*conform(Cplant,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/))
            CplantTTpft(isim,ipft,:,iyr-startyr)=0
            Cplanthrpft(isim,ipft,:,iyr-startyr)=0
            Cplantpft(isim,ipft,:,iyr-startyr)=dim_sum_n(Cplant*conform(Cplant,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/))*1.e-15

            Clitteragepft(isim,ipft,:,iyr-startyr)=dim_sum_n(Clitterage*Clitter*conform(Clitter,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/)) \
                                                         /dim_sum_n(Clitter*conform(Clitter,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/))
            ClitterTTpft(isim,ipft,:,iyr-startyr)=dim_sum_n(Clitterage*Clitterhr*conform(Clitter,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/)) \
                                                         /dim_sum_n(Clitterhr*conform(Clitterhr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/))
            Clitterhrpft(isim,ipft,:,iyr-startyr)=dim_sum_n(Clitterhr*conform(Clitterhr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/)) \
                                                         /sum(conform(patchfractmp,area,1)*patchfractmp)
            Clitterpft(isim,ipft,:,iyr-startyr)=dim_sum_n(Clitter*conform(Clitter,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/))*1.e-15
            Csoilagepft(isim,ipft,:,iyr-startyr)=dim_sum_n(Csoilage*Csoil*conform(Csoil,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/)) \
                                                       /dim_sum_n(Csoil*conform(Csoil,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/))
            CsoilTTpft(isim,ipft,:,iyr-startyr)=dim_sum_n(Csoilage*Csoilhr*conform(Csoil,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/)) \
                                                       /dim_sum_n(Csoilhr*conform(Csoil,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/))
            Csoilhrpft(isim,ipft,:,iyr-startyr)=dim_sum_n(Csoilhr*conform(Csoilhr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/)) \
                                                       /sum(conform(patchfractmp,area,1)*patchfractmp)
            Csoilpft(isim,ipft,:,iyr-startyr)=dim_sum_n(Csoil*conform(Csoil,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/))*1.e-15
            Ctotagepft(isim,ipft,iyr-startyr) = sum(Cplantage*Cplant*conform(Cplant,conform(patchfractmp,area,1)*patchfractmp,(/1,2/))    \ 
                                                   +Clitterage*Clitter*conform(Clitter,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)) \
                                                   +Csoilage*Csoil*conform(Csoil,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)))      \
                                                   /sum((Cplant+Clitter+Csoil)*conform(Csoil,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)))
            CtotTTpft(isim,ipft,iyr-startyr) = (sum(Cplantage*Cplanthr*conform(Cplanthr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/))) \
                                               +sum(Clitterage*Clitterhr*conform(Clitterhr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/))) \
                                               +sum(Csoilage*Csoilhr*conform(Csoilhr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)))) \
                                               /sum((Cplanthr+Clitterhr+Csoilhr)*conform(Csoilhr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)))
            Ctothrpft(isim,ipft,iyr-startyr) =(sum(Cplanthr*conform(Cplanthr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/))) \
                                              +sum(Clitterhr*conform(Clitterhr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/))) \
                                              +sum(Csoilhr*conform(Csoilhr,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)))) \
                                              /sum(conform(patchfractmp,area,1)*patchfractmp)
            ;print(dim_sum_n(Cplant*conform(Cplant,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/)))
            ;print(dim_sum_n(Clitter*conform(Clitter,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)),(/1,2/)))
            Ctotcwtpft(isim,ipft,iyr-startyr) =(sum(Cplant*Cplant/Cplantto*conform(Cplant,conform(patchfractmp,area,1)*patchfractmp,(/1,2/))) \
                                               +sum(Clitter*Clitter/Clitterto*conform(Clitter,conform(patchfractmp,area,1)*patchfractmp,(/1,2/))) \
                                               +sum(Csoil*Csoil/Csoilto*conform(Csoil,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)))) \
                                               /sum((Cplant+Clitter+Csoil)*conform(Csoil,conform(patchfractmp,area,1)*patchfractmp,(/1,2/)))
         end do
         Cplantageglobal(isim,:,iyr-startyr)  = dim_sum_n(Cplantage*Cplant*conform(Cplant,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/)) \
                                                  /dim_sum_n(Cplant*conform(Cplant,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/))
         Clitterageglobal(isim,:,iyr-startyr) = dim_sum_n(Clitterage*Clitter*conform(Clitter,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/)) \
                                                  /dim_sum_n(Clitter*conform(Clitter,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/))
         Csoilageglobal(isim,:,iyr-startyr) = dim_sum_n(Csoilage*Csoil*conform(Csoil,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/)) \
                                                  /dim_sum_n(Csoil*conform(Csoil,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/))
         Ctotageglobal(isim,iyr-startyr) = (sum(Cplantage*Cplant*conform(Cplant,conform(patchfrac,area,1)*patchfrac,(/1,2/))) \
                                               +sum(Clitterage*Clitter*conform(Clitter,conform(patchfrac,area,1)*patchfrac,(/1,2/))) \
                                               +sum(Csoilage*Csoil*conform(Csoil,conform(patchfrac,area,1)*patchfrac,(/1,2/)))) \
                                               /sum((Cplant+Clitter+Csoil)*conform(Csoil,conform(patchfrac,area,1)*patchfrac,(/1,2/)))
         CplantTTglobal(isim,:,iyr-startyr)  = 0
         ClitterTTglobal(isim,:,iyr-startyr) = dim_sum_n(Clitterage*Clitterhr*conform(Clitterhr,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/)) \
                                                  /sum((Clitterhr+Csoilhr)*conform(Clitterhr,conform(patchfrac,area,1)*patchfrac,(/1,2/)))
         CsoilTTglobal(isim,:,iyr-startyr) = dim_sum_n(Csoilage*Csoilhr*conform(Csoilhr,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/)) \
                                                  /sum((Clitterhr+Csoilhr)*conform(Csoilhr,conform(patchfrac,area,1)*patchfrac,(/1,2/)))
         CtotTTglobal(isim,iyr-startyr) =(sum(Cplantage*Cplanthr*conform(Cplanthr,conform(patchfrac,area,1)*patchfrac,(/1,2/))) \
                                               +sum(Clitterage*Clitterhr*conform(Clitterhr,conform(patchfrac,area,1)*patchfrac,(/1,2/))) \
                                               +sum(Csoilage*Csoilhr*conform(Csoilhr,conform(patchfrac,area,1)*patchfrac,(/1,2/)))) \
                                               /sum((Cplanthr+Clitterhr+Csoilhr)*conform(Csoilhr,conform(patchfrac,area,1)*patchfrac,(/1,2/)))
         Cplanthrglobal(isim,:,iyr-startyr)  = 0
         Clitterhrglobal(isim,:,iyr-startyr) = dim_sum_n(Clitterhr*conform(Clitterhr,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/)) \
                                                  /sum(conform(patchfrac,area,1)*patchfrac)
         Csoilhrglobal(isim,:,iyr-startyr) = dim_sum_n(Csoilhr*conform(Csoilhr,conform(patchfrac,area,1)*patchfrac,(/1,2/)),(/1,2/)) \
                                                  /sum(conform(patchfrac,area,1)*patchfrac)
         Ctothrglobal(isim,iyr-startyr) =(sum(Cplanthr*conform(Cplanthr,conform(patchfrac,area,1)*patchfrac,(/1,2/))) \
                                         +sum(Clitterhr*conform(Clitterhr,conform(patchfrac,area,1)*patchfrac,(/1,2/))) \
                                         +sum(Csoilhr*conform(Csoilhr,conform(patchfrac,area,1)*patchfrac,(/1,2/)))) \
                                         /sum(conform(patchfrac,area,1)*patchfrac)
         Ctotcwtglobal(isim,iyr-startyr) = (sum(Cplant*Cplant/Cplantto*conform(Cplant,conform(patchfrac,area,1)*patchfrac,(/1,2/))) \
                                           +sum(Clitter*Clitter/Clitterto*conform(Clitter,conform(patchfrac,area,1)*patchfrac,(/1,2/))) \  
                                           +sum(Csoil*Csoil/Csoilto*conform(Csoil,conform(patchfrac,area,1)*patchfrac,(/1,2/)))) \
                                           /sum((Cplant+Clitter+Csoil)*conform(Csoil,conform(patchfrac,area,1)*patchfrac,(/1,2/)))
      end do
   end do
   asciiwrite("Cplantagepft.txt",Cplantagepft)
   asciiwrite("Clitteragepft.txt",Clitteragepft)
   asciiwrite("Csoilagepft.txt",Csoilagepft)
   asciiwrite("Ctotagepft.txt",Ctotagepft)
   asciiwrite("Cplantageglobal.txt",Cplantageglobal)
   asciiwrite("Clitterageglobal.txt",Clitterageglobal)
   asciiwrite("Csoilageglobal.txt",Csoilageglobal)
   asciiwrite("Ctotageglobal.txt",Ctotageglobal)

   asciiwrite("CplantTTpft.txt",CplantTTpft)
   asciiwrite("ClitterTTpft.txt",ClitterTTpft)
   asciiwrite("CsoilTTpft.txt",CsoilTTpft)
   asciiwrite("CtotTTpft.txt",CtotTTpft)
   asciiwrite("CplantTTglobal.txt",CplantTTglobal)
   asciiwrite("ClitterTTglobal.txt",ClitterTTglobal)
   asciiwrite("CsoilTTglobal.txt",CsoilTTglobal)
   asciiwrite("CtotTTglobal.txt",CtotTTglobal)

   asciiwrite("Ctotcwtpft.txt",Ctotcwtpft)
   asciiwrite("Ctotcwtglobal.txt",Ctotcwtglobal)

   asciiwrite("Cplanthrpft.txt",Cplanthrpft)
   asciiwrite("Clitterhrpft.txt",Clitterhrpft)
   asciiwrite("Csoilhrpft.txt",Csoilhrpft)
   asciiwrite("Ctothrpft.txt",Ctothrpft)
   asciiwrite("Cplanthrglobal.txt",Cplanthrglobal)
   asciiwrite("Clitterhrglobal.txt",Clitterhrglobal)
   asciiwrite("Csoilhrglobal.txt",Csoilhrglobal)
   asciiwrite("Ctothrglobal.txt",Ctothrglobal)
end if
      

wks = gsn_open_wks("pdf","TTmeanagepft"+startyr+"-"+endyr)
gsn_define_colormap(wks,"rainbow+white+gray")    ; Change color map
  res                          = True

  res@tiMainFontHeightF                = 0.025
  res@tiXAxisString            = "simulation years"
;  res@xyLineColors             = (/"black","black","black","blue","black","green","blue","black","green"/)
  res@xyLineColors             = (/"blue4","red4","green4","blue4","red4","green4","blue4","red4","green4"/)
;  res@xyDashPatterns           = (/0,2,16,16,16,16,2,2,2/)
  res@xyDashPatterns           = (/0,0,0,2,2,2,2,2,2/)
  res@xyLineThicknesses        = (/4.0,4.0,4.0,2.0,2.0,2.0,4.0,4.0,4.0/)
  res@trXMinF                  = startyr
  res@trXMaxF                  = endyr
;  res@trYMinF                  = -1000
;  res@trYMaxF                  = 5000
  res@tmLabelAutoStride                = True

  res@gsnMaximize      = True
  res@vpWidthF         = 0.8            ; Change the aspect ratio, but 
  res@vpHeightF        = 0.6            ; make plot as large as possible.

  time=ispan(startyr,endyr,1)

  do ipft=0,npft-1
     res@tiMainString=PFTName(ipft)+" leaf age"
     plot=gsn_csm_xy(wks,time,Cplantagepft(:,ipft,0,:),res)
     res@tiMainString=PFTName(ipft)+" root age"
     plot=gsn_csm_xy(wks,time,Cplantagepft(:,ipft,1,:),res)
     res@tiMainString=PFTName(ipft)+" wood age"
     plot=gsn_csm_xy(wks,time,Cplantagepft(:,ipft,2,:),res)
     res@tiMainString=PFTName(ipft)+" litter metabolic age"
     plot=gsn_csm_xy(wks,time,Clitteragepft(:,ipft,0,:),res)
     res@tiMainString=PFTName(ipft)+" litter structure age"
     plot=gsn_csm_xy(wks,time,Clitteragepft(:,ipft,1,:),res)
     res@tiMainString=PFTName(ipft)+" litter cwd age"
     plot=gsn_csm_xy(wks,time,Clitteragepft(:,ipft,2,:),res)
     res@tiMainString=PFTName(ipft)+" soil fast pool age"
     plot=gsn_csm_xy(wks,time,Csoilagepft(:,ipft,0,:),res)
     res@tiMainString=PFTName(ipft)+" soil slow pool age"
     plot=gsn_csm_xy(wks,time,Csoilagepft(:,ipft,1,:),res)
     res@tiMainString=PFTName(ipft)+" soil passive pool age"
     plot=gsn_csm_xy(wks,time,Csoilagepft(:,ipft,2,:),res)
     res@tiMainString=PFTName(ipft)+" mean age"
     plot=gsn_csm_xy(wks,time,Ctotagepft(:,ipft,:),res)

     res@tiMainString=PFTName(ipft)+" leaf TT"
     plot=gsn_csm_xy(wks,time,CplantTTpft(:,ipft,0,:),res)
     res@tiMainString=PFTName(ipft)+" root TT"
     plot=gsn_csm_xy(wks,time,CplantTTpft(:,ipft,1,:),res)
     res@tiMainString=PFTName(ipft)+" wood TT"
     plot=gsn_csm_xy(wks,time,CplantTTpft(:,ipft,2,:),res)
     res@tiMainString=PFTName(ipft)+" litter metabolic TT"
     plot=gsn_csm_xy(wks,time,ClitterTTpft(:,ipft,0,:),res)
     res@tiMainString=PFTName(ipft)+" litter structure TT"
     plot=gsn_csm_xy(wks,time,ClitterTTpft(:,ipft,1,:),res)
     res@tiMainString=PFTName(ipft)+" litter cwd TT"
     plot=gsn_csm_xy(wks,time,ClitterTTpft(:,ipft,2,:),res)
     res@tiMainString=PFTName(ipft)+" soil fast pool TT"
     plot=gsn_csm_xy(wks,time,CsoilTTpft(:,ipft,0,:),res)
     res@tiMainString=PFTName(ipft)+" soil slow pool TT"
     plot=gsn_csm_xy(wks,time,CsoilTTpft(:,ipft,1,:),res)
     res@tiMainString=PFTName(ipft)+" soil passive pool TT"
     plot=gsn_csm_xy(wks,time,CsoilTTpft(:,ipft,2,:),res)
     res@tiMainString=PFTName(ipft)+" transient time"
     plot=gsn_csm_xy(wks,time,CtotTTpft(:,ipft,:),res)

     res@tiMainString=PFTName(ipft)+" leaf hr"
     plot=gsn_csm_xy(wks,time,Cplanthrpft(:,ipft,0,:),res)
     res@tiMainString=PFTName(ipft)+" root hr"
     plot=gsn_csm_xy(wks,time,Cplanthrpft(:,ipft,1,:),res)
     res@tiMainString=PFTName(ipft)+" wood hr"
     plot=gsn_csm_xy(wks,time,Cplanthrpft(:,ipft,2,:),res)
     res@tiMainString=PFTName(ipft)+" litter metabolic hr"
     plot=gsn_csm_xy(wks,time,Clitterhrpft(:,ipft,0,:),res)
     res@tiMainString=PFTName(ipft)+" litter structure hr"
     plot=gsn_csm_xy(wks,time,Clitterhrpft(:,ipft,1,:),res)
     res@tiMainString=PFTName(ipft)+" litter cwd hr"
     plot=gsn_csm_xy(wks,time,Clitterhrpft(:,ipft,2,:),res)
     res@tiMainString=PFTName(ipft)+" soil fast pool hr"
     plot=gsn_csm_xy(wks,time,Csoilhrpft(:,ipft,0,:),res)
     res@tiMainString=PFTName(ipft)+" soil slow pool hr"
     plot=gsn_csm_xy(wks,time,Csoilhrpft(:,ipft,1,:),res)
     res@tiMainString=PFTName(ipft)+" soil passive pool hr"
     plot=gsn_csm_xy(wks,time,Csoilhrpft(:,ipft,2,:),res)
     res@tiMainString=PFTName(ipft)+" total hr"
     plot=gsn_csm_xy(wks,time,Ctothrpft(:,ipft,:),res)

     res@tiMainString=PFTName(ipft)+" total C weighted turnover time"
     plot=gsn_csm_xy(wks,time,Ctotcwtpft(:,ipft,:),res)
  end do
  res@tiMainString="global leaf mean age"
  plot=gsn_csm_xy(wks,time,Cplantageglobal(:,0,:),res)
  res@tiMainString="global root mean age"
  plot=gsn_csm_xy(wks,time,Cplantageglobal(:,1,:),res)
  res@tiMainString="global wood mean age"
  plot=gsn_csm_xy(wks,time,Cplantageglobal(:,2,:),res)
  res@tiMainString="global litter metabolic age"
  plot=gsn_csm_xy(wks,time,Clitterageglobal(:,0,:),res)
  res@tiMainString="global litter structure age"
  plot=gsn_csm_xy(wks,time,Clitterageglobal(:,1,:),res)
  res@tiMainString="global litter passive age"
  plot=gsn_csm_xy(wks,time,Clitterageglobal(:,2,:),res)
  res@tiMainString="global soil fast pool age"
  plot=gsn_csm_xy(wks,time,Csoilageglobal(:,0,:),res)
  res@tiMainString="global soil slow pool age"
  plot=gsn_csm_xy(wks,time,Csoilageglobal(:,1,:),res)
  res@tiMainString="global soil passive pool age"
  plot=gsn_csm_xy(wks,time,Csoilageglobal(:,2,:),res)
  
  res@tiMainString="global C mean age"
  plot=gsn_csm_xy(wks,time,Ctotageglobal,res)
  plot0=plot
  
  res@tiMainString="global leaf mean transient time"
  plot=gsn_csm_xy(wks,time,CplantTTglobal(:,0,:),res)
  res@tiMainString="global root mean transient time"
  plot=gsn_csm_xy(wks,time,CplantTTglobal(:,1,:),res)
  res@tiMainString="global wood mean transient time"
  plot=gsn_csm_xy(wks,time,CplantTTglobal(:,2,:),res)
  res@tiMainString="global litter metabolic transient time"
  plot=gsn_csm_xy(wks,time,ClitterTTglobal(:,0,:),res)
  res@tiMainString="global litter structure transient time"
  plot=gsn_csm_xy(wks,time,ClitterTTglobal(:,1,:),res)
  res@tiMainString="global litter passive transient time"
  plot=gsn_csm_xy(wks,time,ClitterTTglobal(:,2,:),res)
  res@tiMainString="global soil fast pool transient time"
  plot=gsn_csm_xy(wks,time,CsoilTTglobal(:,0,:),res)
  res@tiMainString="global soil slow pool transient time"
  plot=gsn_csm_xy(wks,time,CsoilTTglobal(:,1,:),res)
  res@tiMainString="global soil passive pool transient time"
  plot=gsn_csm_xy(wks,time,CsoilTTglobal(:,2,:),res)
  
  res@tiMainString="global C transient time"
  plot=gsn_csm_xy(wks,time,CtotTTglobal,res)
  plot1=gsn_csm_xy(wks,time,reshape((/CtotTTglobal,Cecoglobal(0,:,:)/HRglobal(0,:,:)/),(/6,nyr/)),res)
  
  res@tiMainString="global leaf hr"
  plot=gsn_csm_xy(wks,time,Cplanthrglobal(:,0,:),res)
  res@tiMainString="global root hr"
  plot=gsn_csm_xy(wks,time,Cplanthrglobal(:,1,:),res)
  res@tiMainString="global wood hr"
  plot=gsn_csm_xy(wks,time,Cplanthrglobal(:,2,:),res)
  res@tiMainString="global litter metabolic hr"
  plot=gsn_csm_xy(wks,time,Clitterhrglobal(:,0,:),res)
  res@tiMainString="global litter structure hr"
  plot=gsn_csm_xy(wks,time,Clitterhrglobal(:,1,:),res)
  res@tiMainString="global litter passive hr"
  plot=gsn_csm_xy(wks,time,Clitterhrglobal(:,2,:),res)
  res@tiMainString="global soil fast pool hr"
  plot=gsn_csm_xy(wks,time,Csoilhrglobal(:,0,:),res)
  res@tiMainString="global soil slow pool hr"
  plot=gsn_csm_xy(wks,time,Csoilhrglobal(:,1,:),res)
  res@tiMainString="global soil passive pool hr"
  plot=gsn_csm_xy(wks,time,Csoilhrglobal(:,2,:),res)
  
  res@tiMainString="global C hr"
  plot=gsn_csm_xy(wks,time,Ctothrglobal,res)
  res@trYMaxF = 120
  res@trYMinF = 20
  plot3=gsn_csm_xy(wks,time,HRglobal(0,:,:),res)
  delete(res@trYMaxF)
  delete(res@trYMinF)

  res@tiMainString="global C weighted turnover time(CWT)"
  plot=gsn_csm_xy(wks,time,Ctotcwtglobal,res)

  res@tiMainString="global CWT vs MA"
  plot=gsn_csm_xy(wks,time,onedtond(ndtooned((/Ctotcwtglobal,Ctotageglobal/)),(/nsim*2,nyr/)),res)
 
  plot2=gsn_csm_xy(wks,time,NPPglobal(0,:,:),res)

  gsn_panel(wks,(/plot0,plot1,plot2,plot3/),(/2,2/),False)
end
