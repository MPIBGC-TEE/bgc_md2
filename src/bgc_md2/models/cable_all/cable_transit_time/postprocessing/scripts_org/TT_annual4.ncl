load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
print("start ncl")
startyr=(/1901,1901,1901,1901,1901,1901/)
endyr=(/2100,2100,2100,2100,2100,2100/)
;endyr=(/2070,2100,2100,2100,2100,2100/)
;Active=(/True,True,True,True,True,True/)
;Active=(/True,False,False,False,False,False/)
;Active=(/False,True,False,False,False,False/)
;Active=(/False,False,True,False,False,False/)
Active=(/False,False,False,True,False,False/)
;Active=(/False,False,False,False,True,False/)
;Active=(/False,False,False,False,False,True/)
npatch=10
nland=5656
npool=9
patchstart=0
patchend=0
landstart=2000
landend=2000
ntime=365
FilePath="/datastore/lu02b/CABLE-traceability/"
SimName=(/"C/S1/","C/S2/","C/S3/","CN/S1/","CN/S2/","CN/S3/"/)
nsim=dimsizes(SimName)

Xss=new((/npool,npatch,nland/),double)
A=new((/npool,npool,npatch,nland/),double)
A41=new((/ntime,npatch,nland/),double)
A51=new((/ntime,npatch,nland/),double)
A42=new((/ntime,npatch,nland/),double)
A52=new((/ntime,npatch,nland/),double)
A63=new((/ntime,npatch,nland/),double)
NPP=new((/ntime,npatch,nland/),double)
fracCalloc=new((/ntime,3,npatch,nland/),double)
xktemp=new((/ntime,npatch,nland/),double)
xkwater=new((/ntime,npatch,nland/),double)
xkNlimiting=new((/ntime,npatch,nland/),double)
kplant=new((/ntime,3,npatch,nland/),double)
C=new((/npool,npatch,nland/),double)
B=new((/ntime,npool/),double)
Xit=new((/npool,npool/),double)
Xiw=new((/npool,npool/),double)
Xi=new((/npool,npool/),double)
CplantTT=new((/3,npatch,nland/),double)
ClitterTT=new((/3,npatch,nland/),double)
CsoilTT=new((/3,npatch,nland/),double)
;tau=(/3.76,25,10,0.04,0.23,0.824,0.137,5,222/)
;xkoptlitter=0.4
;xkoptsoil=(/0.33,0.6,0.15,0.6,0.16,0.4,0.3,0.2,0.2,0.25,1,0.65,0.5,2,0.5,0.5,0.5/)

;crop:(4,2107)
;initialisation

do isim=0,nsim-1
   print((/SimName(isim)/))
   if(Active(isim))then
      fAC=addfile(FilePath+SimName(isim)+"outAC.nc","r")
      A=flt2dble(fAC->A)
      do ipool=0,npool-1
         C(ipool,:,:)=fAC->C(ipool,ipool,:,:)
         ;print(any(C(ipool,:,:) .eq. 0))
      end do
         ;print(C(3,4,2107))
      do iyr=startyr(isim),endyr(isim)
;         if(iyr .eq. startyr(isim))then
;            if(isim .le. 2)then
;               finage=addfile(FilePath+"/C/S2/annual_TT/out_ncar_0_age.nc","r")
;            else
;               finage=addfile(FilePath+"/CN/S2/annual_TT/out_ncar_0_age.nc","r")
;            end if
;         end if   
;         Cplantage=finage->Cplant
;         Clitterage=finage->Clitter
;         Csoilage=finage->Csoil
               
         print((/iyr/))
         if(iyr .eq. startyr(isim))then
            fin=addfile(FilePath+SimName(isim)+"output/out_ncar_"+iyr+"_ndep.nc","r")
            Cplantnext=flt2dble(fin->Cplant(0,:,:,:))
            Clitternext=Cplantnext
            Csoilnext=Cplantnext
         else
            fin=finnext
         end if
         finnext=addfile(FilePath+SimName(isim)+"output/out_ncar_"+(iyr+1)+"_ndep.nc","r")

         Cplant=where(flt2dble(fin->Cplant) .eq. 0, 1.e-5, flt2dble(fin->Cplant))
         Cplantnext=where(flt2dble(finnext->Cplant(0,:,:,:)) .eq. 0, 1.e-5, flt2dble(finnext->Cplant(0,:,:,:)))
         Clitter=where(flt2dble(fin->Clitter) .eq. 0, 1.e-5, flt2dble(fin->Clitter))
         Clitternext=(/where(flt2dble(finnext->Clitter(0,:,:,:)) .eq. 0,1.e-5,flt2dble(finnext->Clitter(0,:,:,:)))/)
         Csoil=where(flt2dble(fin->Csoil) .eq. 0,1.e-5,flt2dble(fin->Csoil))
         Csoilnext=(/where(flt2dble(finnext->Csoil(0,:,:,:)) .eq. 0,1.e-5,flt2dble(finnext->Csoil(0,:,:,:)))/)

         Cplant=where(flt2dble(fin->Cplant) .eq. 0, 1.e-5, flt2dble(fin->Cplant))
         Clitter=where(flt2dble(fin->Clitter) .eq. 0, 1.e-5, flt2dble(fin->Clitter))
         Csoil=where(flt2dble(fin->Csoil) .eq. 0,1.e-5,flt2dble(fin->Csoil))

         kplant(0:363,:,:,:)=fin->kplant(1:364,:,:,:)
         kplant(364,:,:,:)=finnext->kplant(0,:,:,:)
         kplant=where(abs(kplant) .le. 1.e-5,1.e-5,kplant)

         xktemp(0:363,:,:)=fin->xktemp(1:364,:,:)
         xktemp(364,:,:)=finnext->xktemp(0,:,:)
         xktemp=where(xktemp .eq. 0,1.e-5,xktemp)
         xkwater(0:363,:,:)=fin->xkwater(1:364,:,:)
         xkwater(364,:,:)=finnext->xkwater(0,:,:)
         xkwater=where(xkwater .eq. 0,1.e-5,xkwater)
         xkNlimiting(0:363,:,:)=fin->xkNlimiting(1:364,:,:)
         xkNlimiting(364,:,:)=finnext->xkNlimiting(0,:,:)
         xkNlimiting=where(xkNlimiting .eq. 0,1.e-5,xkNlimiting)
         if(isim .le. 2)then
            xkNlimiting=1.0
         end if

         A74=A(6,3,:,:)
         A75=A(6,4,:,:)
         A76=A(6,5,:,:)
         A85=A(7,4,:,:)
         A86=A(7,5,:,:)
         A87=A(7,6,:,:)
         A97=A(8,6,:,:)
         A98=A(8,7,:,:)
   
         CplantTT=0
         klitterannual=dim_sum_n(conform(Clitter,C(3:5,:,:),(/1,2,3/))*conform(Clitter,xktemp*xkwater*xkNlimiting,(/0,2,3/))*Clitter,0)
         ClitterTT(0,:,:)=klitterannual(0,:,:)*(1-A74)
         ClitterTT(1,:,:)=klitterannual(1,:,:)*(1-A75-A85)
         ClitterTT(2,:,:)=klitterannual(2,:,:)*(1-A76-A86)

         ksoilannual=dim_sum_n(conform(Csoil,C(6:8,:,:),(/1,2,3/))*conform(Csoil,xktemp*xkwater,(/0,2,3/))*Csoil,0)
         CsoilTT(0,:,:)=ksoilannual(0,:,:)*(1-A87-A97)
         CsoilTT(1,:,:)=ksoilannual(1,:,:)*(1-A98)
         CsoilTT(2,:,:)=ksoilannual(2,:,:)

         system("cp -p out_template_annual.nc "+FilePath+SimName(isim)+"annual_TT/out_ncar_"+iyr+"_annualTT.nc")
         fout1=addfile(FilePath+SimName(isim)+"annual_TT/out_ncar_"+iyr+"_annualTT.nc","w")
         fout1->Cplant=0
         fout1->Clitter=(/dble2flt(ClitterTT)/)
         fout1->Csoil=(/dble2flt(CsoilTT)/)
      end do
   end if
end do

end
