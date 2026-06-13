load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
print("start ncl")
startyr=(/1901,1985,1901,1901,1901,1901/)
endyr=(/2100,2092,2100,2100,2100,2100/)
Active=(/False,False,False,False,True,False/)
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

A=new((/npool,npool,npatch,nland/),float)
A41=new((/ntime,npatch,nland/),float)
A51=new((/ntime,npatch,nland/),float)
A42=new((/ntime,npatch,nland/),float)
A52=new((/ntime,npatch,nland/),float)
A63=new((/ntime,npatch,nland/),float)
NPP=new((/ntime,npatch,nland/),float)
fracCalloc=new((/ntime,3,npatch,nland/),float)
xktemp=new((/ntime,npatch,nland/),float)
xkwater=new((/ntime,npatch,nland/),float)
xkNlimiting=new((/ntime,npatch,nland/),float)
kplant=new((/ntime,3,npatch,nland/),float)
C=new((/npool,npatch,nland/),float)
B=new((/ntime,npool/),float)
Xit=new((/npool,npool/),float)
Xiw=new((/npool,npool/),float)
Xi=new((/npool,npool/),float)
Cplantannual=new((/3,npatch,nland/),float)
Clitterannual=new((/3,npatch,nland/),float)
Csoilannual=new((/3,npatch,nland/),float)
Cplantage0=new((/3,npatch,nland/),float)
Clitterage0=new((/3,npatch,nland/),float)
Csoilage0=new((/3,npatch,nland/),float)
Cplantage=new((/365,3,npatch,nland/),float)
Clitterage=new((/365,3,npatch,nland/),float)
Csoilage=new((/365,3,npatch,nland/),float)

do isim=0,nsim-1
   print((/SimName(isim)/))
   if(Active(isim))then
      fAC=addfile(FilePath+SimName(isim)+"outAC.nc","r")
      A=fAC->A
      do ipool=0,npool-1
         C(ipool,:,:)=fAC->C(ipool,ipool,:,:)
      end do
      do iyr=startyr(isim),endyr(isim)
         if(iyr .eq. startyr(isim) .and. startyr(isim) .lt. 2000)then
            system("dmget "+FilePath+SimName(isim)+"output/out_ncar_1???_ndep.nc&")
         else
            if(iyr .eq. 2000 .or. iyr .eq. startyr(isim) .and. startyr(isim) .gt. 2000)then
               system("dmput -r -K 0:40959 "+FilePath+SimName(isim)+"output/out_ncar_1???_ndep.nc&")
               system("dmget "+FilePath+SimName(isim)+"output/out_ncar_2???_ndep.nc&")
            end if
         end if   
               
         print((/iyr/))
         if(iyr .eq. startyr(isim))then
            fin=addfile(FilePath+SimName(isim)+"output/out_ncar_"+iyr+"_ndep.nc","r")
            Cplantnext=fin->Cplant(0,:,:,:)
            Clitternext=Cplantnext
            Csoilnext=Cplantnext
         else
            fin=finnext
         end if
         finnext=addfile(FilePath+SimName(isim)+"output/out_ncar_"+(iyr+1)+"_ndep.nc","r")

         if(iyr .eq. 1901)then
            finage=addfile(FilePath+SimName(isim)+"annual_TT/out_ncar_0_age.nc","r")
         else
            finage=addfile(FilePath+SimName(isim)+"annual_TT/out_ncar_"+(iyr-1)+"_0_age.nc","r")
         end if

         Cplantage0=finage->Cplant
         Clitterage0=finage->Clitter
         Csoilage0=finage->Csoil

         Cplant=where(fin->Cplant .eq. 0, 1.e-5, fin->Cplant)
         Cplantnext=where(finnext->Cplant(0,:,:,:) .eq. 0, 1.e-5, finnext->Cplant(0,:,:,:))
         Clitter=where(fin->Clitter .eq. 0, 1.e-5, fin->Clitter)
         Clitternext=(/where(finnext->Clitter(0,:,:,:) .eq. 0,1.e-5,finnext->Clitter(0,:,:,:))/)
         Csoil=where(fin->Csoil .eq. 0,1.e-5,fin->Csoil)
         Csoilnext=(/where(finnext->Csoil(0,:,:,:) .eq. 0,1.e-5,finnext->Csoil(0,:,:,:))/)

         Cplant(0:363,:,:,:)=Cplant(1:364,:,:,:)
         Cplant(364,:,:,:)=Cplantnext(:,:,:)
         Clitter(0:363,:,:,:)=Clitter(1:364,:,:,:)
         Clitter(364,:,:,:)=Clitternext(:,:,:)
         Csoil(0:363,:,:,:)=Csoil(1:364,:,:,:)
         Csoil(364,:,:,:)=Csoilnext(:,:,:)
         kplant(0:363,:,:,:)=fin->kplant(1:364,:,:,:)
         kplant(364,:,:,:)=finnext->kplant(0,:,:,:)
         kplant=where(abs(kplant) .le. 1.e-5,1.e-5,kplant)

         A41(0:363,:,:)=(/fin->fromLeaftoL(1:364,0,:,:)/)
         A41(364,:,:)=(/finnext->fromLeaftoL(0,0,:,:)/)
         A41C=A41*Cplant(:,0,:,:)*kplant(:,0,:,:)
         A51(0:363,:,:)=(/fin->fromLeaftoL(1:364,1,:,:)/)
         A51(364,:,:)=(/finnext->fromLeaftoL(0,1,:,:)/)
         A51C=A51*Cplant(:,0,:,:)*kplant(:,0,:,:)
         A42(0:363,:,:)=(/fin->fromRoottoL(1:364,0,:,:)/)
         A42(364,:,:)=(/finnext->fromRoottoL(0,0,:,:)/)
         A42C=A42*Cplant(:,2,:,:)*kplant(:,2,:,:)
         A52(0:363,:,:)=(/fin->fromRoottoL(1:364,1,:,:)/)
         A52(364,:,:)=(/finnext->fromRoottoL(0,1,:,:)/)
         A52C=A52*Cplant(:,2,:,:)*kplant(:,2,:,:)
         A63(0:363,:,:)=(/fin->fromWoodtoL(1:364,2,:,:)/)
         A63(364,:,:)=(/finnext->fromWoodtoL(0,2,:,:)/)
         A63C=A63*Cplant(:,1,:,:)*kplant(:,1,:,:)
         NPP(0:363,:,:)=fin->NPP(1:364,:,:)*12.01e-6*3600*24
         NPP(364,:,:)=finnext->NPP(0,:,:)*12.01e-6*3600*24
         NPP=where(NPP .lt. 0,0,NPP)
         fracCalloc(0:363,:,:,:)=fin->fracCalloc(1:364,:,:,:)
         fracCalloc(364,:,:,:)=finnext->fracCalloc(0,:,:,:)
         fracCalloc(:,0,:,:)=where(NPP .lt. 0 .and. Cplant(:,0,:,:) .eq. 0,0,fracCalloc(:,0,:,:))
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

         klitterC=conform(Clitter,C(3:5,:,:),(/1,2,3/))*conform(Clitter,xktemp*xkwater*xkNlimiting,(/0,2,3/))*Clitter
         ksoilC=conform(Csoil,C(6:8,:,:),(/1,2,3/))*conform(Csoil,xktemp*xkwater,(/0,2,3/))*Csoil

         do iday=0,364
               print((/iday/))
            if(iday .eq. 0)then
               Cplantage(iday,0,:,:)=1.0/365+Cplantage0(0,:,:)*(Cplant(iday,0,:,:)-NPP(iday,:,:)*fracCalloc(iday,0,:,:))/Cplant(iday,0,:,:)
               Cplantage(iday,1,:,:)=1.0/365+Cplantage0(1,:,:)*(Cplant(iday,2,:,:)-NPP(iday,:,:)*fracCalloc(iday,2,:,:))/Cplant(iday,2,:,:)
               Cplantage(iday,2,:,:)=1.0/365+Cplantage0(2,:,:)*(Cplant(iday,1,:,:)-NPP(iday,:,:)*fracCalloc(iday,1,:,:))/Cplant(iday,1,:,:)
             
               Clitterage(iday,0,:,:)=1.0/365+Clitterage0(0,:,:)-((Clitterage0(0,:,:)-Cplantage0(0,:,:))*A41C(iday,:,:)+(Clitterage0(0,:,:)-Cplantage0(1,:,:))*A42C(iday,:,:))/Clitter(iday,0,:,:)
               Clitterage(iday,1,:,:)=1.0/365+Clitterage0(1,:,:)-((Clitterage0(1,:,:)-Cplantage0(0,:,:))*A51C(iday,:,:)+(Clitterage0(1,:,:)-Cplantage0(1,:,:))*A52C(iday,:,:))/Clitter(iday,1,:,:)
               Clitterage(iday,2,:,:)=1.0/365+Clitterage0(2,:,:)-((Clitterage0(2,:,:)-Cplantage0(2,:,:))*A63C(iday,:,:))/Clitter(iday,2,:,:)
               Csoilage(iday,0,:,:)=1.0/365+Csoilage0(0,:,:)-((Csoilage0(0,:,:)-Clitterage0(0,:,:))*A74*klitterC(iday,0,:,:)+(Csoilage0(0,:,:)-Clitterage0(1,:,:))*A75*klitterC(iday,1,:,:) \
                                   +(Csoilage0(0,:,:)-Clitterage0(2,:,:))*A76*klitterC(iday,2,:,:))/Csoil(iday,0,:,:)
               Csoilage(iday,1,:,:)=1.0/365+Csoilage0(1,:,:)-((Csoilage0(1,:,:)-Clitterage0(1,:,:))*A85*klitterC(iday,1,:,:)+(Csoilage0(1,:,:)-Clitterage0(2,:,:))*A86*klitterC(iday,2,:,:) \
                                   +(Csoilage0(1,:,:)-Csoilage0(0,:,:))*A87*ksoilC(iday,0,:,:))/Csoil(iday,1,:,:)
               Csoilage(iday,2,:,:)=1.0/365+Csoilage0(2,:,:)-((Csoilage0(2,:,:)-Csoilage0(1,:,:))*A98*ksoilC(iday,1,:,:)\
                               +(Csoilage0(2,:,:)-Csoilage0(0,:,:))*A97*ksoilC(iday,0,:,:))/Csoil(iday,2,:,:)
            else
               Cplantage(iday,0,:,:)=1.0/365+Cplantage(iday-1,0,:,:)*(Cplant(iday,0,:,:)-NPP(iday,:,:)*fracCalloc(iday,0,:,:))/Cplant(iday,0,:,:)
               Cplantage(iday,1,:,:)=1.0/365+Cplantage(iday-1,1,:,:)*(Cplant(iday,2,:,:)-NPP(iday,:,:)*fracCalloc(iday,2,:,:))/Cplant(iday,2,:,:)
               Cplantage(iday,2,:,:)=1.0/365+Cplantage(iday-1,2,:,:)*(Cplant(iday,1,:,:)-NPP(iday,:,:)*fracCalloc(iday,1,:,:))/Cplant(iday,1,:,:)
               Clitterage(iday,0,:,:)=1.0/365+Clitterage(iday-1,0,:,:)-((Clitterage(iday-1,0,:,:)-Cplantage(iday-1,0,:,:))*A41C(iday,:,:) \
                                   +(Clitterage(iday-1,0,:,:)-Cplantage(iday-1,1,:,:))*A42C(iday,:,:))/Clitter(iday,0,:,:)
               Clitterage(iday,1,:,:)=1.0/365+Clitterage(iday-1,1,:,:)-((Clitterage(iday-1,1,:,:)-Cplantage(iday-1,0,:,:))*A51C(iday,:,:) \ 
                                   +(Clitterage(iday-1,1,:,:)-Cplantage(iday-1,1,:,:))*A52C(iday,:,:))/Clitter(iday,1,:,:)
               Clitterage(iday,2,:,:)=1.0/365+Clitterage(iday-1,2,:,:)-((Clitterage(iday-1,2,:,:)-Cplantage(iday-1,2,:,:))*A63C(iday,:,:))/Clitter(iday,2,:,:)
               Csoilage(iday,0,:,:)=1.0/365+Csoilage(iday-1,0,:,:)-((Csoilage(iday-1,0,:,:)-Clitterage(iday-1,0,:,:))*A74*klitterC(iday,0,:,:) \
                                   +(Csoilage(iday-1,0,:,:)-Clitterage(iday-1,1,:,:))*A75*klitterC(iday,1,:,:) \
                                   +(Csoilage(iday-1,0,:,:)-Clitterage(iday-1,2,:,:))*A76*klitterC(iday,2,:,:))/Csoil(iday,0,:,:)
               Csoilage(iday,1,:,:)=1.0/365+Csoilage(iday-1,1,:,:)-((Csoilage(iday-1,1,:,:)-Clitterage(iday-1,1,:,:))*A85*klitterC(iday,1,:,:) \
                                   +(Csoilage(iday-1,1,:,:)-Clitterage(iday-1,2,:,:))*A86*klitterC(iday,2,:,:) \
                                   +(Csoilage(iday-1,1,:,:)-Csoilage(iday-1,0,:,:))*A87*ksoilC(iday,0,:,:))/Csoil(iday,1,:,:)
               Csoilage(iday,2,:,:)=1.0/365+Csoilage(iday-1,2,:,:)-((Csoilage(iday-1,2,:,:)-Csoilage(iday-1,1,:,:))*A98*ksoilC(iday,1,:,:)\
                                +(Csoilage(iday-1,2,:,:)-Csoilage(iday-1,0,:,:))*A97*ksoilC(iday,0,:,:))/Csoil(iday,2,:,:)
            end if
               print((/"Cplant age"/))
               print((/Cplantage(iday,:,0,3000)/))
               print((/"Cplant"/))
               print((/Cplant(iday,:,0,3000)/))
               print((/"NPP fra"/))
               print((/NPP(iday,0,3000)*fracCalloc(iday,:,0,3000)/))
         end do
               
               
         system("cp -p out_template_annual.nc "+FilePath+SimName(isim)+"annual_TT/out_ncar_"+iyr+"_annual_age.nc")
         system("cp -p out_template_annual.nc "+FilePath+SimName(isim)+"annual_TT/out_ncar_"+iyr+"_0_age.nc")
         fout1=addfile(FilePath+SimName(isim)+"annual_TT/out_ncar_"+iyr+"_annual_age.nc","w")
         fout2=addfile(FilePath+SimName(isim)+"annual_TT/out_ncar_"+iyr+"_0_age.nc","w")
       
         fout1->Cplant=(/dim_sum_n(Cplantage*Cplant,0)/dim_sum_n(Cplant,0)/)
         fout1->Clitter=(/dim_sum_n(Clitterage*Clitter,0)/dim_sum_n(Clitter,0)/)
         fout1->Csoil=(/dim_sum_n(Csoilage*Csoil,0)/dim_sum_n(Csoil,0)/)

         fout1->fromLeaftoL(0,:,:)=(/dim_sum_n(A41C,0)/)
         fout1->fromLeaftoL(1,:,:)=(/dim_sum_n(A51C,0)/)
         fout1->fromLeaftoL(2,:,:)=0
         fout1->fromWoodtoL(0,:,:)=0
         fout1->fromWoodtoL(1,:,:)=0
         fout1->fromWoodtoL(2,:,:)=(/dim_sum_n(A63C,0)/)
         fout1->fromRoottoL(0,:,:)=(/dim_sum_n(A42C,0)/)
         fout1->fromRoottoL(1,:,:)=(/dim_sum_n(A52C,0)/)
         fout1->fromRoottoL(2,:,:)=0
         fout1->fromMettoS(0,:,:) =(/dim_sum_n(klitterC(:,0,:,:)*conform(Clitter(:,0,:,:),A74,(/1,2/)),0)/)
         fout1->fromMettoS(1,:,:) =0
         fout1->fromMettoS(2,:,:) =0
         fout1->fromStrtoS(0,:,:) =(/dim_sum_n(klitterC(:,1,:,:)*conform(Clitter(:,1,:,:),A75,(/1,2/)),0)/)
         fout1->fromStrtoS(1,:,:) =(/dim_sum_n(klitterC(:,1,:,:)*conform(Clitter(:,1,:,:),A85,(/1,2/)),0)/)
         fout1->fromStrtoS(2,:,:) =0
         fout1->fromCWDtoS(0,:,:) =(/dim_sum_n(klitterC(:,2,:,:)*conform(Clitter(:,2,:,:),A76,(/1,2/)),0)/)
         fout1->fromCWDtoS(1,:,:) =(/dim_sum_n(klitterC(:,2,:,:)*conform(Clitter(:,2,:,:),A86,(/1,2/)),0)/)
         fout1->fromCWDtoS(2,:,:) =0
         fout1->fromSOMtoSOM(0,:,:) =(/dim_sum_n(ksoilC(:,0,:,:)*conform(Csoil(:,0,:,:),A87,(/1,2/)),0)/)
         fout1->fromSOMtoSOM(1,:,:) =(/dim_sum_n(ksoilC(:,0,:,:)*conform(Csoil(:,0,:,:),A97,(/1,2/)),0)/)
         fout1->fromSOMtoSOM(2,:,:) =(/dim_sum_n(ksoilC(:,1,:,:)*conform(Csoil(:,1,:,:),A98,(/1,2/)),0)/)
         fout1->fracCalloc(0,:,:) =(/dim_sum_n(NPP(:,:,:)*fracCalloc(:,0,:,:),0)/)
         fout1->fracCalloc(1,:,:) =(/dim_sum_n(NPP(:,:,:)*fracCalloc(:,1,:,:),0)/)
         fout1->fracCalloc(2,:,:) =(/dim_sum_n(NPP(:,:,:)*fracCalloc(:,2,:,:),0)/)
      
         fout2->Cplant=(/Cplantage(364,:,:,:)/)
         fout2->Clitter=(/Clitterage(364,:,:,:)/)
         fout2->Csoil=(/Csoilage(364,:,:,:)/)
      end do
   end if
end do

end
