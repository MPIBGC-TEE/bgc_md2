load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
startyr=1901
endyr=2100
nyr=endyr-startyr+1
nsim=9
FileC="CecosystemC.txt"
FileCN="CecosystemCN.txt"
FileNPPC="NPPCtotal.txt"
FileHRC="HeteroRespCtotal.txt"
FileNPPCN="NPPCNtotal.txt"
FileHRCN="HeteroRespCNtotal.txt"
dataC=asciiread(FileC,(/nsim,nyr/),"float")
dataCN=asciiread(FileCN,(/nsim,nyr/),"float")
dataNPPC=asciiread(FileNPPC,(/nsim/3,nyr/),"float")
dataHRC=asciiread(FileHRC,(/nsim/3,nyr/),"float")
dataNPPCN=asciiread(FileNPPCN,(/nsim/3,nyr/),"float")
dataHRCN=asciiread(FileHRCN,(/nsim/3,nyr/),"float")

wks = gsn_open_wks("pdf","Figure_2")
gsn_define_colormap(wks,"rainbow+white+gray")    ; Change color map
  res                          = True

  res@tiMainFontHeightF        = 0.025
  res@tiXAxisString            = "simulation years"
  res@xyLineColors             = (/"black","black","black","blue","red","green","blue","red","green"/)
  res@xyDashPatterns           = (/0,2,16,16,16,16,2,2,2/)
  res@xyLineThicknesses        = (/4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0/)
  res@trXMinF                  = startyr
  res@trXMaxF                  = endyr
  res@trYMinF                  = -1000
  res@trYMaxF                  = 5000
  res@tmLabelAutoStride                = True

  res@gsnMaximize      = True
  res@vpWidthF         = 0.8            ; Change the aspect ratio, but 
  res@vpHeightF        = 0.6            ; make plot as large as possible.

  time=ispan(startyr,endyr,1)

  lgres                    = True
  lgres@lgLineColors       = (/"blue","red","green","blue","red","green","blue","red","green"/)
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .050           ; legend label font thickness
  lgres@vpWidthF           = 0.225          ; width of legend (NDC)
  lgres@vpHeightF          = 0.225          ; height of legend (NDC)
  lgres@lbPerimOn          = False
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoDashIndex    = False
  lgres@lgDashIndexes      = (/0,0,0,16,16,16,2,2,2/)
  lgLineThicknesses        = (/32.0,32.0,32.0,32.0,32.0,32.0,32.0,32.0,32.0/)

  labels = (/"C S1 X","C S2 X","C S3 X","C S1 delta","C S2 delta","C S3 delta","C S1 E","C S2 E","C S3 E"/)
  nitems = dimsizes(labels)                 ; number of legend items

; Create legend
  lbid = gsn_create_legend(wks,nitems,labels,lgres)

;  plotC=gsn_csm_xy(wks,time,dataC((/0,1,3,4,6,7/),:),res)
  Cpoolstart=dim_avg((dataC(0,0:9)+dataCN(0,0:9)))/2

  res@tmYUseLeft = False
  res@tmYRMode = "Explicit"
  res@tmYRValues=(/Cpoolstart-3000,Cpoolstart-2000,Cpoolstart-1000,Cpoolstart,Cpoolstart+1000,Cpoolstart+2000,Cpoolstart+3000/)
  print(res@tmYRValues)
  res@tmYRLabelsOn = True
  res@tmYRLabelFontHeightF = 0.023
  res@tmYRLabels=(/"-3000","-2000","-1000","0","1000","2000","3000"/)
  res@tmYRMinorValues = fspan(Cpoolstart-3000,Cpoolstart+3000,31)
  print(res@tmYRMinorValues)
  res@trYMinF                  = -1000
  res@trYMaxF                  = 4000

  plotC=gsn_csm_xy(wks,time,dataC,res)
  tmp=dataC
  tmp((/1,2,7,8/),:)=tmp((/1,2,7,8/),:)-conform_dims((/4,nyr/),dim_avg(tmp((/1,2,7,8/),0:9)),0)+Cpoolstart
  tmp((/0,3,6/),:)=tmp@_FillValue
  plot1=gsn_csm_xy(wks,time,tmp,res)

  res@tmYUseLeft = True
  amres = True
  amres@amParallelPosF   = -0.31           ; move legend to the right
  amres@amOrthogonalPosF = -0.322            ; move the legend down
;  annoid2 = gsn_add_annotation(plotC,lbid,amres)

  plotCN=gsn_csm_xy(wks,time,dataCN,res)

 gsn_panel(wks,(/plotC,plotCN/),(/2,1/),False)

  res@trYMinF                  = 0
  res@trYMaxF                  = 3000

  plot=gsn_csm_xy(wks,time,reshape((/dataCN(0:2,:),dataCN(0:2,:)/),(/6,nyr/)),res)
  
  print(Cpoolstart)

  res@trYMinF                  = -1000
  res@trYMaxF                  = 4000
  tmp1=dataC
  tmp1((/0,6/),:)=tmp1((/2,8/),:)-conform_dims((/2,nyr/),dim_avg(tmp1((/2,8/),0:9)),0)+Cpoolstart
  tmp1((/1,7/),:)=dataCN((/2,8/),:)-conform_dims((/2,nyr/),dim_avg(dataCN((/2,8/),0:9)),0)+Cpoolstart
  tmp1(3,:)=tmp1(5,:)
  tmp1(4,:)=dataCN(5,:)
  tmp1((/2,5,8/),:)=tmp1@_FillValue
  res@tmYUseLeft = False
  plot=gsn_csm_xy(wks,time,tmp1,res)
  
  gsn_panel(wks,(/plot1,plot/),(/2,1/),False)

  res@trYMinF                  = -600
  res@trYMaxF                  = 1000
  tmp1=dataC
  tmp1((/0,6/),:)=tmp1((/0,6/),:)-conform_dims((/2,nyr/),dim_avg(tmp1((/0,6/),0:9)),0);+Cpoolstart
  tmp1((/1,7/),:)=dataCN((/0,6/),:)-conform_dims((/2,nyr/),dim_avg(dataCN((/0,6/),0:9)),0);+Cpoolstart
  tmp1(3,:)=tmp1(3,:)
  tmp1(4,:)=dataCN(3,:)
  tmp1((/2,5,8/),:)=tmp1@_FillValue
  res@tmYUseLeft = True
  plot=gsn_csm_xy(wks,time,tmp1,res)
  
  gsn_panel(wks,(/plot1,plot/),(/2,1/),False)

  res@trYMinF                  = 30
  res@trYMaxF                  = 110

  tmp1(0:2,:)=dataNPPCN
  plot1=gsn_csm_xy(wks,time,tmp1(0:2,:),res)

  tmp1(0:2,:)=dataHRCN
  plot=gsn_csm_xy(wks,time,tmp1(0:2,:),res)


  gsn_panel(wks,(/plot1,plot/),(/1,2/),False)

  res@trYMinF  = 0
  res@trYMaxF  = 3000
  tmp1=dataC

  plot11=gsn_csm_xy(wks,time,tmp1(0:2,:),res)

  res@trYMinF  = -600
  res@trYMaxF  = 1000
  tmp1=dataC - conform(dataC,dataC(:,0),0)
  plot21=gsn_csm_xy(wks,time,tmp1((/0,3,6/),:),res)

  res@trYMinF  = 0
  res@trYMaxF  = 3000
  tmp1=dataCN

  plot12=gsn_csm_xy(wks,time,tmp1(0:2,:),res)

  res@trYMinF  = -600
  res@trYMaxF  = 1000
  tmp1=dataCN - conform(dataCN,dataCN(:,0),0)
  plot22=gsn_csm_xy(wks,time,tmp1((/0,3,6/),:),res)

  gsn_panel(wks,(/plot11,plot12,plot21,plot22/),(/2,2/),False)

  

end
