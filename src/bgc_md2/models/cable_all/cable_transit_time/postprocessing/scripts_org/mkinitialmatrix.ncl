load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
FilePath="/datastore/lu02b/CABLE-traceability/"
SimName=(/"C/S1","C/S2","C/S3","CN/S1","CN/S2","CN/S3"/)
nsim=dimsizes(SimName)
year=1901
tau=(/(/3.72837    ,10,20,0.04,0.23,0.824,0.137,5,222/),\
      (/1.65467    ,10,30,0.04,0.23,0.824,0.137,5,222/),\
      (/0.52343    ,10,20,0.04,0.23,0.824,0.137,5,222/),\
      (/0.50679    ,10,10,0.04,0.23,0.824,0.137,5,222/),\
      (/1.44000    , 2, 4,0.04,0.23,0.824,0.137,5,222/),\
      (/0.2910     , 0.28918, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/0.21420    , 0.21404, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/0.54065    , 0.54030, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/0.28935    , 0.28935, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/       0.37, 0.37000, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/          1, 1, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/          1, 1, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/          1, 1, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/ 0.43293   , 2, 5,0.04,0.23,0.824,0.137,5,222/),\
      (/          1, 1, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/          1, 1, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/          1, 1, 1,0.04,0.23,0.824,0.137,5,222/),\
      (/-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999/)/)
    
xkoptlitter=(/(/1,1,1,0.4,0.4,0.4,1,1,1/),\
              (/1,1,1,0.4,0.4,0.4,1,1,1/),\
              (/1,1,1,0.3,0.3,0.3,1,1,1/),\
              (/1,1,1,0.6,0.6,0.6,1,1,1/),\
              (/1,1,1,0.3,0.3,0.3,1,1,1/),\
              (/1,1,1,0.3,0.3,0.3,1,1,1/),\
              (/1,1,1,0.3,0.3,0.3,1,1,1/),\
              (/1,1,1,0.2,0.2,0.2,1,1,1/),\
              (/1,1,1,0.2,0.2,0.2,1,1,1/),\
              (/1,1,1,0.4,0.4,0.4,1,1,1/),\
              (/1,1,1,0.4,0.4,0.4,1,1,1/),\
              (/1,1,1,0.4,0.4,0.4,1,1,1/),\
              (/1,1,1,0.4,0.4,0.4,1,1,1/),\
              (/1,1,1,2.0,2.0,2.0,1,1,1/),\
              (/1,1,1,0.4,0.4,0.4,1,1,1/),\
              (/1,1,1,0.4,0.4,0.4,1,1,1/),\
              (/1,1,1,0.4,0.4,0.4,1,1,1/),\
              (/-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999/)/)

xkoptsoil=(/(/1,1,1,1,1,1,0.40,0.40,0.40/),\
            (/1,1,1,1,1,1,0.40,0.40,0.40/),\
            (/1,1,1,1,1,1,0.30,0.30,0.30/),\
            (/1,1,1,1,1,1,0.60,0.60,0.60/),\
            (/1,1,1,1,1,1,0.30,0.30,0.30/),\
            (/1,1,1,1,1,1, 0.3, 0.3, 0.3/),\
            (/1,1,1,1,1,1, 0.3, 0.3, 0.3/),\
            (/1,1,1,1,1,1, 0.2, 0.2, 0.2/),\
            (/1,1,1,1,1,1,0.25, 0.3, 0.3/),\  ;crop *1.25;1.5;1.5 of original number
            (/1,1,1,1,1,1,0.25,0.25,0.25/),\
            (/1,1,1,1,1,1,   1,   1,   1/),\
            (/1,1,1,1,1,1,0.65,0.65,0.65/),\
            (/1,1,1,1,1,1, 0.5, 0.5, 0.5/),\
            (/1,1,1,1,1,1,   2,   2,   2/),\
            (/1,1,1,1,1,1, 0.5, 0.5, 0.5/),\
            (/1,1,1,1,1,1, 1.0, 1.0, 1.0/),\
            (/1,1,1,1,1,1, 1.0, 1.0, 1.0/),\
            (/-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999/)/)

fracLigninplant=(/(/0,0,0,0,0.25,0,0,0,0/),\
                  (/0,0,0,0,0.2,0,0,0,0/),\
                  (/0,0,0,0,0.2,0,0,0,0/),\
                  (/0,0,0,0,0.2,0,0,0,0/),\
                  (/0,0,0,0,0.2,0,0,0,0/),\
                  (/0,0,0,0,0.1,0,0,0,0/),\
                  (/0,0,0,0,0.1,0,0,0,0/),\
                  (/0,0,0,0,0.1,0,0,0,0/),\
                  (/0,0,0,0,0.1,0,0,0,0/),\
                  (/0,0,0,0,0.1,0,0,0,0/),\
                  (/0,0,0,0,0.15,0,0,0,0/),\
                  (/0,0,0,0,0.15,0,0,0,0/),\
                  (/0,0,0,0,0.15,0,0,0,0/),\
                  (/0,0,0,0,0.15,0,0,0,0/),\
                  (/0,0,0,0,0.15,0,0,0,0/),\
                  (/0,0,0,0,0.25,0,0,0,0/),\
                  (/0,0,0,0,0.1,0,0,0,0/),\
                  (/-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999/)/)
silt=(/0,0,0,0,0,0,0.33,0,0/)
clay=(/0,0,0,0,0,0,0.3,0,0/)
npatch=10
nland=5656
npool=9

A=new((/npool,npool,npatch,nland/),float)
C=new((/npool,npool,npatch,nland/),float)

A=0
C=0
do isim=0,nsim-1
   print((/SimName(isim)/))
   fin=addfile(FilePath+SimName(isim)+"/output/out_ncar_"+year+"_ndep.nc","r")
   iveg=where(ismissing(fin->iveg),18,fin->iveg)
;   npatchveg=dim_num_n(.not. ismissing(fin->iveg),0)
   do ipool=0,npool-1
      print((1-0.75*(silt(ipool)+clay(ipool))))
      print(any(where(tau(ndtooned(iveg-1),ipool) .eq. -9999,C@_FillValue,1.0/tau(ndtooned(iveg-1),ipool))/365.0 .eq. 0))
      print(any(where(xkoptlitter(ndtooned(iveg-1),ipool) .eq. -9999,C@_FillValue,xkoptlitter(ndtooned(iveg-1),ipool)) .eq. 0))
      print(any(where(xkoptsoil(ndtooned(iveg-1),ipool) .eq. -9999,C@_FillValue,xkoptsoil(ndtooned(iveg-1),ipool)) .eq. 0))
      print(any(where(fracLigninplant(ndtooned(iveg-1),ipool) .eq. -9999,C@_FillValue,exp(-3.0*fracLigninplant(ndtooned(iveg-1),ipool))).eq. 0))
      ivegoned=(ndtooned(iveg))
;      print(ivegoned(ind0))
;      print(
      tmp=exp(-3.0*where(fracLigninplant(ndtooned(iveg-1),ipool) .eq. -9999,C@_FillValue,fracLigninplant(ndtooned(iveg-1),ipool)))
      C(ipool,ipool,:,:)=onedtond(where(tau(ndtooned(iveg-1),ipool) .eq. -9999,C@_FillValue,1.0/tau(ndtooned(iveg-1),ipool))/365.0 \
                         *where(xkoptlitter(ndtooned(iveg-1),ipool) .eq. -9999,C@_FillValue,xkoptlitter(ndtooned(iveg-1),ipool))  \
                         *where(xkoptsoil(ndtooned(iveg-1),ipool) .eq. -9999,C@_FillValue,xkoptsoil(ndtooned(iveg-1),ipool)) \
                         *where(fracLigninplant(ndtooned(iveg-1),ipool) .eq. -9999,C@_FillValue,exp(-3.0*fracLigninplant(ndtooned(iveg-1),ipool))) \
                         *(1-0.75*(silt(ipool)+clay(ipool))),(/npatch,nland/)) 
      A(ipool,ipool,:,:)=-1
      print((/ipool/))
      print(any(C(ipool,ipool,:,:) .eq. 0))
;      print(C(ipool,ipool,:,:))
   end do
   A(3:5,0,:,:)= (/fin->fromLeaftoL (1,:,:,:)/)
   A(3:5,1,:,:)= (/fin->fromRoottoL (1,:,:,:)/)
   A(3:5,2,:,:)= (/fin->fromWoodtoL (1,:,:,:)/)
   A(6:8,3,:,:)= (/fin->fromMettoS  (1,:,:,:)/)
   A(6:8,4,:,:)= (/fin->fromStrtoS  (1,:,:,:)/)
   A(6:8,5,:,:)= (/fin->fromCWDtoS  (1,:,:,:)/)
   A(7  ,6,:,:)= (/fin->fromSOMtoSOM(1,0,:,:)/)
   A(8  ,6,:,:)= (/fin->fromSOMtoSOM(1,1,:,:)/)
   A(8  ,7,:,:)= (/fin->fromSOMtoSOM(1,2,:,:)/)
   A@pool_name  = (/"leaf,root,wood,metabolic,structure,CWD,fast,slow,passive"/)
   C@pool_name  = (/"leaf,root,wood,metabolic,structure,CWD,fast,slow,passive"/)
   
   system("if [ -f "+FilePath+SimName(isim)+"/outAC.nc ];then rm "+FilePath+SimName(isim)+"/outAC.nc;fi")
   fout  = addfile (FilePath+SimName(isim)+"/outAC.nc", "c")  ; open output file

    ;===================================================================
    ; explicitly declare file definition mode. Improve efficiency.
    ;===================================================================
        setfileoption(fout,"DefineMode",True)

    ;===================================================================
    ; create global attributes of the file
    ;===================================================================
        fAtt               = True
	fAtt@creation_date = systemfunc ("date")        
	fileattdef( fout, fAtt )            ; copy file attributes    

    ;===================================================================
    ; predefine the coordinate variables and their dimensionality
    ; Note: to get an UNLIMITED record dimension, we set the dimensionality
    ; to -1 (or the actual size) and set the dimension name to True.
    ;===================================================================
        dimNames = (/"poolx", "pooly", "patch", "land"/)  
	dimSizes = (/  npool,   npool,  npatch,  nland/) 
	dimUnlim = (/  False,   False,   False,  False/)   
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
       filevardef(fout,   "A"    ,"float"      ,(/"poolx","pooly","patch","land"/))
       filevardef(fout,   "C"    ,"float"      ,(/"poolx","pooly","patch","land"/))
                                                              ; different from name on script
    ;===================================================================
    ; explicitly exit file definition mode. **NOT REQUIRED**
    ;===================================================================
        setfileoption(fout,"DefineMode",False)

    ;===================================================================
    ; output only the data values since the dimensionality and such have
    ; been predefined. The "(/", "/)" syntax tells NCL to only output the
    ; data values to the predefined locations on the file.
    ;====================================================================
       fout->A       = (/A/)
       fout->C       = (/C/)

end do
end
