; The goals of this script are the following:
; - adapt the original script mkinitialmatrix.ncl to the paths where our cable
;   output is stored:
;   The original script produces a file 'outAC.nc' which is referenced in 
;   the computation of the fluxes 'GetFluxFromCABLEoutput.txt'
; - Generalize the script to compute the result for different years to compare:
;   Where the original script produces the one file 'outAC.nc' using the first 
;   year of output (1901) we now produce a file for every year from 1901 to 2004 
begin
DirPathStr="/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4/"

npatch=10
nland=5656
npool=9

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

;The following loop over the years should actually not be necessarry since the results should all be equal, 
;It is however a test we should do
; alternatively we set year = 1901
do year=1901,2004 
	print(year)
	;==== remove outfile if it exists
	outFilePath=DirPathStr+"mm_"+year+"_outAC.nc"
	command="if [ -f "+outFilePath+" ];then rm "+outFilePath+";fi"
	;print(command)
	system(command)
	
	fin = addfile(DirPathStr+"/out_ncar_"+year+"_ndep.nc","r")
	fout = addfile (outFilePath, "c")  ; open output file
	
	A=new((/npool,npool,npatch,nland/),float)
	C=new((/npool,npool,npatch,nland/),float)
	A=0
	do ipool=0,npool-1
	      A(ipool,ipool,:,:)=-1
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
	
	
	C=0
	iveg=where(ismissing(fin->iveg),18,fin->iveg) 
	do ipool=0,npool-1
		C(ipool,ipool,:,:)=onedtond(\
	      		 where(\
				tau(ndtooned(iveg-1),ipool) .eq. -9999,\
				C@_FillValue,\
				1.0/tau(ndtooned(iveg-1),ipool))/365.0 \
	        	 *where( \
			 	xkoptlitter(ndtooned(iveg-1),ipool) .eq. -9999, \
				C@_FillValue, \ 
				xkoptlitter(ndtooned(iveg-1),ipool) \
			 )  \
	        	 *where( \
			 	xkoptsoil(ndtooned(iveg-1),ipool) .eq. -9999, \
				C@_FillValue,\
				xkoptsoil(ndtooned(iveg-1),ipool)\
			) \
	        	 *where(\
			 	fracLigninplant(ndtooned(iveg-1),ipool) .eq. -9999,\
				C@_FillValue, \
				exp(-3.0*fracLigninplant(ndtooned(iveg-1),ipool))\
			) \
	        	*(1-0.75*(silt(ipool)+clay(ipool))),\
			(/npatch,nland/)\
		) 
	end do
	C@pool_name  = (/"leaf,root,wood,metabolic,structure,CWD,fast,slow,passive"/)
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
