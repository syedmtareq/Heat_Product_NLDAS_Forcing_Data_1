      program NLDAS_T_HI_NDHS_grib

!---- This computer program was developed by: 
!---- Dr. Mohammad Al-Hamdan, PhD
!---- USRA at NASA/MSFC
!---- mohammad.alhamdan@nasa.gov
!---- June 15, 2010 
!---- Modified by Bill Crosson, January 4, 2012

        USE MSIMSL
	    USE MSFLIB					
	    USE PORTLIB

	    parameter (npts = 103936)
!                       = number of grid cells

	    parameter (maxdays = 185)
!                       = maximum number of days for which model is run

	    parameter (num_zones = 4)
!                       = number of time zones in area of interest (here US)

	    character fname*39,yearname*4,filename*100, year_next_name*4
		character*2 dayname, monthssname, day_next_name, month_next_name 
	    character*200 option1,option2,option3,option4,option5,option6
	    character*250 optionss, options3
		character dirin*27       !! Must match string length exactly

		real SpecHum_2D(npts,24+num_zones), T_K_2D(npts,24+num_zones), AtmPres_2D(npts,24+num_zones)
		real uwind_2D(npts,24+num_zones), vwind_2D(npts,24+num_zones), Precip_2D(npts,24+num_zones)

	    real Tmax_F(npts,maxdays),Tmin_F(npts,maxdays),Tmean_F(npts,maxdays),HI_max(npts,maxdays)
		real Time_Tmin(npts,maxdays),Time_Tmax(npts,maxdays),Time_HImax(npts,maxdays) 
	    real RH_time_of_Tmax(npts,maxdays),RH_time_of_Tmin(npts,maxdays) 
	    real WSmax(npts,maxdays),Time_WSmax(npts,maxdays),WD_time_of_WSmax(npts,maxdays) 
		real NDHS(npts,maxdays)
		real Precip_daily(npts,maxdays)
		  
		real t_off(npts)
!            = offset in hours between local standard time and UTC

		integer n,febdays,year,jday,hour,monthss
		integer day_of_year

	    write(*,*) 'Enter Year'	  
	    read (*,*) year					  
	    write(yearname,'(i4.4)') year
		
	    write(*,*) 'Enter first month to process'	  
	    read (*,*) month_first					  
		
	    write(*,*) 'Enter last month to process'	  
	    read (*,*) month_last					  
		
	    febdays = 28
!if(mod(year,4) .eq. 0) febdays = 29
		if(mod(year,4) .eq. 0 .and. (mod(year,100) .ne. 0 .or. mod(year,400) .eq. 0)) febdays = 29
	   
	    day_of_year = 0

		dirin = 'C:\NLDAS_Forcing_Data\    \'										 

		write(dirin(23:26),'(a4)') yearname

		open(66,file='NLDAS_Centroids_ESRI_TimeOffset.dat',status='old')
        read(66,*)

	    do 55 i=1,npts
		   read(66,*,end=200) j,t_off(i)  
		   if (t_off(i) .ge. -4.) t_off(i) = -5.  ! -4 corresponds to Atlantic time, not in CONUS region 
 55		end do

 200	close(66)

	    do 25 monthss = month_first,month_last
	    
	     write(monthssname,'(i2.2)') monthss
		 month_next = monthss
		 iyear_next = year

		 n = 31

		 if (monthss.eq.2)  n=febdays 
		 if (monthss.eq.4)  n=30
		 if (monthss.eq.6)  n=30
		 if (monthss.eq.9)  n=30
		 if (monthss.eq.11) n=30

		 do 35 jday=1,n
		    day_of_year = day_of_year + 1

!           Determine day and month for following day

			if (jday .lt. n) then
			   jday_next = jday + 1
			else
			   if (monthss .lt. 12) then
			      jday_next = 1
			      month_next = monthss + 1
			   else
			      jday_next = 1
			      month_next = 1
				  iyear_next = year + 1
			   end if
			end if
			 
		    write(dayname,'(i2.2)') jday

		    write(day_next_name,'(i2.2)') jday_next
		    write(month_next_name,'(i2.2)') month_next
		    write(year_next_name,'(i4.4)') iyear_next

			write(*,*) 'Month, Day: ',monthss, jday


!           Next line gets directory of all files for current and next day to enable calculation
!           of daily values for local time, since local day crosses into next day in UTC time.

	        optionss='dir '//dirin//'NLDAS_FORA0125_H.A'//yearname//monthssname//dayname//'.*.002.grb  ' //dirin//'NLDAS_FORA0125_H.A'//year_next_name//month_next_name//day_next_name//'.*.002.grb /b >folded.txt'
			write(*,*) optionss
	        iresults=system(optionss)

		    open(1,file='folded.txt',status='unknown')
      		    
			read(1,*) fname   	! skip over first 5 files 
			read(1,*) fname   	! since hour 0 in EST = hour 5 in UTC
			read(1,*) fname   
			read(1,*) fname   
			read(1,*) fname   
		    
			do 11 hour = 1,25 + (num_zones-1)  ! I know, this could be written as 24+num_zones!

! Spanning 25 hours since max/min values can occur at end or beginning time in 24-hour period
! plus 3 more hours to span U.S. time zones.    

			   write(*,*) 'Month, Day, Hour: ',monthss, jday, hour
			   write(*,*)
			   read(1,'(a39)',end=111) fname   !! make sure this character string matches file names

			   option1=''//dirin//fname//' -d 2 -text -o '//dirin//fname//'.out1'	 !Specific Humidity at 2 m (kg/kg)
			   result1=runqq('wgrib',option1)
			
			   option2=''//dirin//fname//' -d 1 -text -o '//dirin//fname//'.out2'	 !Temperature at 2 m (Kelvin)
			   result2=runqq('wgrib',option2)
			
			   option3=''//dirin//fname//' -d 3 -text -o '//dirin//fname//'.out3'	 !Surface Pressure (Pa)
			   result3=runqq('wgrib',option3)
			
			   option4=''//dirin//fname//' -d 4 -text -o '//dirin//fname//'.out4'	 !u-wind component (m/s)
			   result4=runqq('wgrib',option4)
			
			   option5=''//dirin//fname//' -d 5 -text -o '//dirin//fname//'.out5'	 !v-wind component (m/s)
			   result5=runqq('wgrib',option5)
			
			   option6=''//dirin//fname//' -d 10 -text -o '//dirin//fname//'.out6'	 !Precipitation (mm for previous hour)
			   result6=runqq('wgrib',option6)
			
11		    continue

 111	    optionss='dir '//dirin//'NLDAS_FORA0125_H.A'//yearname//monthssname//dayname//'.*.002.grb.out1  ' //dirin//'NLDAS_FORA0125_H.A'//year_next_name//month_next_name//day_next_name//'.*.002.grb.out1  /b >folded_nldas_rawdata_SpecificHumidity.txt'
	        iresults=system(optionss)

	        optionss='dir '//dirin//'NLDAS_FORA0125_H.A'//yearname//monthssname//dayname//'.*.002.grb.out2  ' //dirin//'NLDAS_FORA0125_H.A'//year_next_name//month_next_name//day_next_name//'.*.002.grb.out2  /b >folded_nldas_rawdata_Temperature_Kelv.txt'
	        iresults=system(optionss)

	        optionss='dir '//dirin//'NLDAS_FORA0125_H.A'//yearname//monthssname//dayname//'.*.002.grb.out3  ' //dirin//'NLDAS_FORA0125_H.A'//year_next_name//month_next_name//day_next_name//'.*.002.grb.out3  /b >folded_nldas_rawdata_AtmPressure_Pa.txt'
	        iresults=system(optionss)
			 
	        optionss='dir '//dirin//'NLDAS_FORA0125_H.A'//yearname//monthssname//dayname//'.*.002.grb.out4  ' //dirin//'NLDAS_FORA0125_H.A'//year_next_name//month_next_name//day_next_name//'.*.002.grb.out4  /b >folded_nldas_rawdata_u_wind.txt'
	        iresults=system(optionss)
			 
	        optionss='dir '//dirin//'NLDAS_FORA0125_H.A'//yearname//monthssname//dayname//'.*.002.grb.out5  ' //dirin//'NLDAS_FORA0125_H.A'//year_next_name//month_next_name//day_next_name//'.*.002.grb.out5  /b >folded_nldas_rawdata_v_wind.txt'
	        iresults=system(optionss)
			 
	        optionss='dir '//dirin//'NLDAS_FORA0125_H.A'//yearname//monthssname//dayname//'.*.002.grb.out6  ' //dirin//'NLDAS_FORA0125_H.A'//year_next_name//month_next_name//day_next_name//'.*.002.grb.out6  /b >folded_nldas_rawdata_Precipitation.txt'
	        iresults=system(optionss)
			
		    call arrays_2D(dirin,num_zones,npts,SpecHum_2D,T_K_2D,AtmPres_2D,uwind_2D,vwind_2D,Precip_2D)

			write(*,*) 'calling stats_heat'
		    call stats_heat(num_zones,npts,maxdays,day_of_year,t_off,SpecHum_2D,T_K_2D,AtmPres_2D,Tmax_F,Tmin_F,HI_max,NDHS,Tmean_F,Time_Tmin,Time_Tmax,Time_HImax,RH_time_of_Tmax,RH_time_of_Tmin)

			write(*,*) 'calling stats_wind'
		    call stats_wind(num_zones,npts,maxdays,day_of_year,t_off,uwind_2D,vwind_2D,WSmax,Time_WSmax,WD_time_of_WSmax)

			write(*,*) 'calling stats_precip'
		    call stats_precip(num_zones,npts,maxdays,day_of_year,t_off,Precip_2D,Precip_daily)
 
			options3 = 'del ' //dirin//'*.out1 ' //dirin//'*.out2 '//dirin//'*.out3 '//dirin//'*.out4 '//dirin//'*.out5 '//dirin//'*.out6 '
			write(*,*) 'about to delete files! ',options3
	        result1=system(options3)

	        close(1)

35		 continue
25	    continue
	
	    write(*,*) 'total of ',day_of_year,' days processed'

        write(filename,'(a18,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Max_HI','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1121,file=dirin//filename,status='unknown')

        write(filename,'(a20,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Max_temp','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1122,file=dirin//filename,status='unknown')

        write(filename,'(a20,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Min_temp','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1123,file=dirin//filename,status='unknown')

        write(filename,'(a28,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Time_of_Max_temp','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1124,file=dirin//filename,status='unknown')

        write(filename,'(a28,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Time_of_Min_temp','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1125,file=dirin//filename,status='unknown')

        write(filename,'(a21,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Mean_temp','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1126,file=dirin//filename,status='unknown')

        write(filename,'(a26,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Time_of_Max_HI','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1127,file=dirin//filename,status='unknown')

        write(filename,'(a16,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_NDHS','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1130,file=dirin//filename,status='unknown')

        write(filename,'(a24,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Total_Precip','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1140,file=dirin//filename,status='unknown')

        write(filename,'(a34,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_RH_at_time_of_Max_Temp','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1151,file=dirin//filename,status='unknown')

        write(filename,'(a34,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_RH_at_time_of_Min_Temp','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1152,file=dirin//filename,status='unknown')

        write(filename,'(a26,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Max_Wind_Speed','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1161,file=dirin//filename,status='unknown')

        write(filename,'(a34,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Time_of_Max_Wind_Speed','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1162,file=dirin//filename,status='unknown')

        write(filename,'(a46,a1,i4,a1,i2,a1,i2,a8)') 'Array_NLDAS_Wind_Dir_at_Time_of_Max_Wind_Speed','_',year,'_',month_first,'-',month_last,'_LST.txt'
	    open (1163,file=dirin//filename,status='unknown')

        do j=1,day_of_year
           write(*,*) 'Day of the Year',j
 	       do i=1,npts
 	          if(HI_max(i,j).lt.0.) HI_max(i,j) = -999.	  
	          if(Tmin_F(i,j).ge.1000.) Tmin_F(i,j) = -999.
	          if(Tmax_F(i,j).ge.1000.) Tmax_F(i,j) = -999.
			  if(Tmax_F(i,j).eq.-999.) NDHS(i,j) = -999.
	          if(Precip_daily(i,j) .ge. 1000.) Precip_daily(i,j) = -999.
	       end do
        end do

118	    do j=1,npts
           write(1121,'(i6,1x,366(F6.1,1x))') j, (HI_max(j,k),k=1,day_of_year)	
	       write(1122,'(i6,1x,366(F6.1,1x))') j, (Tmax_F(j,k),k=1,day_of_year)	
	       write(1123,'(i6,1x,366(F6.1,1x))') j, (Tmin_F(j,k),k=1,day_of_year)	
	       write(1124,'(i6,1x,366(F6.1,1x))') j, (Time_Tmax(j,k),k=1,day_of_year)	
	       write(1125,'(i6,1x,366(F6.1,1x))') j, (Time_Tmin(j,k),k=1,day_of_year)	
	       write(1126,'(i6,1x,366(F6.1,1x))') j, (Tmean_F(j,k),k=1,day_of_year)	
	       write(1127,'(i6,1x,366(F6.1,1x))') j, (Time_HImax(j,k),k=1,day_of_year)	
	       write(1130,'(i6,1x,366(F6.1,1x))') j, (NDHS  (j,k),k=1,day_of_year)	
	       write(1140,'(i6,1x,366(F6.1,1x))') j, (Precip_daily(j,k),k=1,day_of_year)	
	       write(1151,'(i6,1x,366(F6.1,1x))') j, (RH_time_of_Tmax(j,k),k=1,day_of_year)	
	       write(1152,'(i6,1x,366(F6.1,1x))') j, (RH_time_of_Tmin(j,k),k=1,day_of_year)	
	       write(1161,'(i6,1x,366(F6.1,1x))') j, (WSmax(j,k),k=1,day_of_year)	
	       write(1162,'(i6,1x,366(F6.1,1x))') j, (Time_WSmax(j,k),k=1,day_of_year)	
	       write(1163,'(i6,1x,366(F6.1,1x))') j, (WD_time_of_WSmax(j,k),k=1,day_of_year)	
	    end do

        close(1121)
        close(1122)
        close(1123)
        close(1124)
        close(1125)
        close(1126)
		close(1127)	
        close(1130)
        close(1140)
        close(1151)
        close(1152)
        close(1161)
        close(1162)
        close(1163)

	    result3=system("del folded*.txt")
		   
		stop
		end
	
! --------------------------------------------------------------------------------------------------------
		subroutine arrays_2D(dirin,num_zones,npts,SpecHum_2D,T_K_2D,AtmPres_2D,uwind_2D,vwind_2D,Precip_2D)

!---- This subroutine was developed by Bill Crosson, September 6, 2011
!---- Last modified January 4, 2012

		character*44 fname1,fname2,fname3,fname4,fname5,fname6

		integer hour

		character dirin*27       !! Must match string length exactly

		real SpecHum_2D(npts,24+num_zones), T_K_2D(npts,24+num_zones), AtmPres_2D(npts,24+num_zones)
		real uwind_2D(npts,24+num_zones), vwind_2D(npts,24+num_zones), Precip_2D(npts,24+num_zones)

		open(1,file='folded_nldas_rawdata_SpecificHumidity.txt',status='unknown')
		open(2,file='folded_nldas_rawdata_Temperature_Kelv.txt',status='unknown')
		open(3,file='folded_nldas_rawdata_AtmPressure_Pa.txt',status='unknown')
		open(4,file='folded_nldas_rawdata_u_wind.txt',status='unknown')
		open(5,file='folded_nldas_rawdata_v_wind.txt',status='unknown')
		open(6,file='folded_nldas_rawdata_Precipitation.txt',status='unknown')

		do 11 hour = 1,24 + num_zones

		   read(1,'(a44)',end=1000) fname1
 1000	   read(2,'(a44)',end=2000) fname2
 2000	   read(3,'(a44)',end=3000) fname3
 3000	   read(4,'(a44)',end=4000) fname4
 4000	   read(5,'(a44)',end=5000) fname5
 5000	   read(6,'(a44)',end=100)  fname6

 100	   open(2275,file=''//dirin//fname1,status='unknown')
		   open(2276,file=''//dirin//fname2,status='unknown')
		   open(2277,file=''//dirin//fname3,status='unknown')
		   open(2278,file=''//dirin//fname4,status='unknown')
		   open(2279,file=''//dirin//fname5,status='unknown')
		   open(2280,file=''//dirin//fname6,status='unknown')

		   read(2275,*,end=111) 	 ! reading header lines
		   read(2276,*,end=111) 
	       read(2277,*,end=111) 
	       read(2278,*,end=111) 
	       read(2279,*,end=111) 
	       read(2280,*,end=111) 

	       do 51 i=1,npts
		      read(2275,*,end=1515) SpecHum_2D(i,hour)
1515	      read(2276,*,end=1525) T_K_2D(i,hour)
1525	      read(2277,*,end=1535) AtmPres_2D(i,hour)		  
1535	      read(2278,*,end=1545) uwind_2D(i,hour)		  
1545	      read(2279,*,end=1555) vwind_2D(i,hour)		  
1555	      read(2280,*,end=51)   Precip_2D(i,hour)		  
 51		   end do

11		end do

111   close(2275)
      close(2276)
      close(2277)
      close(2278)
      close(2279)
      close(2280)

	  close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      close(6)

	  return
      end
		        
! --------------------------------------------------------------------------------------------------------
	  subroutine stats_heat(num_zones,npts,maxdays,numday,t_off,SpecHum_2D,T_K_2D,AtmPres_2D,Tmax_F,Tmin_F,HI_max,NDHS,Tmean_F,Time_Tmin,Time_Tmax,Time_HImax,RH_time_of_Tmax,RH_time_of_Tmin)

!---- This subroutine was developed by: 
!---- Dr. Mohammad Al-Hamdan, PhD
!---- June 15, 2010 
!---- Modified by Bill Crosson, January 4, 2012

!     Computes statistics for temperature variables - Tmax, Tmin, Tmean, time of Tmax and Tmin, RH at time of Tmax and Tmin, HI_max and NDHS (Net Daily Heat Stress)
	
	  real Tmax_F(npts,maxdays),Tmin_F(npts,maxdays),Tmean_F(npts,maxdays),HI_max(npts,maxdays),NDHS(npts,maxdays)   
	  real Time_Tmin(npts,maxdays),Time_Tmax(npts,maxdays),Time_HImax(npts,maxdays) 
	  real RH_time_of_Tmax(npts,maxdays),RH_time_of_Tmin(npts,maxdays) 

	  real RH_2D(npts,24+num_zones)
      real SpecHum_2D(npts,24+num_zones), T_K_2D(npts,24+num_zones), AtmPres_2D(npts,24+num_zones), HI_2D(npts,24+num_zones), TF_2D(npts,24+num_zones)
      real t_off(npts)
!          = offset in hours between local standard time and UTC

	  HI_hot = 90.
	  HI_cool = 75.

	  do ihr = 1,25   ! looping over 25 hours to span hour 0 to hour 24, inclusive, since max/min can occur at these times.
	    do k=1,npts
		   jhr = ihr - num_zones - 1 - t_off(k)	 ! indexing based on time offset; for EST, jhr=1,25; for PST, jhr=4,28
           if(T_K_2D(k,jhr).gt.10**10) T_K_2D(k,jhr)=10**9
           if(AtmPres_2D(k,jhr).gt.10**10) AtmPres_2D(k,jhr)=10**9
		   if(SpecHum_2D(k,jhr).gt.10**10) SpecHum_2D(k,jhr)=10**9
		   TF_2D(k,jhr)=10**9

		   if (T_K_2D(k,jhr) .lt. 1000. .and. AtmPres_2D(k,jhr) .lt. 1e6 .and. SpecHum_2D(k,jhr) .lt. 10.) then
			 RH_2D(k,jhr)=100*(SpecHum_2D(k,jhr)/((0.622*10**(-2937.4/T_K_2D(k,jhr)-4.9283*log10(T_K_2D(k,jhr))+23.547))/(AtmPres_2D(k,jhr)/100.)))   ! Multiplied by a 100 to convert to Percent as needed in the HI equation.
			 TF_2D(k,jhr)=((T_K_2D(k,jhr)-273.15)*1.8)+32.
	         if (TF_2D(k,jhr).ge.80.) then
		        HI_2D(k,jhr)= -42.379+2.04901523*TF_2D(k,jhr)+10.14333127*RH_2D(k,jhr)-0.22475541*TF_2D(k,jhr)*RH_2D(k,jhr)-6.83783/1000*TF_2D(k,jhr)**2-5.481717/100*RH_2D(k,jhr)**2+1.22874/1000*TF_2D(k,jhr)**2*RH_2D(k,jhr)+8.5282/10000*TF_2D(k,jhr)*RH_2D(k,jhr)**2-1.99/1000000*TF_2D(k,jhr)**2*RH_2D(k,jhr)**2
	         else 
		        HI_2D(k,jhr)= -999.
	         end if	
		  else
		     TF_2D(k,jhr) = -999.
			 HI_2D(k,jhr) = -999.
		  end if 

	    end do

	  end do

11	  do 77 j=1,npts
	   Tmax_F(j,numday)=-999.
	   Tmin_F(j,numday)= 999.
	   HI_max(j,numday)=-999.
	   NDHS(j,numday) = 0.0
	   Time_Tmax(j,numday)=-999.
	   Time_HImax(j,numday)=-999.
	   Time_Tmin(j,numday)=-999.
	   RH_time_of_Tmax(j,numday)=-999.
	   RH_time_of_Tmin(j,numday)=-999.

!      Compute Tmax, Tmin, time of Tmax and Tmin, HI_max, time of HI_max, and RH at time of Tmax and Tmin

	   do ihr=1,25 
		  jhr = ihr - num_zones - 1 - t_off(j)	 
		  if(HI_2D(j,jhr).gt.HI_max(j,numday)) then	
		     HI_max(j,numday) = HI_2D(j,jhr)
			 Time_HImax(j,numday) = ihr-1
		  end if

		  if(TF_2D(j,jhr).gt.Tmax_F(j,numday)) then				
		     Tmax_F(j,numday) = TF_2D(j,jhr)
			 Time_Tmax(j,numday) = ihr-1
			 RH_time_of_Tmax(j,numday) = amin1(100.0,RH_2D(j,jhr))	
		  end if

		  if(TF_2D(j,jhr).lt.Tmin_F(j,numday) .and. TF_2D(j,jhr).gt. -900.) then	
		     Tmin_F(j,numday) = TF_2D(j,jhr)					 
			 Time_Tmin(j,numday) = ihr-1
			 RH_time_of_Tmin(j,numday) = amin1(100.0,RH_2D(j,jhr))	
		  end if

       end do

!  For NDHS, only loop over 24 hours since this is an accumulating variable

	   do ihr=1,24 
		  jhr = ihr - num_zones - 1 - t_off(j)	 
		  if(HI_2D(j,jhr) .gt. HI_hot) then	
		     NDHS(j,numday) = NDHS(j,numday) + HI_2D(j,jhr)-HI_hot
		  end if

		  if(TF_2D(j,jhr) .lt. HI_cool) then	
		     NDHS(j,numday) = NDHS(j,numday) + TF_2D(j,jhr)-HI_cool
		  end if

		  NDHS(j,numday) = amax1(0.,NDHS(j,numday)) 

       end do

       NDHS(j,numday) = amax1(0.,NDHS(j,numday)) 

!  Calculate mean temperature over 24 hours

	   Tsum = 0.
	   do ihr=1,24 
		  jhr = ihr - num_zones - 1 - t_off(j)	 
		  Tsum = Tsum + TF_2D(j,jhr)
       end do

       Tmean_F(j,numday) = Tsum/24. 

77	end do

        return
        end

 ! --------------------------------------------------------------------------------------------------------
	  subroutine stats_wind(num_zones,npts,maxdays,numday,t_off,uwind_2D,vwind_2D,WSmax,Time_WSmax,WD_time_of_WSmax)

!---- This subroutine was developed by: 
!---- Bill Crosson, January 4, 2012

!     Computes statistics for wind variables - WSmax, time of WSmax, wind direction at time of WSmax 
	
	  real uwind_2D(npts,24+num_zones), vwind_2D(npts,24+num_zones),WS_2D(npts,24+num_zones),WD_2D(npts,24+num_zones)
	  real WSmax(npts,maxdays),Time_WSmax(npts,maxdays),WD_time_of_WSmax(npts,maxdays)   

      real t_off(npts)
!          = offset in hours between local standard time and UTC

	  do ihr = 1,25   ! looping over 25 hours to span hour 0 to hour 24, inclusive, since max can occur at these times.
	    do k=1,npts
		   jhr = ihr - num_zones - 1 - t_off(k)	 ! indexing based on time offset; for EST, jhr=1,25; for PST, jhr=4,28
		   if(uwind_2D(k,jhr).gt.10**10) uwind_2D(k,jhr)=10**9
		   if(vwind_2D(k,jhr).gt.10**10) vwind_2D(k,jhr)=10**9
		   WS_2D(k,jhr)=10**9
		   WD_2D(k,jhr)=10**9

		   if (uwind_2D(k,jhr) .lt. 1000. .and. vwind_2D(k,jhr) .lt. 1000.) then
			 WS_2D(k,jhr)= sqrt(uwind_2D(k,jhr)**2 + vwind_2D(k,jhr)**2)    
			 WD_2D(k,jhr)= amod(270. - atan2d(vwind_2D(k,jhr),uwind_2D(k,jhr)),360.) 
		   else
		     WS_2D(k,jhr) = -999.
			 WD_2D(k,jhr) = -999.
		   end if 

	    end do
	  end do

11	  do 77 j=1,npts
	   WSmax(j,numday)=-999.
	   Time_WSmax(j,numday)=-999.
	   WD_time_of_WSmax(j,numday)=-999.

!      Compute WSmax, Time of WSmax, and WD at time of WSmax

	   do ihr=1,25 
		  jhr = ihr - num_zones - 1 - t_off(j)	 
		  if(WS_2D(j,jhr).gt.WSmax(j,numday)) then	
		     WSmax(j,numday) = WS_2D(j,jhr)
			 Time_WSmax(j,numday) = ihr-1
			 WD_time_of_WSmax(j,numday) = WD_2D(j,jhr)	
		  end if
       end do

77	end do

        return
        end

! --------------------------------------------------------------------------------------------------------
	  subroutine stats_precip(num_zones,npts,maxdays,numday,t_off,Precip_2D,Precip_daily)

!---- This subroutine was developed by: 
!---- Bill Crosson, January 4, 2012

!     Computes daily total precipitation

	  real Precip_daily(npts,maxdays) 

	  real Precip_2D(npts,24+num_zones)
      real t_off(npts)
!          = offset in hours between local standard time and UTC

	  do ihr = 1,25   ! looping over 25 hours to span hour 0 to hour 24, inclusive, since max/min can occur at these times.
	    do k=1,npts
		   jhr = ihr - num_zones - 1 - t_off(k)	 ! indexing based on time offset; for EST, jhr=1,25; for PST, jhr=4,28
           if(Precip_2D(k,jhr).gt.10**10) Precip_2D(k,jhr)=10**9
	    end do
	  end do

11	  do 77 j=1,npts
	     Precip_daily(j,numday)=0.

!  Only loop over 24 hours since this is an accumulating variable

	   do ihr=1,24 
		  jhr = ihr - num_zones - 1 - t_off(j)	 
		  Precip_daily(j,numday) = Precip_daily(j,numday) + Precip_2D(j,jhr)
       end do

77	end do

        return
        end

