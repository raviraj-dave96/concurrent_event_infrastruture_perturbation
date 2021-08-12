#!--radians


## Model for simulating percolation of rain and snowmelt water through
## the unsaturated zone to a perched groundwater body.
## Subsurface flow of groundwater along the slope 
## Simulation in a catchment
## Interception by vegetation
## Interception by evaporation modelled as constant value
## Loss of water through roots of vegetation
## Increase in surface runoff contributed by vegetation
## Decrease in infiltration by roots
## Increase in shear strength 
## Calculation of slope stability with the infinite slope model
## Copyright Siva Theo Fan @ 2020
## Distributed for educational scientific purposes only
## Any use of this code should be cited properly
## All rights reserved 

binding
## Routing
 RainStations=raingui.map;          	# map with location of rainstations
 Dem=dtm.map;				# DEM with __m resolution
 Ldd=ldd.map;				    	# local drain direction map created from DEM
 #Veg=veg.map;					# Map showing vegetation cover
 #SoilType=soiltype.map;                 # soil map
 RunOff=runoff;                      	# reported stack of maps with
 RunoffTimeSeries=runoff.tss;        	# reported timeseries with runoff
 Outlet=outlets.map; 			# map with 8 runoff sampling locations
 TotArea=totarea.map;			# scalar map with 1 values in the area of interest
 Clone=clone.map;            		# clone map
 #SoilInfiltrationTable=infilcap.tbl; 	# table with infiltr. cap. of soil types
 #IntercTable=interc.tbl;        		# table with loss of rain through interception by vegetation

## Numerical timestep 
   T=1;                      		# timestep in day 
    
##############################################################
## Meteorological data input
 # Rain=100;              			# constant rain input (cm per day)       
 Rain=rainfall2018.tss;    		# rainfall and snowmelt history (cm/day.) 
 # avgtemp=avgtemp.tss;			# daily average air temperature (degree celcius)
 # avgtemp=avgtemp.tss;			# daily average air temperature (degree celcius)
 # DDF=0.2;					# degree-day factor for snowmelt (cm/day)
 # snowdepth=snowdepth.tss;		# measured snowdepth using snowgauge
 Evap=0;                     		# evapotranspiration cm/day switched off here
 ##############################################################
 
## kinematic function for Hortonian overland flow 
 Beta=scalar(0.6);
 N=0.100; 					# Manning's N;
 ChannelN=0.100; 				# channel's Manning value
 ChannelWidth=5; 				# accuflux.map;

## Soil water characteristic curves of soil layers/deposits
 Tetamax1=0.46;              		# maximal moisture content unsaturated layers
 Tetamax2=0.46;
 Tetamax3=0.46;
 Tetamax4=0.46;
 Tetar=0.08;                 		# minimal moisture content common for all 4 layers
 ha=0.07;                               # Air entry value in meters (1m=10kPa)
 alpha=0.1;                             # Slope of SWCC 
 Frh1=0.3;                   		# fractional depth of layer 1 unsaturated zone
 Frh2=0.3;					# fractional depth of layer 2 unsaturated zone
 Frh3=0.2;					# fractional depth of layer 3 unsaturated zone
 Frh4=0.2;					# fractional depth of layer 4 unsaturated zone
 Loss=0;                     		# loss of groundwater to rock in cm per day, impermeable if 0

## Mechanical propertis of soil/deposits for slope stability analysis
 TanPhi=0.66;               		# tangens of the angle of internal;  phi=35 degrees
 VarPhi=0.5;                   		# Variation of tangens phi
 Cohes=2.5;               	# soil cohesion  kN/m2
 VarC=0.5;                     		# Variation of cohesion
 Bulk=26;                    		# unit weight of soil kN/m3

## Statistical parameters
 C1=1.2533141;                 		# Costant for cumulative Probability function
 C2=0.5792933;                 		# Costant for Gaussian distribution
 Pi=3.1415927;                 		# Pi value

## Reporting
 Waterheight=waterh;        		# maps with waterheights in cm
 ProbFail=probf;            		# maps with probability of failure  
 CumFail=cumf;              		# number of days P(F<1)>0.5
 Unstab=unstab;             		# cells in map becoming unstable for more than 1 timestep
 LsMonitor=lssample1.map;   		# 3 monitor points for groundwaterheight
 LsSample1=lssample2.map;     		# monitor point for counting of unstable area
 LsSample2=lssample3.map;     		# monitor point for counting of unstable area
areamap
TotArea;

timer
1 365 1;                      		# time zero, total duration and timesteps 
anim=1,1+1..endtime;			# animation time to call history at specified intervals	

initial
## Initial conditions for hydrology and mechanical modules
## Distribution of cohesion and angle of internal friction
 #TanPhi=lookupscalar(friction.tbl, SoilType);
 #Cohes=lookupscalar(cohesion.tbl, SoilType);

## DEM derived parameters
 S=slope.map;                           # create a slope map from the DEM
 # S=if(Dem gt 350, S+5,S); 
 B=celllength();                        # pixelsize m
 CL=celllength();				      	# celllenght (cm) 
 CA=celllength()*celllength();				      	# celarea cm2
 H=soildepth.map; 				# if(S>0.7,700,soildepth.map);       		# soil depth cm (is 0 if slope higher as soil friction angle)
 TA=areatotal(TotArea,Clone);           # Calculating the total amount of cells
 CumFail=1;                            	# Days the slope is unstable;

## Timesteps to divide rainfall 
 Raindur=1/T;                        	# number of timesteps with rain during one day 
 Tn=0;                                  # counting timesteps
## Hydrological conditions
 #InterCept=lookupscalar(IntercTable,Veg); # Interception loss by vegetation
 # Ksat1=lookupscalar(SoilInfiltrationTable, SoilType); #creating a Ks map (mm/day) with landcover map
 Ksat1=150; # if(H le 2000, 15, Ksat1); 
 Ksat2=100; 					# lookupscalar(SoilInfiltrationTable, SoilType);
 Ksat3=75; 					# lookupscalar(SoilInfiltrationTable, SoilType); 
 Ksat4=50; 					# lookupscalar(SoilInfiltrationTable, SoilType);

## Conditions for upstream contributed runoff 
## term for Alpha
 AlphaFact=(ChannelN)/(sqrt(S))**Beta;
 AlphaPower=(2/3)*Beta; #Power for Alpha;
 WH=0.0000001;
 FlowWidth= ChannelWidth;
 Q=0.0001;
 nrTimeSlices=1;
 DCL=max(downstreamdist(Ldd),CL);
 Q0=scalar(0.00000001);
 
## Coverage of meteorological stations for the whole area
 RainZones=rainzone1.map; #spreadzone(RainStations, 0, 1);

## Moisture content at beginning of the simulation volumetric (cm3/cm3)
 Moisturecont1=0.35; 			# theta1_ini.map;			# Values derived from 10 years simulation
 Moisturecont2=0.35; 			# theta2_ini.map;			# Values derived from 10 years simulation
 Moisturecont3=0.35; 			# theta3_ini.map;			# Values derived from 10 years simulation
 Moisturecont4=0.35; 			# theta4_ini.map;			# Values derived from 10 years simulation
 Sin=sin(atan(S));
 Cos=cos(atan(S));
 CumFail=0;                  		# initial value of unstable days
 #Waterheight=inwaterh.map;  		# map with initial waterheight
 Waterheight=2;
 Pure_Waterheight=2;

dynamic

## Calculation of new depth (cm) of unsaturated layers (0 if H=0)
    H1=max(5,(H-Waterheight)*Frh1);
    H2=max(5,(H-Waterheight)*Frh2);
    H3=max(5,(H-Waterheight)*Frh3);
    H4=max(5,(H-Waterheight)*Frh4);

## Counting timesteps
    Tn=Tn+1;

## Atmopshere-ground interactions
## air temperature 
 # ta=timeinputscalar(avgtemp,RainZones);
 # tcrit=1;					# calibrated temperature threshold for precipitation to fall as snow 
## precipitation in cm per timestep (include rainfall and snowfall)
 # Precip=timeinputscalar(rainint.tss,boolean(1))*T;
## calculate and report maps with rainfall at each timestep (mm/day)
 Precip=10*timeinputscalar(Rain,RainZones); # if(ta>0,timeinputscalar(Rain,RainZones)*1.25,0); 
 # Precip=if(Precip==0, 10, Precip); #Rain; # if(Tn+1<Raindur,Rain*T,0);
 # Snow=if(ta<=tcrit,timeinputscalar(Rain,RainZones),0);	# precipitation as snow if air temperature <=1 degree celcius
 # apot=if(ta>tcrit,ta*DDF,0); 		# potential snowmelt using degree day approach
 # SD=timeinputscalar(snowdepth, RainZones);
 # SWE=SD*(0.5/0.98)+Snow;			# Snow-water equivalent from snowdepth, density of snow and pure water
 # ground temperature estimation
 # tg0=ta*1.15;				# estimation of initial ground temperature from air temperature
 # fs=0.2;					# empirical damping parameter
 # tg=tg0*exp(-fs*SD/100); 		# estimation of soil temperature beneath snow cover
 tg=36;
 # aact=3.75*min(apot,SWE);		# actual snowmelt
 # asw=if(tg>0,Precip+aact,0);		# available surface water from rainfall and snowmelt Above 0 degrees, raining on snow not considered
 asw=Precip;
 # Calculation of net rain
 NetPrecip=if(asw>0,asw-T,0);

 ################################
 # output runoff at each timestep
 RunOff=if(NetPrecip<Ksat1,0,max(NetPrecip-Ksat1,0));
 # compute both runoff and actual infiltration
 # RunOff= accuthresholdflux(Ldd,asw,Ksat1);
 # Infiltration=accuthresholdstate(Ldd,asw,Ksat1);
 # runoff calculation using one of kimematic
 ALPHA = AlphaFact*((FlowWidth+2*WH)**AlphaPower);
 QIn=RunOff*CA/T/DCL;
 Qr = kinematic(Ldd,Q0,QIn,ALPHA,Beta,nrTimeSlices,T,DCL);
 V=Qr/CA; # flow velocity(cm3/s), CA is cell area
 WH=if(FlowWidth >0.001, (ALPHA*(Qr**Beta))/FlowWidth,0); #wh in cm unit
 ################################

## Fluxes in soil
 #Rain which can infiltrate into the soil
 #Upstream runoff which can infiltrate into the soil
  Infil_runoffupstream=if(WH<Ksat1,WH,0);
  Raininput=asw+Infil_runoffupstream-RunOff; # if(H eq 0,0,NetPrecip);
 #calculation of downwards percolation in unsatured layers in cm per timestep
  Percolation1=if(tg<0,0,Ksat1*
  ((Moisturecont1-Tetar)/(Tetamax1-Tetar))**8.5)*T;
  Percolation2=if(Moisturecont2<0.25,0,Ksat2*
  ((Moisturecont2-Tetar)/(Tetamax2-Tetar))**8.5)*T;
  Percolation3=if(Moisturecont3<0.25,0,Ksat3*
  ((Moisturecont3-Tetar)/(Tetamax3-Tetar))**8.5)*T;
  Percolation4=if(Moisturecont4<0.25,0,Ksat4*
  ((Moisturecont4-Tetar)/(Tetamax4-Tetar))**8.5)*T;

## Calculation of new moisture content in layers 

  Moisturecont1=if(H1 eq 0,0,max(min(Tetamax1,Moisturecont1+(Raininput-Evap*T-Percolation1)/H1),Tetar));
  Moisturecont2=if(H2 eq 0,0,max(min(Tetamax2,Moisturecont2+(Percolation1-Percolation2)/H2),0.25));
  Moisturecont3=if(H3 eq 0,0,max(min(Tetamax3,Moisturecont3+(Percolation2-Percolation3)/H3),0.25));
  Moisturecont4=if(H4 eq 0,0,max(min(Tetamax4,Moisturecont4+(Percolation3-Percolation4)/H4),0.25));

  Tetae1=(H1*Moisturecont1-Tetar)/(H1*Tetamax1-Tetar);
  Tetae2=(H2*Moisturecont1-Tetar)/(H2*Tetamax1-Tetar);
  Tetae3=(H3*Moisturecont1-Tetar)/(H3*Tetamax1-Tetar);
  Tetae4=(H4*Moisturecont1-Tetar)/(H4*Tetamax1-Tetar);

  AverMoistCont=Frh1*Moisturecont1+ Frh2*Moisturecont2+Frh3*Moisturecont3+ Frh4*Moisturecont4;
  AverTetae=Frh1*Tetae1+ Frh2*Tetae2+Frh3*Tetae3+ Frh4*Tetae4;
  # AverTetae=(Tetae1+Tetae2+Tetae3+Tetae4)/4;

##ROUTING GROUNDWATER
 #discharge in cm3 of water out of pixel per timestep
  # Q=Ksat4*Sin*B*Waterheight*T;
 #change waterheight through outflow pixel
  # Deltawaterheight=Q/(B*B*(Tetamax4+0.001-Moisturecont4));
 #calculation inflow to pixel in terms of waterheight change
  # Inflowdeltawaterheight=Deltawaterheight; # Inflowdeltawaterheight=upstream(Deltawaterheight);
 #New waterheight by percolation and inflow,outflow of groundwater
  # Waterheight=if(H eq 0,0,max(0,Waterheight+Inflowdeltawaterheight-Deltawaterheight+(Percolation4-Loss*T)/(Tetamax4+0.01-Moisturecont4)));
  # Waterheight=min(H,Waterheight);
  # Waterheight=max(0,Waterheight);

##ROUTING GROUNDWATER Santy X-Y
 # Fraction of water to be routed towards x and 
 # y direction (-)(FractToX pos to N and E) 
 FractToX = Sin / (abs(Sin)+abs(Cos)); 
 FractToY = Cos / (abs(Sin)+abs(Cos));
 # discharge in cm3 of pure water out of pixel per timestep
 Q=(Ksat4*Sin*Waterheight*CL*T)/100;		#moving volume of water m3

 #Here we determine the amount of pure water cm3 flowing 
 #out of the central cell in a X and Y direction 
 Q_x=Q*FractToX;
 Q_y=Q*FractToY;

 #Calculate new pure water height in the central cell due 
 #to out and inflow (cm) 
 Pure_Waterheight=Pure_Waterheight- abs(Q_x)/CA - abs(Q_y)/CA + 
 max(0, shift0(Q_x, 0, -1))/CA + max(0, shift0(-Q_x, 0, 1))/CA + 
 max(0, shift0(Q_y, 1, 0))/CA + max(0, shift0(-Q_y,-1, 0))/CA; 
 Pure_Waterheight =if(Pure_Waterheight le 0,0, Pure_Waterheight);

 #Calculation of new water height in soil 
 Waterheight=max(0,Pure_Waterheight/(Tetamax4+0.01-AverMoistCont));
 #New soil water height caused by percolation and loss of water in underground
 Waterheight=if(H eq 0,0,max(0, Waterheight+WH+(Percolation4-Loss*T)/(Tetamax4+0.01-Moisturecont4)));

 Waterheight=min(H,Waterheight);
 Waterheight=max(0,Waterheight);

 ## Calculation of Slope stability (dimension is mtrs now)
  #matric suction_layer 1
  hs1=ha*exp(alpha*(1-Tetae1));
  phib1=57.3*TanPhi*((Moisturecont1-Tetae1)/(Tetamax1-Tetar));
  deltac1=hs1*(phib1); # *-10;				# Apparent cohesion in m
  #matric suction_layer 2
  hs2=ha*exp(alpha*(1-Tetae2));
  phib2=57.3*TanPhi*((Moisturecont2-Tetae2)/(Tetamax2-Tetar));
  deltac2=hs2*(phib2); # *-10;				# Apparent cohesion in m
  #matric suction_layer 3
  hs3=ha*exp(alpha*(1-Tetae3));
  phib3=57.3*TanPhi*((Moisturecont3-Tetae3)/(Tetamax3-Tetar));
  deltac3=hs3*(phib3); # *-10;				# Apparent cohesion in m
  #matric suction_layer 4
  hs4=ha*exp(alpha*(1-Tetae4));
  phib4=57.3*TanPhi*((Moisturecont4-Tetae4)/(Tetamax4-Tetar));
  deltac4=hs4*(phib4); # *-10;	
  hs=hs1+hs2+hs3+hs4; 					# Apparent cohesion in m
  deltac=1-((deltac1+deltac2+deltac3+deltac4));  
  #porepressure
  Porepr=(((Waterheight)/100)*10*sqr(cos(atan(S))));	
  #Safety Factor
  # Safety=if(H eq 0,2,if(S gt 0.4,1.5,min(1.5,(Cohes+(deltac)+(((H/100)*Bulk*sqr(cos(atan(S))))-Porepr)*TanPhi)/((H/100)*Bulk*sin(atan(S))*cos(atan(S))))));
  Safety=min(1.5,(Cohes+deltac+(((H/100)*Bulk*sqr(cos(atan(S))))-Porepr)*TanPhi)/((H/100)*Bulk*sin(atan(S))*cos(atan(S))));
  # Safety=if(Tn gt 190, Safety-0.2, Safety); 
 ##Calculation of Probability of Failure
  # Mean capicity
  Mcap=(Cohes+(((H/100)*Bulk*sqr(cos(atan(S))))-Porepr)*TanPhi);
  #whereby effective normal pressure P':
  P=((H/100)*Bulk*sqr(cos(atan(S))))-Porepr;
  # Mean demand
  TM=((H/100)*Bulk*sin(atan(S))*cos(atan(S)));
  #Mean safety factor(maximum 2)
  MSF=if(TM>0,min(2,Mcap/TM,2));
  #Probability of failure
  MSM=Mcap-TM;
  VTPhi=P**2*VarPhi;
  SSMP=sqrt(VarC+VTPhi);
  ZScore=if(H eq   0,0,MSM/SSMP);
  Distribution=scalar(atan(abs(ZScore)*(C1+C2*ZScore**2)));
  Distribution=Distribution/Pi+0.5;
  #Probability F<1 (0,1)
  ProbFail=if(H eq 0,0,if(ZScore>0,1-Distribution,Distribution));

  #number of days with unstable slope segment(P(F<1)>0.5)
  CumFail=if(ProbFail>0.5,(CumFail+T),CumFail);
  #cells which become unstable for more than 1 timestep
  Unstab=if(CumFail>T,1,0);
  #percentage of area becoming unstable
  PercUnstab=(areatotal(Unstab,Clone)/TA)*100;
   
##Reports maps
report waterh.map=Waterheight;
report probf.map=ProbFail;
report cumf.map=CumFail;
report unstab.map=Unstab;
report PercUnstab.map=PercUnstab;
#report Slope=S;
report MeanSafety.map=MSF;
report MeanCapacity.map=Mcap;
report MeanDemand.map=TM;
report MeanSafetyMargin.map=MSM;
report SSMP.map=SSMP;
report ZScore.map=ZScore;
report Porepr.map=Porepr;
report effectivestress.map=P;
report Ksat1.map=Ksat1;
report rainfall=Precip;
report(anim) fos.map=Safety;
report H.map=H;
report(anim) fos = Safety;
report(anim) probf =ProbFail;
report rainzone.map=RainZones;

# report tg=tg;
# report snow=Snow;
# report snowmelt=aact;
# report SWE=SWE;
#report asw=asw;
#report runoff=RunOff;
#report percol=Percolation1;
#report percol2=Percolation2;
#report percol3=Percolation3;
#report theta1=Moisturecont1;
#report theta2=Moisturecont2;
#report avgtheta=AverTetae;
#report thetae1=Tetae1;
#report thetae2=Tetae2;
#report thetae3=Tetae3;

##Reports timeseries
report Qr.tss=timeoutput(Outlet,Qr);
report V.tss=timeoutput(Outlet,V);
report WH.tss=timeoutput(Outlet,WH);
 #report tg1.tss=timeoutput(LsMonitor,tg);
 # report tg2.tss=timeoutput(LsSample1,tg);
# report tg3.tss=timeoutput(LsSample2,tg);
# report snow1.tss=timeoutput(LsMonitor,Snow);
# report snow2.tss=timeoutput(LsSample1,Snow);
# report snow3.tss=timeoutput(LsSample2,Snow);
# report snowmelt1.tss=timeoutput(LsMonitor,aact);
# report snowmelt2.tss=timeoutput(LsSample1,aact);
# report snowmelt3.tss=timeoutput(LsSample2,aact);
# report SWE1.tss=timeoutput(LsMonitor,SWE);
# report SWE2.tss=timeoutput(LsSample1,SWE);
# report SWE3.tss=timeoutput(LsSample2,SWE);
#report asw1.tss=timeoutput(LsMonitor,asw);
#report asw2.tss=timeoutput(LsSample1,asw);
#report asw3.tss=timeoutput(LsSample2,asw);
#report runoff1.tss=timeoutput(LsMonitor,RunOff);
#report runoff2.tss=timeoutput(LsSample1,RunOff);
#report runoff3.tss=timeoutput(LsSample2,RunOff);
#report waterh1.tss=timeoutput(LsMonitor,min(Waterheight,H));
#report waterh2.tss=timeoutput(LsSample1,min(Waterheight,H));
#report waterh3.tss=timeoutput(LsSample2,min(Waterheight,H));
#report p_unstab1.tss=timeoutput(LsMonitor,PercUnstab);
#report p_unstab2.tss=timeoutput(LsSample1,PercUnstab);
#report p_unstab3.tss=timeoutput(LsSample2,PercUnstab);
#report precip1.tss=timeoutput(LsMonitor,Precip);
#report precip2.tss=timeoutput(LsSample1,Precip);
#report precip3.tss=timeoutput(LsSample2,Precip);
#report precip4.tss=timeoutput(RainStations,Precip);
#report pwp1.tss=timeoutput(LsMonitor,Porepr*10);			# PWP - Multiply by 10 unit conversion m to kPa
#report pwp2.tss=timeoutput(LsSample1,Porepr*10);			# PWP - Multiply by 10 unit conversion m to kPa
#report pwp3.tss=timeoutput(LsSample2,Porepr*10);			# PWP - Multiply by 10 unit conversion m to kPa
#report Q1.tss=timeoutput(LsMonitor,Q);
#report Q2.tss=timeoutput(LsSample1,Q);
#report Q3.tss=timeoutput(LsSample2,Q);
report fos1.tss=timeoutput(LsMonitor,Safety);
report fos2.tss=timeoutput(LsSample1,Safety);
report fos3.tss=timeoutput(LsSample2,Safety);
report probf1.tss=timeoutput(LsMonitor,ProbFail);
report probf2.tss=timeoutput(LsSample1,ProbFail);
report probf3.tss=timeoutput(LsSample2,ProbFail);
#report percol1.1.tss=timeoutput(LsMonitor,Percolation1);
#report percol1.2.tss=timeoutput(LsMonitor,Percolation2);
#report percol1.3.tss=timeoutput(LsMonitor,Percolation3);
#report percol1.4.tss=timeoutput(LsMonitor,Percolation4);
#report percol2.1.tss=timeoutput(LsSample1,Percolation1);
#report percol2.2.tss=timeoutput(LsSample1,Percolation2);
#report percol2.3.tss=timeoutput(LsSample1,Percolation3);
#report percol2.4.tss=timeoutput(LsSample1,Percolation4);
#report percol3.1.tss=timeoutput(LsSample2,Percolation1);
#report percol3.2.tss=timeoutput(LsSample2,Percolation2);
#report percol3.3.tss=timeoutput(LsSample2,Percolation3);
#report percol3.4.tss=timeoutput(LsSample2,Percolation4);
#report h1.1.tss=timeoutput(LsMonitor,H1);
#report h2.1.tss=timeoutput(LsSample1,H1);
#report h3.1.tss=timeoutput(LsSample2,H1);
#report h1.2.tss=timeoutput(LsMonitor,H2);
#report h2.2.tss=timeoutput(LsSample1,H2);
#report h3.2.tss=timeoutput(LsSample2,H2);
#report h1.3.tss=timeoutput(LsMonitor,H3);
#report h2.3.tss=timeoutput(LsSample1,H3);
#report h3.3.tss=timeoutput(LsSample2,H3);
#report hs1.tss=timeoutput(LsMonitor,hs*10);			# Matric suction - Multiply by 10 unit conversion m to kPa
#report hs2.tss=timeoutput(LsSample1,hs*10);			# Matric suction - Multiply by 10 unit conversion m to kPa
#report hs3.tss=timeoutput(LsSample2,hs*10);			# Matric suction - Multiply by 10 unit conversion m to kPa
#report deltac2.tss=timeoutput(LsSample1,deltac);			# Apparent vohesion - Multiple by 10 unit conversion m to kPa
#report deltac3.tss=timeoutput(LsSample2,deltac);			# Apparent vohesion - Multiple by 10 unit conversion m to kPa
report theta1.1.tss=timeoutput(LsMonitor,Moisturecont1);
report theta1.2.tss=timeoutput(LsMonitor,Moisturecont2);
report theta1.3.tss=timeoutput(LsMonitor,Moisturecont3);
#report theta1.4.tss=timeoutput(LsMonitor,Moisturecont4);
#report theta2.1.tss=timeoutput(LsSample1,Moisturecont1);
#report theta2.2.tss=timeoutput(LsSample1,Moisturecont2);
#report theta2.3.tss=timeoutput(LsSample1,Moisturecont3);
#report theta2.4.tss=timeoutput(LsSample1,Moisturecont4);
#report theta3.1.tss=timeoutput(LsSample2,Moisturecont1);
#report theta3.2.tss=timeoutput(LsSample2,Moisturecont2);
#report theta3.3.tss=timeoutput(LsSample2,Moisturecont3);
#report theta3.4.tss=timeoutput(LsSample2,Moisturecont4);
report avgtheta1.tss=timeoutput(LsMonitor,AverTetae);
report avgtheta2.tss=timeoutput(LsSample1,AverTetae);
report avgthetas3.tss=timeoutput(LsSample2,AverTetae);
#report thetae1.1.tss=timeoutput(LsMonitor,Tetae1);
#report thetae1.2.tss=timeoutput(LsMonitor,Tetae2);
#report thetae1.3.tss=timeoutput(LsMonitor,Tetae3);
#report thetae1.4.tss=timeoutput(LsMonitor,Tetae4);
#report thetae2.1.tss=timeoutput(LsSample1,Tetae1);
#report thetae2.2.tss=timeoutput(LsSample1,Tetae2);
#report thetae2.3.tss=timeoutput(LsSample1,Tetae3);
#report thetae2.4.tss=timeoutput(LsSample1,Tetae4);
#report thetae3.1.tss=timeoutput(LsSample2,Tetae1);
#report thetae3.2.tss=timeoutput(LsSample2,Tetae2);
#report thetae3.3.tss=timeoutput(LsSample2,Tetae3);
#report thetae3.4.tss=timeoutput(LsSample2,Tetae4);
# report Ksat11.tss=timeoutput(LsMonitor,Ksat1);
# report Ksat12.tss=timeoutput(LsMonitor,Ksat1);
# report Ksat21.tss=timeoutput(LsMonitor,Ksat2);
# report Ksat22.tss=timeoutput(LsMonitor,Ksat2);