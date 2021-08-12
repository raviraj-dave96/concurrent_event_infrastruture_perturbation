## Model for landslide spreading scouring and deposition. 
## c2011  Theo van Asch; c2021 Srikrishnan Siva Subramanian             
## All rights reserved                           
## Distributed for educational scientific purposes only 
## Any use of this code should be cited properly    
## From rain input, overlandflow, scouring to debris flows 
## Mass conservative, no complete momentum conservation.  
## Infinite slope model
## Apparent friction angle 
## Flexible timestep based on CFL condition.
## Erosion like Hungr
## Recalculating density and concentration 
## At certain concentration from manning to Bingahm  
## TOEVOEGEN SCOURING IN ROOT
## Initial moisture content, percolation, saturation, saturation exceeded runoff
binding 
######################################################
# /Model inputs (either as constants or as maps)/ 
# #####################################################
## Map bindings 
#timestep (seconds) 
TS              = 1;                        # time step duration in seconds
#Dem            = dtm_shuida_large.map;     # digital elevation model map (m) 
#Dem             = dtm_maliu_basin.map;      # Dem map
Dem             = dtm.map;
#Storm           = rain_maliu_basin.map;     # Map indicating  area with rain
Storm             =rain_gui.map;
#Storm          = rain_shuida_basin.map;
#ScourA         = slide_shuida05.map;       # erodable area 
#ScourA         = debris_shuida_basin.map;  # erodable area 
#Dfdraw          = debris.map;



#Monitor       = samplepoint.map;     # Local points for Monitoring impacts 
#Clone         = clone7.map;          # Clone 
Alpha          = alpha.map;               # Slope map base
Beta           = beta.map;                # Slope map debris flow  
CBed           = 0.65;                # concentration solids in bed (m/m) 
DGrain         = 0.0019;              # mean diameter grains (m)
RhoSolid       = 2600;                # mass density solids kg/m3
RhoWater       = 1000;                # mass density water kg/m3
PhiBed         = 25;                  # friction angle bed material 

##Input and model parameters 

Evap=0;                     			# evapotranspiration cm/day switched off here
              
 Tetamax1=0.46;              			# maximal moisture content unsaturated layers
 Tetamax2=0.46;
 Tetamax3=0.46;
 Tetamax4=0.46;
 Tetar=0.02;                 			# minimal moisture content common for all 4 layers
 ha=0.07;                                	# Air entry value in meters (1m=10kPa)
 alpha=0.1;                             	# Slope of SWCC 
 Frh1=0.3;                   			# fractional depthof layers unsaturated zone
 Frh2=0.3;
 Frh3=0.2;
 Frh4=0.2;
 Loss=0;                     			# loss of groundwater to rock in cm per day

N            = 0.04;         # Manning's n 
Klat         = 1;            # lateral stress coefficient (constant)
Grav         = 9.8;          # Gravity accelleration               #
GammaBed     = 18;           # Unit weight of the bed material +water (kN.m-3)   
GammaWat     = 10;           # Unit weight of water(kN m-3)
GammaSolid   = 26;           # Unit weight of debris solids from bed (kN-3) 
VerySmall    = 0.0012;       # threshold to avoid dividing by very small 
CritH        = 0.1;          # Critical height for scorerate an deposition 
###############################################################################
ErConst        = 0.1;     # Scour constant in erosion rate equation
DepConst       = 0.0001;      # Deposition constant in deposition rate equation
p              = 0.67;       # inertial constant in deposition rate equation
sb             = 1;         # degree of saturation bed material
Rain           = rainint.tss;         # Rain intensity in (m/hour).
 
RainDur      = 86400;         # Rain duration in sec 18 hr
    
                               
KsSoil        = 0.000333333;         # Hydraulic conductivity soil (m/hour)
##############################################################################
######################################
# /Numerical stability control/
######################################

CFLimsup    = 0.5;   #Higher value of Courant-Friedrichs-Levy               
CFLiminf    = 0.3;   #Lower value of Courant-Friedrichs-Levy
MinNLoops   = 1;     #Minimum number of internal loops
MaxNLoops   = 124;   #Maximum number of internal loops
InitialLoops= 1;     #inititial number of internal loops

##Reporting
 
  LsMonitor=lssample1.map;   		# 3 monitor points for groundwaterheight
  LsSample1=lssample2.map;     		# monitor point for counting of unstable area
  LsSample2=lssample3.map;     		# monitor point for counting of unstable area

areamap 
#dem2.map;
#dtm_maliu_basin.map; 
#dtm_shuida_basin.map; 
dtm.map;
timer 

1 86400 1;                      # number of timesteps 
rep_a = 1, 3600 + 3600..endtime;	    # reporting intervals (at t=1hr., then each
rep_b = 1, 3600+3600..endtime;
rep_c = 30000+60..endtime; 
#rep_b = 750 + 1..endtime;	    # n+n until endtime)
end   =  endtime;               #can be used for reporting only end map   

initial 

# Soil depth
H_s=soildepth.map;

#moisturecontent at beginning of the simulation volumetric (cm3/cm3)
 Moisturecont1=0.05;			# Values derived from 10 years simulation
 Moisturecont2=0.05;			# Values derived from 10 years simulation
 Moisturecont3=0.05;			# Values derived from 10 years simulation
 Moisturecont4=0.05;			# Values derived from 10 years simulation
 Tetae1=(Moisturecont1-Tetar)/(Tetamax1-Tetar); # Initial degree of saturation
 Sin=sin(atan(slope(Dem)));
 Cos=cos(atan(slope(Dem)));

#Initial ground temperature
 tg=36;

#Rain intensity m/hr    
 #  Rain =0.050;

 # coverage of meteorological stations for the whole area
 RainZones=spreadzone(Storm, 0, 1);

#Initial unit weight of flow (water) ( kN.m-3)
   GamDf=10;            
#Volumetric concentration of solids in flow (m/m)  
   CSF=spatial (0.001);                
#Initial thickness of excess rain
   Hwat=0;                   #introduir un valor inicial de 0.1 m
   Waterheight=0;
   Pure_Waterheight=0;
       
#Initial cumulative excess Rain for massbalance
   ExRcum=spatial(0);
#Initial thickness of total flow (debris) flow (m)  
   H =scalar(0.0001);
#Initial thickness of total flow (debris) flow (m)  
   HDf =0.00;
#Total discharge of sum of sediment and water 
    QTot=0;
#initial map with depth of erodible material 
    BedDepth          = scmat.map; 
    #BedDepth          = scmaterial.map;            
    #SliDep           = slide_rain.map;    
    #SliDep           = gullaccu.map;
    #SliDep           = material.map;
   # SliDepIni         =if(Dem gt 970, BedDepth,0);
  SliDepIni         =BedDepth;
 #total volume of  material at start simulation in the catchment 
   # MatStart=maptotal(if(Dem gt 970,SliDepIni*cellarea(),0));
  MatStart=maptotal(SliDepIni*cellarea());

#Initial scourrate per time step (m) 
   Scourat=0;
#Initial routed scoured material    
   Solids=0;
#Cumulative scour for mass balance 
   CumScou=spatial(0);
#Initial Velocity and discharge
   Vel = spatial(scalar(0));
   VelVis=0;
   VelVoe=0;
# Give the cell lenght of a pixel  
   CL = celllength();
   CA = celllength()*celllength();
#Number of internal loops
   Loops=InitialLoops;
#Initial impact of flow in timestep 0 
   HitOld=0;
#Counter
    n=0;
#Initial correction factor for mass balance 
    CorrFact=1;


################################################################################


# The rheology changes when it arrives in the plane due to consolidation effects 
    
    PhiDf=if(Dem gt 100,40,35); #Bingham
    CohesDf=if(Dem gt 100,1,5);#Bingham
    ViscosDf=if(Dem gt 100,1,1);#Bingham
    TanPhiDf=tan(PhiDf); #Bingham
    #Turcoef=if(Dem gt 930,20,2);#Voellmy
    #PhiDf=if(Dem gt 930,2,20);#Voellmy
    #TanPhiDf=tan(PhiDf); #Voelmy
    #Phibed=36; #Egishira erosion
    #TanPhiBed=tan(Phibed);#Egishira erosion
################################################################################          

 
dynamic 

## Calculation of new depth (cm) of unsaturated layers (0 if H=0)
    H1=max(5,(H_s-Waterheight)*Frh1);
    H2=max(5,(H_s-Waterheight)*Frh2);
    H3=max(5,(H_s-Waterheight)*Frh3);
    H4=max(5,(H_s-Waterheight)*Frh4);

#################################
# /Topographical issues /
#################################
#Shift the debris map   

#Absolute height flow surface  
    Z=Dem+H;

#Basal surface gradient

    #report 
    report (end) Alpha = min(35,scalar(atan(slope(Dem))));
    # report Alpha=if( Dem gt 980 and Alpha lt 20,20,Alpha); 
    #Alpha= if (Dem lt 960,0,Alpha); 
    TanAlpha = tan(Alpha);
    SinAlpha = sin(Alpha);
    CosAlpha = sin(Alpha);
  


# Flow surface gradient and aspect 
    #Beta = min(20,scalar(atan(slope(Z))));          # slope map () 
     #report 
     report(end) Beta    = min(35,scalar(atan(slope(Z))));
     SinBeta  = sin(Beta); 
     CosBeta  = cos(Beta); 
     TanBeta  = tan(Beta); 

# Aspect direction of the free surface 
    Theta    = aspect(Z);               
    Theta    = if(nodirection(Theta),1,Theta);
    SinTheta = sin(Theta); 
    CosTheta = cos(Theta); 
   
# Fraction of material to be routed towards x and 
# y direction (-)(FractToX pos to N and E) 
    FractToX = SinTheta / (abs(SinTheta)+abs(CosTheta)); 
    FractToY = CosTheta / (abs(SinTheta)+abs(CosTheta)); 
    
    #FractToX = gradx(Z)/ (abs(gradx(Z))+abs(grady(Z))); 
    #FractToY = grady(Z)/ (abs(gradx(Z))+abs(grady(Z))); 


# The tangens of the friction of the bed material  
     TanPhiBed=tan(PhiBed);

# Counter for Rain duration in seconds   
 n=n+TS;
 


##################################################################   
# /Repeat {  }until statement internal loop subdivide TS (timestep) 
# into a number of smaller timesteps to ensure stability/
###################################################################   
    
    Loops=InitialLoops;
    repeat{     
    
    CFLmax=0;
  # Initiate values internal loops ( wander wether we need extension loop)  
  # we did it not for all parameters   
    #H_loop=H;
    #Vel_loop=Vel;
    #Scourat_loop=Scourat;
   
#########################################################################
# /Calculation of amount of waterdepth for surface run-off (m, vertically)
######################################################################### 

# Local Rain intensity in m/hr (vertically)
Precip=timeinputscalar(Rain,Storm);   
#Precip=if(Storm,Rain,0);
#Excess rain m/sec    
    ExceR=max(0.000000001,(Precip-KsSoil)/3600);
 #ExceR=(Precip-KsSoil);
    #ExceR=if(n lt RainDur and Storm,(Precip-KsSoil)/3600,0);
#Cumultive excess rain for mass balance (no routing)  
    ExRcum=ExRcum+ExceR*TS/Loops;
# Water height  rain per timestep (routed)  
    Hwat=Hwat+ExceR*TS/Loops;
# Rainfall after runoff
    Raininput=(Precip)-ExceR;

#calculation of downwards percolation in unsatured layers in m per timestep
  Percolation1=if(tg<0,0,KsSoil*
  ((Moisturecont1-Tetar)/(Tetamax1-Tetar))**8.5)*TS;
  Percolation2=if(Moisturecont2<0.25,0,KsSoil*
  ((Moisturecont2-Tetar)/(Tetamax2-Tetar))**8.5)*TS;
  Percolation3=if(Moisturecont3<0.25,0,KsSoil*
  ((Moisturecont3-Tetar)/(Tetamax3-Tetar))**8.5)*TS;
  Percolation4=if(Moisturecont4<0.25,0,KsSoil*
  ((Moisturecont4-Tetar)/(Tetamax4-Tetar))**8.5)*TS;


## Calculation of new moisture content in layers 

  Moisturecont1=if(H1 eq 0,0,max(min(Tetamax1,Moisturecont1+(Raininput-Evap*TS-Percolation1)/H1),Tetar));
  Moisturecont2=if(H2 eq 0,0,max(min(Tetamax2,Moisturecont2+(Percolation1-Percolation2)/H2),0.25));
  Moisturecont3=if(H3 eq 0,0,max(min(Tetamax3,Moisturecont3+(Percolation2-Percolation3)/H3),0.25));
  Moisturecont4=if(H4 eq 0,0,max(min(Tetamax4,Moisturecont4+(Percolation3-Percolation4)/H4),0.25));

  Tetae1=(H1*Moisturecont1-Tetar)/(H1*Tetamax1-Tetar);
  Tetae2=(H2*Moisturecont1-Tetar)/(H2*Tetamax1-Tetar);
  Tetae3=(H3*Moisturecont1-Tetar)/(H3*Tetamax1-Tetar);
  Tetae4=(H4*Moisturecont1-Tetar)/(H4*Tetamax1-Tetar);

  # AverMoistCont=Frh1*Moisturecont1+ Frh2*Moisturecont2+Frh3*Moisturecont3+ Frh4*Moisturecont4;
  AverMoistCont=(Moisturecont1+Moisturecont2+Moisturecont3+Moisturecont4)/4;  
  AverTetae=Frh1*Tetae1+ Frh2*Tetae2+Frh3*Tetae3+ Frh4*Tetae4;
  

##ROUTING GROUNDWATER Santy X-Y

 # Fraction of water to be routed towards x and 
 # y direction (-)(FractToX pos to N and E) 
  FractToX_w = Sin / (abs(Sin)+abs(Cos)); 
  FractToY_w = Cos / (abs(Sin)+abs(Cos));
 #discharge in cm3 of pure water out of pixel per timestep
 Q=(KsSoil*Sin*Waterheight*CL*TS);		#moving volume of water m3

 #Here we determine the amount of pure water cm3 flowing 
 #out of the central cell in a X and Y direction 
 Q_x=Q*FractToX_w;
 Q_y=Q*FractToY_w;

#Calculate new pure water height in the central cell due 
 #to out and inflow (m) 
 Pure_Waterheight=Pure_Waterheight- abs(Q_x)/CA - abs(Q_y)/CA + 
 max(0, shift0(Q_x, 0, -1))/CA + max(0, shift0(-Q_x, 0, 1))/CA + 
 max(0, shift0(Q_y, 1, 0))/CA + max(0, shift0(-Q_y,-1, 0))/CA; 
 Pure_Waterheight =if(Pure_Waterheight le 0,0, Pure_Waterheight);

 #Calculation of new water height in soil 
 Waterheight=max(0,Pure_Waterheight/(Tetamax4+0.01-AverMoistCont));
 #New soil water height caused by percolation and loss of water in underground
 Waterheight=if(H_s eq 0,0,max(0, Waterheight+Hwat+(Percolation4-Loss*TS)/(Tetamax4+0.01-Moisturecont4)));

 Waterheight=min(H_s,Waterheight);
 Waterheight=max(0,Waterheight);

##############################################################    
## Scour rate  
###########################################################
#Maximum concentration according to Takahashi 
    MaxCSF=max(0,min(0.9*CBed,(RhoWater*TanAlpha)/
     ((RhoSolid-RhoWater)*(TanPhiBed-TanAlpha))));
     report(end) maxC.map=MaxCSF;
   
         #report(rep_b) MaxCSF=(RhoWater*TanAlpha)/
         #((RhoSolid-RhoWater)*(TanPhiBed-TanAlpha));
  
         #report(rep_b) MaxCSF=(TanPhiBed-TanAlpha);
####################################################################### 
#Scourat according to Hungr                                            #
   #Scourat=max(0,ErConst*Vel*H);                                      #
                                                                       #
#Scourate in m of solids according to Takahashi                        #
    Scourat=max(0,ErConst*(MaxCSF-CSF)*QTot/((CBed-MaxCSF)*DGrain));  #
                                                                       #
#Scourate in m of solids according to Takahashi                        #
    #Rho=(CSF*(RhoSolid-RhoWater)+RhoWater)/1000;                       #
    #YStr=CSF*(RhoSolid-RhoWater)/1000*Grav*H*CosTheta**2*TanPhiBed;    #
    #DrStr=Rho*Grav*H*SinTheta;                                         #
    #Scourat=max(0,ErConst*(DrStr-YStr));                               #
########################################################################
# Scouring only in scour area with loose material and above CritH 
# No scouring when material thickness is 0
     Scourat=if(H lt CritH,0,Scourat);
     Scourat=if(BedDepth le 0,0,Scourat);
     Scourat=if(CSF gt MaxCSF,0,Scourat);
# Critical velocity for deposition 
    Rho=CSF*(RhoSolid-RhoWater)+RhoWater;    #density of debris flow
    k1=CSF*(RhoSolid-RhoWater)*TanPhiBed/Rho;
    k2=(CBed/CSF)**(1/3)-1;
    report (rep_b) Vcrit=2/(5*DGrain)*(10*k1*Rho/(0.02*RhoSolid))**0.5*k2*H**(1.5);
#Rate of deposition is gt 0 if Vel lt Vcrit and CSF lt MaxCSF
    Deporat=DepConst*(1-Vel/(p*Vcrit))*(MaxCSF-CSF)/CBed*Vel;

#Deposition condition If condition fulfilled Deporat is negative 
    
     Deporat=if(Vel lt p*Vcrit and CSF gt MaxCSF,Deporat,0);
     report (rep_b) Deporat=if(H lt CritH,0,Deporat);
#Net erosion 
     NetEr=Scourat+Deporat;   #deposition rate always negative ?
#############################################################################     
#########################################################################
#Cumulative Scour for mass balance ( no routing)  
    #CumScou=max(0,CumScou+TS/Loops*NetEr);
#Bed change inclusive porosity due to scour in one timestep
    #BedChange=TS/Loops*NetEr*(GammaSolid/GammaBed);
               #SoilLoss=TS/Loops*Scourat*1/ConSolbed;
#New Solids in m 
    #Solids=min(H*0.9*CBed,max(0,Solids+TS/Loops*NetEr));
#New  depth of bed material after  scouring and (or) deposition
    #BedDepth=max(0,BedDepth-BedChange);
#Lax operation for BedDepth
    #AverBedDepth=(4*BedDepth+shift0(BedDepth,0,-1)+shift0(BedDepth,0,1)
    #+shift0(BedDepth,-1,0)+shift0(BedDepth,1,0))/8;
    #BedDepth=AverBedDepth; 
###################################################################################
    #new calculation of solid content (m) due to erosion        
    Solids=Solids+NetEr*(CSF+(1-CSF)*(AverMoistCont/0.46))*TS/Loops;
    Hold      =H;
	#new calculation of volumetric concentration of solids (-) in the flow 
    CSF       = if(H lt VerySmall, 0,(CSF*Hold+(Scourat*CSF-Deporat*CSF) 
         *TS/Loops)/H);
    # cumulative erosion (positive or negative ) in m  
    CumScou 	 = max(0,CumScou + TS/Loops*NetEr);
	#What is left of loose erodable material due to net erosion 
    BedDepth = max(0, BedDepth - TS/Loops*NetEr);


##########################################################################
#############################################################################









##############################################################    
## Scour rate Egashira in meters solid material per second 
## measured perpendicular. Switch off all comments under head "scour rate "above  
###########################################################
#Maximum concentration according to Takahashi 
    #MaxCSF=max(0,min(0.9*CBed,(RhoWater*TanAlpha)/
     #((RhoSolid-RhoWater)*(TanPhiBed-TanAlpha))));
     #report(end) maxC.map=MaxCSF;

#Defining theta(e)=BetaE from Egishira
    #DensDif=(RhoSolid-RhoWater)*CSF;
    #TanPhibed=tan(Phibed);
    #AlphaE=atan((DensDif/(DensDif+RhoWater))*TanPhibed);
    #TanAlphaE=tan(AlphaE);
 #Scouretate (m/sec) (negative sedimentation)     
     #Scourat=if( H lt CritH,0,
     #if(BedDepth le 0,0, ErConst*Vel*CBed*(TanAlpha-TanAlphaE)));
#No scouring any more above critial concentration
    #Scourat=if(CSF gt MaxCSF,0,Scourat); 
#Rename for reporting (see below)    
    #Scourate=NetEr;
#Cumulative Scour for mass balance ( no routing)  
    #CumScou=CumScou+TS/Loops*Scourat;
#soil loss inclusive porosity due to scour in one timestep
    #BedChange=TS/Loops*Scourat*(GammaSolid/GammaBed);
#Solids in m per timestep (routed) 
    #Solids=Solids+TS/Loops*Scourat;
#New Soil depth due to scouring
    #BedDepth=max(0,BedDepth-BedChange);

####################################################
# /Calculation of velocity  
################################################  
##Coulomb-Viscous approach with acceleration term or Manning; 

     
    Sf=if(H lt VerySmall,0,CosAlpha*CosAlpha*
    TanPhiDf+(1/(GamDf*H))
    *(3/2*CohesDf+(3*ViscosDf/H)*Vel));   
  
    VelVis=min(10,max(0,if(H eq 0,0, 
    Vel+TS/Loops*Grav*(SinAlpha*CosAlpha+Klat*TanBeta-Sf))));
    VelMan = (H**2/3*SinAlpha**1/2)/N;
    WF=min(1,CSF/0.4);#Weightfactor; 0.4 is maximum CSF for manning =0 
             #Vel=if (CSF lt 0.2, VelMan,VelVis); 
    Vel= (1-WF)*VelMan+WF*VelVis;
    Vel=if(Vel lt 0.01,0,min(10,Vel));
    Vel=if(H gt VerySmall,Vel,0);

     
##Voellmy approach with accelleration term  or Manning 
    #VelVoeOld=VelVoe;
    #Sf=TanPhiDf+(VelVoe*VelVoe)/(Turcoef*H);
    #VelVoe=max(0,VelVoe +TS/Loops*Grav*(CosAlpha*SinAlpha+TanBeta-Sf));
    #VelVoe =if(H gt VerySmall,VelVoe,VelVoeOld);    
    #VelMan = (H**2/3*SinAlpha**1/2)/N;
         #WF=CSF/0.9*CBed; #Weightfactor; 0.9*CBed is maximum CSF 
         #Vel= (1-WF)*VelMan+WF*VelVoe;
    #WF=min(1,CSF/0.4); #Weightfactor; 
    #Vel= (1-WF)*VelMan+WF*VelVoe;
         #Vel= if(CSF lt 0.2,VelMan,Vel);
         #Vel=min(10,max(0, Vel));
         #Vel=min(0.5*CL*Loops/TS,max(0, Vel));
    #Vel=if(H gt VerySmall,Vel,0);

#Lax operation for velocity
    AverVel=(4*Vel+shift0(Vel,0,-1)+shift0(Vel,0,1)
    +shift0(Vel,-1,0)+shift0(Vel,1,0))/8;
    Vel=AverVel; 
   

################################################################
#/ Mass balance and Routing and new concentration of the material
################################################################   

## routing of  Hwat in a timestep
    #H_Volwat = if(Hwat eq 0,0,Hwat*Vel*TS/Loops); 
    H_Volwat = Hwat*Vel*TS/Loops; 

# routing of the moving part: split into x- and y- directions (m) 
    H_Volwat_x  = H_Volwat * FractToX; 
    H_Volwat_y  = H_Volwat * FractToY; 
# move the flow: mass balance (m) 
   Hwat = Hwat -  
     abs(H_Volwat_x)/CL - 
     abs(H_Volwat_y)/CL + 
     max(0, shift0(H_Volwat_x, 0, -1))/CL + max(0, shift0(-H_Volwat_x, 0, 1))/CL + 
     max(0, shift0(H_Volwat_y, 1,  0))/CL + max(0, shift0(-H_Volwat_y,-1, 0))/CL; 
#H not negative      
     Hwat=max(0,Hwat); 


##routing of  Solids in a timestep
    H_VolSolids = if(Solids eq 0,0,Solids*Vel*TS/Loops); 
# routing of the moving part: split into x- and y- directions (m) 
    H_VolSolids_x  = H_VolSolids * FractToX; 
    H_VolSolids_y  = H_VolSolids * FractToY; 
# move the flow: mass balance (m) 
    Solids = Solids -  
     abs(H_VolSolids_x)/CL - 
     abs(H_VolSolids_y)/CL + 
     max(0, shift0(H_VolSolids_x, 0, -1))/CL + max(0, shift0(-H_VolSolids_x, 0, 1))/CL + 
     max(0, shift0(H_VolSolids_y, 1,  0))/CL + max(0, shift0(-H_VolSolids_y,-1, 0))/CL; 
#H not negative      
     Solids=max(0,Solids); 


#New total flow height and corresponding mass balance correction
    H =Solids+Hwat;
    VolRain   =maptotal(ExRcum*cellarea());
    VolSolids =maptotal(Solids*cellarea());
    Volend    =maptotal(H*cellarea()); 
    CorrFac=(VolRain+VolSolids)/Volend;
    Solids=CorrFac*Solids;
    Hwat=CorrFac*Hwat;
    H=Solids+Hwat;
    QTot=H*Vel;

 
#Lax operation for H and HWat
    A=shift0(H,0,-1);
    B=shift0(H,0,1);
    C=shift0(H,-1,0);
    D=shift0(H,1,0);
    
        A=if(A lt VerySmall,B,A);
        B=if(B lt VerySmall,A,B); 
        C=if(C lt VerySmall,D,C); 
        D=if(D lt VerySmall,C,D);
        
        Sum=A+B+C+D;  
        #AverH=if(H lt 0.01,H,(4*H+Sum)/8);
   #AverH=(4*H+Sum)/8;
   #H = AverH;

    A=shift0(Hwat,0,-1);
    B=shift0(Hwat,0,1);
    C=shift0(Hwat,-1,0);
    D=shift0(Hwat,1,0);
    
        A=if(A lt VerySmall,B,A);
        B=if(B lt VerySmall,A,B); 
        C=if(C lt VerySmall,D,C); 
        D=if(D lt VerySmall,C,D);
        
        Sum=A+B+C+D;  
        #AverHwat=if(Hwat lt 0.01,Hwat,(4*Hwat+Sum)/8);
   #AverHwat=(4*Hwat+Sum)/8;
   #H = AverHwat;

                    
###########################################################
#Calculation of  concentration and Specific weight 
###########################################################

#New CSF flow
   report(rep_b) CSF=if(H lt VerySmall,0.001,max(0.001,min(0.9*CBed,Solids/H)));
#New Gamma   flow
    GamDf=if(H lt VerySmall,10,
    (Solids*GammaBed+Hwat*GammaWat)/H);

 
####################################
# /Numerical stability analyses/ 
####################################

#calculation of Cour Fr Levy condition 
     CFL=TS/Loops*sqrt(2)*mapmaximum(Vel)/CL;
     CFMap=TS/Loops*sqrt(2)*(Vel)/CL;
     CFLmax=max(CFLmax,CFL);
#update number of internal loops 
     IsUnstable=if(CFLmax ge CFLimsup,boolean(1),boolean(0));
#try wit more loops if unstable 
      Loops=if(IsUnstable,Loops+1,Loops); 
# amount of loops not more than going under minimum CFL 
      Loops=if(CFLmax le CFLiminf,max(Loops-1,MinNLoops),Loops);   
#exit loop
      } until not IsUnstable or Loops gt MaxNLoops;


#####################################      
# /End of internal loop/ 
####################################

#reports 
     report(rep_b) neter=NetEr;
     #report(rep_b) solids=Solids;
     #report (rep_b) hwat=Hwat;   # depth of  fluid component 
     report (rep_b) h=H;   # depth of flow of any concentration 
     report(rep_c) vel=Vel; #velocity report 
      # time report map with remaining loose in situ aterial
     report(rep_b) bedz=BedDepth; 
     # Thickness of debris flow in case of deposition of solid materialin flow  
     # assuming a certain porosity 
      Dfdepo=if(CSF gt 0.2,CSF*H*1/CBed,0); 
      report(end) dfdepo.map=Dfdepo;
              #report(end) phi.map=scalar(PhiDf);
              #report Dfdep=if(CSF gt 0.2,CSF*H*GammaSolis/GammaBed,0);
     
     #Thickness of Df debris flow including all the water 
     report (rep_c) Df=if(CSF gt 0.2,H,0);
     
     #material map with erodible material after stop of debris flow which is used 
     #as input for the next event
     report (end) ermat.map=Dfdepo+BedDepth;
     
     #How many internal loops per timestep? 
     report loops.tss=Loops;      

#Make a boolean map of the debris flow 
       
       report (end) dfdraw.map=scalar(if(H gt 0 and CSF gt 0.2,1,0));

###################################################################
#/Maximum Impact and tss Impact for monitor points 
#####################################################################
#Impact pressure kPa 
   #kConst=1.261*exp(CSF);
   #ImPact=kConst*0.1*GamDf*Vel*Vel;
#Static pressure 
    #ImPact=GamDf*H;
    #Hit=windowaverage(ImPact,3*celllength());
    #HitMax=max(HitOld,Hit);
    #report HitMax.tss=timeoutput(Monitor,HitMax);
    #HitOld=HitMax;
    

#################
#/Mass balance 
####################

#Calculation of total volume accumulated rain and non routed scoured material 
#plus total volume of routed material plus fan deposit volumes   

VolRain   =maptotal(ExRcum*cellarea());
VolSolids =maptotal(Solids*cellarea());
Volend    =maptotal(H*cellarea()); 

 volrain.tss=VolRain;
 solids.tss=VolSolids;
 volend.tss=Volend;
 report massbal.tss=(VolRain+VolSolids)/Volend;

####################
#Deposition 
###################

#consolidated deposition in different zones 
#Calculated thickness is hypothetical only valid at end of simulation 
#But we need a tss plot to show when the dposition becomed constant
#For the routing the value of Df is used. 

DepoDf1 =if(Dem gt 400 and Dem le 500,Dfdepo,0);
DepoDf2 =if(Dem gt 300 and Dem le 400,Dfdepo,0);
DepoDf3 =if(Dem gt 200 and Dem le 300,Dfdepo,0);
DepoDf4 =if(Dem gt 100 and Dem le 200,Dfdepo,0);
DepoDf5 =if(Dem gt 20 and Dem le 100,Dfdepo,0);
DepoDf6 =if(Dem gt 0 and Dem le 500,Dfdepo,0);
DepoDf7 =if(Dem gt 0 and Dem le 1000, Dfdepo, 0); 
DepoDf8 =if(Dem gt 0 and Dem le 2000, Dfdepo, 0);
DepoDf9 =if(Dem gt 0, Dfdepo, 0);



DepoVol1 =maptotal (DepoDf1*cellarea());
DepoVol2 =maptotal (DepoDf2*cellarea());
DepoVol3 =maptotal (DepoDf3*cellarea());
DepoVol4 =maptotal (DepoDf4*cellarea());
DepoVol5 =maptotal (DepoDf5*cellarea());
DepoVol6 =maptotal (DepoDf6*cellarea());
DepoVol7 =maptotal (DepoDf7*cellarea()); 
DepoVol8 =maptotal (DepoDf8*cellarea()); 
DepoVol9 =maptotal (DepoDf9*cellarea()); 
#report VolGul.tss = maptotal (gullaccu.map*cellarea());

 #MatEnd= maptotal(if(Dem gt 970,BedDepth*cellarea(),0));  
#total volume of displaced material by debrisflow , which remained in the catchment 
# MatRem= maptotal(if(Dem gt 970,Dfdepo*cellarea(),0));
report FanVol1.tss=DepoVol1;
report FanVol2.tss=DepoVol2;
report FanVol3.tss=DepoVol3;
report FanVol4.tss=DepoVol4;
report FanVol5.tss=DepoVol5;
report FanVol6.tss=DepoVol6;
report FanVol7.tss=DepoVol7;
report FanVol8.tss=DepoVol8;
report FanVol9.tss=DepoVol9;
report FanVol10.tss=timeoutput(catchments.map,maptotal(DepoDf1*cellarea()));
#report FanVol23.tss=DepoVol2+DepoVol3;

#timeplot of net erosion in catchment (positive) during simulation 
 #report DiffEros.tss=MatStart-MatEnd+MatRem;


#Correction factor of mass due to routing over irrigular DEM 
#and to lax operations. 
 #VolRain =maptotal(ExRcum*cellarea());
 #VolScou =maptotal(CumScou*cellarea());
 #Volend  =maptotal(H*cellarea());
#CorrFac=(VolRain+VolScou)/Volend;
#CorrFac= min(max(CorrFac,0.9),1.1);
#Solids=CorrFac*Solids;
#Hwat=CorrFac*Hwat;
#H=Solids+Hwat;

report DfDepo.1.tss=timeoutput(LsMonitor,Dfdepo); 
report DfDepo.2.tss=timeoutput(LsSample1,Dfdepo);
report DfDepo.3.tss=timeoutput(LsSample2,Dfdepo);
report FVol1.tss=timeoutput(LsMonitor,Dfdepo*cellarea()); 
report FVol2.tss=timeoutput(LsSample1,Dfdepo*cellarea());
report FVol3.tss=timeoutput(LsSample2,Dfdepo*cellarea());
#report theta1.1.tss=timeoutput(LsMonitor,Moisturecont1);
#report theta1.2.tss=timeoutput(LsMonitor,Moisturecont2);
#report theta1.3.tss=timeoutput(LsMonitor,Moisturecont3);
#report theta1.4.tss=timeoutput(LsMonitor,Moisturecont4);
#report theta2.1.tss=timeoutput(LsSample1,Moisturecont1);
#report theta2.2.tss=timeoutput(LsSample1,Moisturecont2);
#report theta2.3.tss=timeoutput(LsSample1,Moisturecont3);
#report theta2.4.tss=timeoutput(LsSample1,Moisturecont4);
#report theta3.1.tss=timeoutput(LsSample2,Moisturecont1);
#report theta3.2.tss=timeoutput(LsSample2,Moisturecont2);
#report theta3.3.tss=timeoutput(LsSample2,Moisturecont3);
#report theta3.4.tss=timeoutput(LsSample2,Moisturecont4);
report avgtheta1.tss=timeoutput(LsMonitor,AverMoistCont);
report avgtheta2.tss=timeoutput(LsSample1,AverMoistCont);
report avgtheta3.tss=timeoutput(LsSample2,AverMoistCont);
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
#report Ksat.tss=timeoutput(LsMonitor,KsSoil);
#report Scourat1.tss=timeoutput(LsMonitor,Scourat);
#report Scourat2.tss=timeoutput(LsSample1,Scourat);  
#report Scourat3.tss=timeoutput(LsSample2,Scourat);
#report Deporat1.tss=timeoutput(LsMonitor,Deporat);
#report Deporat2.tss=timeoutput(LsSample1,Deporat);  
#report Deporat3.tss=timeoutput(LsSample2,Deporat);
#report NetEr1.tss=timeoutput(LsMonitor,NetEr);
#report NetEr2.tss=timeoutput(LsSample1,NetEr);  
#report NetEr3.tss=timeoutput(LsSample2,NetEr);
   