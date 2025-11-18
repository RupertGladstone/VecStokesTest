
-- ##########################################################
-- ## Material parameters

yearinsec = 365.25 * 24.0 * 60.0 * 60.0
Pa2MPa = 1.0E-06
MPainPa = 1.0e6

-- ## rhoi_si = 917.0
rhoi_si = 910.0  -- ## ice density in kg/m^3
rhoi    = rhoi_si / (1.0e6 * yearinsec^2)

-- ## ocean water density 
rhoo_si = 1027.0 -- ## kg/m^3
rhoo = rhoo_si / (1.0e6 * yearinsec^2.0)

-- ## fresh water density
rhow = 1000.0 / (1.0e6 * yearinsec^2.0)

gravity = -9.81 * yearinsec^2
gravity_abs = -gravity
density_ratio = 1027.0 / 917.0
gravity_si = -9.81

-- ## Ice fusion latent heat
Lf_si = 335000.0  -- ## Joules per kg
Lf = Lf_si * yearinsec^2.0 -- ## some old Elmer units conversion thing?  Not used I think...

-- ## Sea level elevation
zsl = 0.0

-- ## Specific heat of ocean water (3974.0 Chen's value)(J/kg degCel or Kelvin)
-- ## Depends on ocean temp and salinity. is typically around 3850 J/(kg·K) at standard seawater conditions (35 PSU salinity, 25°C, 1 atm pressure).
-- ## So 4218 is basically fresh water. Lower means more saline or warmer.
cw_si = 3974.0 
cw = cw_si * yearinsec^2.0

-- ## prescribed salinity at ice base for calculating ocean pressure melting temperature.  PSU.
Salinity_b = 35.0 

--  For Glen's flow law (Paterson 2010)
n  = 3.0
ng = n
m  = 1.0/n
A1_SI = 3.985e-13 --- Pa^-3 s^-1 -> SI
A2_SI = 1.916e3
A1 = A1_SI*yearinsec*1.0e18  -- ## 1.2567e13 --- MPa^(-3) a^(-1)
A2 = A2_SI*yearinsec*1.0e18  -- ## 6.0423e28

Q1 = 60.0e3 -- ## J/mol for T > –10°C (warm ice)
Q2 = 139.0e3 -- ## J/mol for T < –10°C (cold ice)
Tlim = -10.0


-- ## GLToleranceInit=1.0e-1
GLTolerance=1.0e-3


--  Temperature of the simulation in Celsius
--  with the formula for A works only if T > -10
Tc=-1.0



-- ##########################################################
-- ## hard coded paths to forcing data

-- VELOCITY_DATA = "/projappl/project_2002875/data/antarctic/antarctica_m2slim.nc"
-- BETA_GUESS    = "/projappl/project_2002875/data/antarctic/aa_v3_e8_l11_beta.nc"
-- SMB_DATA      = "/projappl/project_2002875/data/antarctic/smbref_1995_2014_mar.nc"

datadir = "/g/data/jk72/cxz581/data/antarctic/"
outdir  = "./VTUoutputs"



-- ##########################################################
-- ## functions for material parameters and conditions

-- ## identify the grounding line in a 2D domain based on 
-- ## geometry (this is a crude approximation).
function groundingline(thick,bed,surf)
  if ((surf - thick) > (bed + 100.0) or (surf - thick) <= (bed + 1.0)) then
     gl_mask = -1.0
  else
     gl_mask = 1.0
  end
  return gl_mask
end


-- ## convert a viscosity enhancement factor to a 
-- ## flow enhancement factor, with limits.
function ConvertEF(viscEF,lowerLimit)
  if (viscEF < lowerLimit) then
    viscEF = lowerLimit
  end
  flowEF = 1.0 / viscEF
  if (flowEF < lowerLimit) then
    flowEF = lowerLimit
  end
  return flowEF
end

-- ## get second invariant of the strain rate
function calc_EII(Exx,Eyy,Ezz,Exy,Eyz,Ezx)
   tr = Exx+Eyy+Ezz
   EII = math.sqrt(0.5*( Exx*Exx + Eyy*Eyy + Ezz*Ezz+ 2.0*(Exy*Exy + Eyz*Eyz + Ezx*Ezx)- (1.0/3.0)*tr*tr ) )
   return EII
end

-- ## Paterson rate factor in Elmer/Ice units (MPa^-n a^-1)
function A_pat(T)
  Tc = T - 273.15
  if (Tc < -10.0) then
    A = A1 * math.exp(-Q1/(8.314*T))
  elseif (Tc > 0) then
    A = A2 * math.exp(-Q2/(8.314 * (273.15)))
  else
    A = A2 * math.exp(-Q2/(8.314*T))
  end
  return A
end

-- Build Glen enhancement E so that η_glen == η_inv at current εII
-- Inputs: eII [a^-1], mu [MPa·a], bottomEF [-], T [K]
-- ## Scale the drag coefficient to tune thinning rates...
function GlenE_from_power(mu, bottomEF,T)
  Tc = T - 273.15
  if (Tc < -10.0) then
    A = A1 * math.exp(-Q1/(8.314*T))
  elseif (Tc > 0) then
    A = A2 * math.exp(-Q2/(8.314 * (273.15)))
  else
    A = A2 * math.exp(-Q2/(8.314*T))
  end
  if (bottomEF == 0.0) then
    bottomeEF = 1.0e-3
  end
  n = 3.0
  eta_inv = mu * (bottomEF*bottomEF)   -- your inversion η = mu * EF^2
  if (eta_inv < 1.0e-12) then
    eta_inv = 1.0e-12
  end
  -- A_equiv that reproduces eta_inv at given eII
  E = 1/(2*A*(eta_inv^3))

  -- keep things sane (tune if needed)
  -- if (E < 0.01) then
  --   E = 0.01
  -- end
  if (E > 5.0) then
    E = 5.0
  end
  return E
end

function effective_viscosity_inv(bottomEF, mu, exx,eyy,ezz,exy,exz,eyz)
  -- calcuate the effective (second-invariant) strain rate
  ee = (0.5*(exx^2 + eyy^2 + ezz^2 + 2.0*(exy^2+exz^2+eyz^2)))^0.5
  -- visocity from inversion
  eta_inv = bottomEF^2 * mu
  effec_eta = eta_inv*(ee^(-1/3))
  return effec_eta
end

function ScaleDragCoef(coef)
  scaledCoef = coef*0.1
  return scaledCoef
end


-- ## for imposing a velocity condition based on temperature
-- ## input is temperature relative to pressure melting
function tempCondition(temp,tempCutoff)
  if (temp < tempCutoff) then
    cond = 1.0
  else
    cond = -1.0
  end
  return cond
end

-- ## Impose basal mass balance as lower surface normal velocity 
-- ## under shelf for steady simulations (e.g. inversions).
-- ## Note on sign:
-- ## Normal velocity is taken to be positive "outward" i.e. 
-- ## approx. downward under the ice shelf.
-- ## bmb is taken to be positive for mass gain and negative for 
-- ## mass loss.
-- ## So negative bmb => positive normal velocity
function bmb_as_vel(bmb,gmask)
    -- ##  temp = bmb
-- ##  if (isnan(bmb)) then
-- ##    temp = 0.0
-- ##  end
  if (gmask < -0.5) then
    vel = -1.0 * bmb -- chen: no minus 
  else
    vel = 0.0
  end
  return vel
end

-- ## more accurate identification of grounding line for a body
-- ## force condition if GroundedSolver is present.
-- ## Also checks velocity: allow coarse mesh for low speed GL.
function glCondition(glMask,vel,velCutoff)
  if ( (glMask < 0.5) and (glMask > -0.5) ) then 
    cond = 1.0
  else
    cond = -1.0
  end
  if ( vel < velCutoff ) then
    cond = -1.0
  end
  return cond
end

-- ## function to scale a normal velocity slip coefficient 
-- ## at the ice upper surface to restrict emergence 
-- ## velocity (stronger constraint in slow flowing regions).
function ControlEmergVel(vx,vy,upplim)
  ScalingSpeed = 10.0
  speed = math.sqrt(vx*vx + vy*vy)
  SlipCoef = upplim*(1.0 - math.tanh(speed/ScalingSpeed))
end

-- ## set the distance at GL to non-zero values depending on flow speed...
-- ## (experimental; not currently used)
function distBF(vel)
  if vel > refvel then
    dist = 0.0
  else
    dist = 100000.0 * (refvel - vel)/refvel
  end
  return dist
end

-- ## I have commented out this function as I am trying a modified version in the ASB.lua script. 
-- ## function for setting an upper limit to mesh size based on distance
-- ## (e.g. distance from grounding line)
-- #function refinebydist(distance)
-- #  factor = distance/GLdistlim
-- #  if (factor < 0.0) then
-- #    factor = 0.0
-- #  end    
-- #  if (factor > 1.0) then
-- #    factor = 1.0
-- #  end    
-- #  Mmax = Mmaxclose*(1.0-factor) + Mmaxfar*factor
-- #  return Mmax
-- #end

-- ## function for setting an upper limit for mesh size
function setmaxmesh(gldist,bdist,vel,glmask)
  gldistfactor = gldist/GLdistlim
  if (gldistfactor < 0.0) then
    gldistfactor = 0.0
  end    
  if (gldistfactor > 1.0) then
    gldistfactor = 1.0
  end    
  bdistfactor = bdist/Bdistlim
  if (bdistfactor < 0.0) then
    bdistfactor = 0.0
  end    
  if (bdistfactor > 1.0) then
    bdistfactor = 1.0
  end    
  if (gldistfactor < bdistfactor) then
    distfactor = gldistfactor
  else
    distfactor = bdistfactor
  end
  velfactor = vel/refvel
  if (velfactor > 1.0) then
    velfactor = 1.0
  end
  if velfactor < 0.5 then
    velfactor = 0.5
  end
  Mmax = ( Mmaxclose*(1.0-distfactor) + Mmaxfar*distfactor ) / (velfactor)
  if (glmask < 0.5) then
      if (Mmax > Mmaxshelf) then
      Mmax = Mmaxshelf
    end
  end
  --##if (gldist > distscale) and (Mmax < Mmaxfar) then --##and (bdist > 10000) then
  --##  Mmax = Mmaxfar --##/ distfactor
  --##end 
  return Mmax
end

-- ## function for setting a lower limit for mesh size
-- ## If gldist is smaller than GLdistlim, the effective distance is 0
-- ## which also makes hte distfactor 0, and same is bdist is within
-- ## Bdistlim. Why is it different for the max mesh?? for max mesh to get smoother transition?
function setminmesh(gldist,bdist,glmask)
  effectivedist = gldist - GLdistlim
  if (effectivedist < 0.0) then
    effectivedist = 0.0
  end
  distfactor = effectivedist / distscale
  if (distfactor > 1.0) then
    distfactor = 1.0
  end
  if (bdist < Bdistlim) then
    distfactor = 0
  end
  Mmin = Mminfine*(1.0-distfactor) + Mmincoarse*distfactor
  if (glmask < 0.5) then
    if (Mmin > (Mmaxshelf - 50.0) ) then
      Mmin = Mmaxshelf - 50.0
    end
  end
  --##if (gldist > distscale) and (Mmin < Mmaxfar) then --##and (bdist > 10000) then
  --##   Mmin = Mmaxfar
  --##end
  return Mmin
end

-- ## set the lower surface for a given upper surface and thickness
function getlowsurface(upp_surf,thick,bed)
  if (thick < MINH) then
    thick = MINH
  end
  if ((upp_surf - thick) > bed) then
    low_surf = upp_surf - thick
  else
    low_surf = bed
  end
  return low_surf
end

function getuppersurface(upp_surf)
  uppsurf = upp_surf
  return uppsurf
end      

function correctmask(groundedmask, low_surf)
  mask = groundedmask
  if (low_surf > 0 and mask <1) then
    mask = 1
  end 
    return mask
end

function getlowersurface(low_surf, bed, groundedmask)
  lowsurf = low_surf
  if (low_surf < bed) then
    lowsurf = bed
  end 
  return lowsurf
end         

-- ## set the upper and lower surfaces to floatation
function floatUpper(thick,bed)
--##  if (groundedmask > -1 and surface-bed < MINH) then
--##    upp_surf = bed + MINH
--##  end
  if (thick < MINH) then
    thick = MINH
  end
  if ( (thick*rhoi/rhoo) >= -(bed) ) then
    upp_surf = bed + thick
  else
    upp_surf = thick - thick*(rhoi/rhoo)
  end
  return upp_surf
end

function floatLower(thick,bed)
  if (thick < MINH) then
    thick = MINH
  end
  if ( (thick*rhoi/rhoo) >= -bed ) then
    low_surf = bed
  else
    low_surf = -thick*rhoi/rhoo
  end
  return low_surf
end     


-- ## variable timestepping (TODO: dt_init and dt_max and dt_incr should be passed in) 1.2 dt_max=0.25
function timeStepCalc(nt)
  dt_init = 0.000001 
  dt_max = 0.10
  dt_incr = 1.15
  dt = dt_init * 1.05^nt 
  if ( dt > dt_max ) then
    dt = dt_max
  end
  return dt
end

function timeStepCalc2(nt)
  dt_init = 0.0001
  dt_max = 0.1
  dt_incr = 1.15
  dt = dt_init * 1.05^nt
  if ( dt > dt_max ) then
    dt = dt_max
  end
  return dt
end

-- ## variable timestepping (TODO: dt_init and dt_max and dt_incr should be passed in) 1.2 dt_max=0.25
function timeStepCalc4(nt)
  dt_init = 0.0007 
  dt_max = 0.10
  dt_incr = 1.05
  dt = dt_init * 1.05^nt 
  if ( dt > dt_max ) then
    dt = dt_max
  end
  return dt
end



-- ## variable timestepping for restart 
function timeStepCalc5(nt)
  nt_start = 57 -- ## where the last one ended
  dt_init = 0.0001 
  dt_max = 0.10
  dt_incr = 1.15
  dt = dt_init * 1.05^(nt_start + nt) 
  if ( dt > dt_max ) then
    dt = dt_max
  end
  return dt
end



function timeStepCalc3(nt, dt_init, dt_max, dt_incr)
  dt = dt_init * dt_incr^nt 
  if ( dt > dt_max ) then
    dt = dt_max
  end
  return dt
end

-- ## thermal properties
function conductivity(Tin)
  T = Tin
  if (T > 273.15) then
    T = 273.15
  end
 k=9.828*math.exp(-5.7E-03*T)
  -- Apply bounds
  if (k < 1.0) then
    k = 1.0
  elseif (k > 5.0) then
    k = 5.0
  end
 return k
end

-- ## thermal properties
function conductivitycelc(Tin)
  T = Tin + 273.15
  if (T > 273.15) then
    T = 273.15
  end
 k=9.828*math.exp(-5.7E-03*T)
  -- Apply bounds
  if (k < 1.0) then
    k = 1.0
  elseif (k > 5.0) then
    k = 5.0
  end
 return k
end

-- ## heat capacity
function capacity(Tin)
  T = Tin
  if (T > 273.15) then
    T = 273.15
  end
  c=146.3+(7.253*T)
  if (c < 500.0) then
    c = 500.0
  elseif (c > 2000.0) then
    c = 2000.0
  end
  return c
end

-- ## heat capacity
function capacitycelc(Tin)
  T = Tin + 273.15
  if (T > 273.15) then
    T = 273.15
  end
  c=146.3+(7.253*T)
  if (c < 500.0) then
    c = 500.0
  elseif (c > 2000.0) then
    c = 2000.0
  end
  return c
end

function sw_pressure(z)
  if (z >  0) then
    p=0.0
  else
    p=-rhoo*gravity*z
  end
  return p
end

function init_velo1(v, g1, g2, zs, zb, z)
  gt = math.sqrt(g1*g1 + g2*g2)
  vin1=-v*(g1/gt)*z
  return vin1
end

function init_velo2(v, g1, g2, z)
  gt = math.sqrt(g1*g1 + g2*g2)
  vin2=-v*(g2/gt)*z
  return vin2
end

-- ## inputs: T is temperature in Kelvin.
-- ## returns temperature in Centigrade.
function K2C(T)
  Th= T - 273.15
  return Th
end  

-- ## inputs: TinC is temperature in Centigrade, p is pressure. 
-- ## Returns temperature relative to pressure melting point.
-- ## melting temperature decreases as the pressure increases
function relativetemp(TinC,p)
  pe = p
  if (pe < 0.0) then
    pe = 0.0
  end
  Trel = TinC + 9.8e-08*1.0e06*pe --- 9.8e-08 is the clausius clapeyron slope for ice
  if (Trel > 0.0) then
    Trel = 0.0
  end
  return Trel
end  

function pressuremelting(p)
  pe = p
  if (pe < 0.0) then
    pe = 0.0
  end  
  Tpmp = 273.15 - 9.8e-08*1.0e06*pe
  return Tpmp
end

function pressuremelting_salinity(p)
  pe = p
  if (pe < 0.0) then
    pe = 0.0
  end  

  Tpmp = 273.15 - 5.73e-02*Salinity_b + 9.39e-02 - 7.53e-08*1.0e06*pe

  return Tpmp
end

function initMu(TempRel)
  if (TempRel < Tlim) then
    AF = A1_SI * math.exp( -Q1/(8.314 * (273.15 + TempRel)))
  elseif (TempRel > 0) then
    AF = A2_SI * math.exp( -Q2/(8.314 * (273.15)))
  else
    AF = A2_SI * math.exp( -Q2/(8.314 * (273.15 + TempRel)))
  end
  glen = (2.0 * AF)^(-1.0/n)
  viscosity = glen * yearinsec^(-1.0/n) * Pa2MPa
  return viscosity
end
--  mu = math.sqrt(viscosity)
--  return mu

function limitslc(slc)
  slco = slc
  if (slco < 1.0e-06) then
    slco = 1.0e-06
  end
  return slco
end

-- Condition to impose no flux on the lateral side applied if
-- surface slope in the normal direction positive (should result in inflow)
-- and greater than 50m/km
function OutFlow(N1,N2,G1,G2) 
  cond=N1*G1 + N2*G2 - 0.05
  return cond
end

function evalcost(velx,vx,vely,vy)
  if (math.abs(vx)<1.0e06) then
     Cost1=0.5*(velx-tx(1))*(velx-vxy)
  else
     Cost1=0.0
  end
  if (math.abs(vy)<1.0e06) then
     Cost2=0.5*(vely-vy)*(vely-vxy)
  else
     Cost2=0.0
  end   
  return Cost1+Cost2
end

function evalcostder(vel,v)
  if (abs(v) < 1.0e06) then
    return (vel - v)
  else
    return 0.0
  end
end  

function initbeta(slc)
  dummy = slc + 0.00001
  return  math.log(dummy)
end


function lowerlimit(bed,groundedmask)
  hmin_watercolumn = 1.0
  if (groundedmask < 0) then
    return bed + hmin_watercolumn
  else
    return bed
  end
end

-- EF function using ocean connectivity
function effective_pressure(GroundedMask,thickness,BedElevation)
  if (BedElevation == nan) then
    BedElevation = 0.0
  end

  if (GroundedMask<-0.5) then
    P_effe = 0.01
  else
  P_ice = -rhoi*gravity*thickness
  P_water = -rhoo*gravity*(zsl-BedElevation)
  P_effe = P_ice - P_water
  end

  if (P_effe < 0.01) then
      P_effe = 0.01
  end

  return P_effe
end

-- function to decide the groundedmask based on initial geometry
function mask_ini(surface,thickness,BedElevation)
  watercolumn=surface-thickness-BedElevation
  if (watercolumn<0) then
    mask_initial = 1
  elseif (watercolumn>0) then
    mask_initial = -1
  else
    mask_initial = 0
  end
  return mask_initial
end

function rCoulomb_C_ini(beta,SSAVELOCITY1,SSAVELOCITY2,zb,h,GroundedMask)
  SSAvel_mag = math.sqrt(SSAVELOCITY1^2 + SSAVELOCITY2^2)
  if (GroundedMask<-0.5) then
    haf = 0.001
  elseif (zb >= 0 ) then
    haf = h
  else
    haf = h + zb/(rhoi_si/rhoo_si)
  end
  hth = 75.0
  --if haf < 0.001 then
  --  haf = 0.001
  --end
  lbd = haf/hth
  if (lbd > 1) then 
    lbd = 1.0
  elseif (lbd < 0) then 
    lbd = 1.0
  end
  tau_b = 10.0^beta*SSAvel_mag
  m = 3.0
  u0 = 300
  C = (tau_b/((SSAvel_mag/(SSAvel_mag+u0))^(1/m)))*(1.0/lbd)
  --if (GroundedMask<-0.5) then
  --  C = 0.1
  --end
  return C
end


function rCoulomb_C_NoScaling(beta,SSAVELOCITY1,SSAVELOCITY2)
  SSAvel_mag = math.sqrt(SSAVELOCITY1^2 + SSAVELOCITY2^2)
  tau_b = 10.0^beta*SSAvel_mag
  m = 3.0
  u0 = 300
  C = (tau_b/((SSAvel_mag/(SSAvel_mag+u0))^(1/m)))
  return C
end



function HAF_zc(zb,h,GroundedMask)
  if (GroundedMask<-0.5) then
    haf = 0.001
  elseif (zb >= 0 ) then
    haf = h
  else
    haf = h + zb/(rhoi_si/rhoo_si)
  end
  return haf
end

-- function to deepen the bedrock 
-- minimum watercolumn thickness is 50 m
-- 200m deepening tapering to 100m over 100km
function bed_deepen(surface,thickness,BedElevation,GroundedMask,DistGL)
  --   if (thickness < MINH) then -- this is already done in floatLower or upper.
  --     thickness = MINH
  --   end
    if (BedElevation == nan) then
      BedElevation = 0.0
    end
    watercolumn = surface-thickness-BedElevation
    
    if ( (watercolumn >= 0) and (watercolumn <50)) then 
       watercolumn_new= 50 + (100-50)/100000*DistGL
    elseif (watercolumn < 0) then
      watercolumn_new = 50
    else
       watercolumn_new = watercolumn
    end
    bed_new = BedElevation
    if (GroundedMask > -0.5 and watercolumn < 0) then
      bed_new = surface-thickness
    end
    if (GroundedMask<-0.5) then
      bed_new = surface-thickness - watercolumn_new
    end
    return bed_new
end
-- friction loads
function max(a,b)
  if (a>b) then
    return a
  end
  return b
end

function frictionloads(lx,ly,lz,vx,vy,vz)
  cx = max(-lx*vx,0)
  cy = max(-ly*vy,0)
  cz = max(-lz*vz,0)
  return cx + cy + cz
end

function frictionloads_mask(lx,ly,lz,vx,vy,vz,groundedmask)
--  print("lx:", lx, "ly:", ly, "lz:", lz)
--  print("vx:", vx, "vy:", vy, "vz:", vz)
--  print("groundedmask:", groundedmask)
  -- Ensure variables are not nil
  lx = lx or 0
  ly = ly or 0
  lz = lz or 0
  vx = vx or 0
  vy = vy or 0
  vz = vz or 0
  groundedmask = groundedmask or 0

  cx = max(-lx*vx,0)
  cy = max(-ly*vy,0)
  cz = max(-lz*vz,0)
--  if (GroundedMask<-0.5) then
  if (groundedmask<0.5) then
    cx = 0
    cy = 0
    cz = 0
  end
  return cx + cy + cz
end

function bottomheatflux(ghf,groundedmask,pressure)
   ustar = 0.005
   -- thermal exchange coefficient
   gammaT = 0.006
   Tocean = -1.5
   -- Tice = pressuremelting_salinity(pressure)
   Tice = 272.5
   dT = Tice - 273.15 - Tocean
    -- 3974
    --
   if (groundedmask<-0.5) then
      -- gammaT*rhoo_si*cw_SI*ustar*dT unit: Pa m yr-1
      -- unit of heatflux in Elmer/ice MPa m yr-1
      heatflux = gammaT*rhoo_si*cw_si*ustar*dT*Pa2MPa
   else
      heatflux = ghf*yearinsec*Pa2MPa
   end
   return heatflux
end


function bottomheatflux_roms(ghf,groundedmask,bmb_flux)
   if (groundedmask<-0.5) then
      -- gammaT*rhoo_si*cw_SI*ustar*dT unit: Pa m yr-1
      -- unit of heatflux in Elmer/ice MPa m yr-1
      heatflux = bmb_flux
   else
      heatflux = ghf*yearinsec*Pa2MPa
   end
   return heatflux
end