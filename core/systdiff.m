%% System of Equations for siger_e
if str2double(bnd(1))==1 
    sige_e   = ( sige_pe*dme   - sige_*dme   + sige_*dpe   - sige_me*dpe   )/(2*dpe*dme);
    sige_ee  = ( sige_epe*dme  - sige_e*dme  + sige_e*dpe  - sige_eme*dpe  )/(2*dpe*dme);
    sige_eee = ( sige_eepe*dme - sige_ee*dme + sige_ee*dpe - sige_eeme*dpe )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    sige_e   = sige_epe - sige_eepe*dpe;
    sige_ee  = sige_eepe - sige_eeepe*dpe;
    sige_eee = ( sige_eepe - sige_ee )/dpe;
 
elseif str2double(bnd(1))==2
    sige_e   = sige_eme + sige_eeme*dme;
    sige_ee  = sige_eeme + sige_eeeme*dme;
    sige_eee = ( sige_ee - sige_eeme )/dme;

end
 
if str2double(bnd(2))==1       
    sige_z   = ( sige_pz*dmz   - sige_*dmz   + sige_*dpz   - sige_mz*dpz  )/(2*dpz*dmz);
    sige_zz  = ( sige_zpz*dmz  - sige_z*dmz  + sige_z*dpz  - sige_zmz*dpz )/(2*dpz*dmz);
    sige_zzz = ( sige_zzpz*dmz - sige_zz*dmz + sige_zz*dpz - sige_zzmz*dpz )/(2*dpz*dmz);
     
elseif str2double(bnd(2))==0
    sige_z   = sige_zpz - sige_zzpz*dpz;
    sige_zz  = sige_zzpz - sige_zzzpz*dpz;
    sige_zzz = ( sige_zzpz - sige_zz )/dpz;
    
elseif str2double(bnd(2))==2
    sige_z   = sige_zmz  + sige_zzmz*dmz;
    sige_zz  = sige_zzmz + sige_zzzmz*dmz;
    sige_zzz = ( sige_zz - sige_zzmz )/dmz;
 
end
 
if str2double(bnd(1))==1    
    sige_ze_ = ( sige_zpe*dme - sige_z*dme + sige_z*dpe - sige_zme*dpe)/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    sige_ze_ = sige_ezpe - sige_ezepe*dpe;
    
elseif str2double(bnd(1))==2
    sige_ze_ = sige_ezme + sige_ezeme*dme;
    
end
 
if str2double(bnd(2))==1        
    sige_ez_ = ( sige_epz*dmz - sige_e*dmz + sige_e*dpz - sige_emz*dpz)/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    sige_ez_  = sige_ezpz - sige_ezzpz*dpz;
 
elseif str2double(bnd(2))==2
    sige_ez_  = sige_ezmz + sige_ezzmz*dmz;
 
end
 
if and(str2double(bnd(1))==3,str2double(bnd(2))==3)
    sige_e   = 0;
    sige_ee  = 0;
    sige_eee = 0;
    
    sige_z   = 0;
    sige_zz  = 0;
    sige_zzz = 0;
    
    sige_ez  = 0;
    sige_eze = 0;
    sige_ezz = 0;
 
else
    sige_ez = 1/2*sige_ez_ + 1/2*sige_ze_;
end
 
if str2double(bnd(1))==1
    sige_eze  = ( sige_ezpe*dme  - sige_ez*dme  + sige_ez*dpe  - sige_ezme*dpe  )/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    sige_eze  = ( sige_ezpe  - sige_ez  )/dpe;
 
elseif str2double(bnd(1))==2
    sige_eze  = ( sige_ez  - sige_ezme  )/dme;
end
 
if str2double(bnd(2))==1
    sige_ezz  = ( sige_ezpz*dmz  - sige_ez*dmz  + sige_ez*dpz  - sige_ezmz*dpz  )/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    sige_ezz  = ( sige_ezpz  - sige_ez  )/dpz;
 
elseif str2double(bnd(2))==2
    sige_ezz  = ( sige_ez  - sige_ezmz  )/dmz;
end

%% System of Equations for sigez_z
if str2double(bnd(1))==1 
    sigz_e   = ( sigz_pe*dme   - sigz_*dme   + sigz_*dpe   - sigz_me*dpe   )/(2*dpe*dme);
    sigz_ee  = ( sigz_epe*dme  - sigz_e*dme  + sigz_e*dpe  - sigz_eme*dpe  )/(2*dpe*dme);
    sigz_eee = ( sigz_eepe*dme - sigz_ee*dme + sigz_ee*dpe - sigz_eeme*dpe )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    sigz_e   = sigz_epe - sigz_eepe*dpe;
    sigz_ee  = sigz_eepe - sigz_eeepe*dpe;
    sigz_eee = ( sigz_eepe - sigz_ee )/dpe;
 
elseif str2double(bnd(1))==2
    sigz_e   = sigz_eme + sigz_eeme*dme;
    sigz_ee  = sigz_eeme + sigz_eeeme*dme;
    sigz_eee = ( sigz_ee - sigz_eeme )/dme;

end
 
if str2double(bnd(2))==1       
    sigz_z   = ( sigz_pz*dmz   - sigz_*dmz   + sigz_*dpz   - sigz_mz*dpz  )/(2*dpz*dmz);
    sigz_zz  = ( sigz_zpz*dmz  - sigz_z*dmz  + sigz_z*dpz  - sigz_zmz*dpz )/(2*dpz*dmz);
    sigz_zzz = ( sigz_zzpz*dmz - sigz_zz*dmz + sigz_zz*dpz - sigz_zzmz*dpz )/(2*dpz*dmz);
     
elseif str2double(bnd(2))==0
    sigz_z   = sigz_zpz - sigz_zzpz*dpz;
    sigz_zz  = sigz_zzpz - sigz_zzzpz*dpz;
    sigz_zzz = ( sigz_zzpz - sigz_zz )/dpz;
    
elseif str2double(bnd(2))==2
    sigz_z   = sigz_zmz  + sigz_zzmz*dmz;
    sigz_zz  = sigz_zzmz + sigz_zzzmz*dmz;
    sigz_zzz = ( sigz_zz - sigz_zzmz )/dmz;
 
end
 
if str2double(bnd(1))==1    
    sigz_ze_ = ( sigz_zpe*dme - sigz_z*dme + sigz_z*dpe - sigz_zme*dpe)/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    sigz_ze_ = sigz_ezpe - sigz_ezepe*dpe;
    
elseif str2double(bnd(1))==2
    sigz_ze_ = sigz_ezme + sigz_ezeme*dme;
    
end
 
if str2double(bnd(2))==1        
    sigz_ez_ = ( sigz_epz*dmz - sigz_e*dmz + sigz_e*dpz - sigz_emz*dpz)/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    sigz_ez_  = sigz_ezpz - sigz_ezzpz*dpz;
 
elseif str2double(bnd(2))==2
    sigz_ez_  = sigz_ezmz + sigz_ezzmz*dmz;
 
end
 
if and(str2double(bnd(1))==3,str2double(bnd(2))==3)
    sigz_e   = 0;
    sigz_ee  = 0;
    sigz_eee = 0;
    
    sigz_z   = 0;
    sigz_zz  = 0;
    sigz_zzz = 0;
    
    sigz_ez  = 0;
    sigz_eze = 0;
    sigz_ezz = 0;
 
else
    sigz_ez = 1/2*sigz_ez_ + 1/2*sigz_ze_;
end
 
if str2double(bnd(1))==1
    sigz_eze  = ( sigz_ezpe*dme  - sigz_ez*dme  + sigz_ez*dpe  - sigz_ezme*dpe  )/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    sigz_eze  = ( sigz_ezpe  - sigz_ez  )/dpe;
 
elseif str2double(bnd(1))==2
    sigz_eze  = ( sigz_ez  - sigz_ezme  )/dme;
end
 
if str2double(bnd(2))==1
    sigz_ezz  = ( sigz_ezpz*dmz  - sigz_ez*dmz  + sigz_ez*dpz  - sigz_ezmz*dpz  )/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    sigz_ezz  = ( sigz_ezpz  - sigz_ez  )/dpz;
 
elseif str2double(bnd(2))==2
    sigz_ezz  = ( sigz_ez  - sigz_ezmz  )/dmz;
end

%% System of Equations for sigez
if str2double(bnd(1))==1 
    sigez_e   = ( sigez_pe*dme   - sigez_*dme   + sigez_*dpe   - sigez_me*dpe   )/(2*dpe*dme);
    sigez_ee  = ( sigez_epe*dme  - sigez_e*dme  + sigez_e*dpe  - sigez_eme*dpe  )/(2*dpe*dme);
    sigez_eee = ( sigez_eepe*dme - sigez_ee*dme + sigez_ee*dpe - sigez_eeme*dpe )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    sigez_e   = sigez_epe - sigez_eepe*dpe;
    sigez_ee  = sigez_eepe - sigez_eeepe*dpe;
    sigez_eee = ( sigez_eepe - sigez_ee )/dpe;
 
elseif str2double(bnd(1))==2
    sigez_e   = sigez_eme + sigez_eeme*dme;
    sigez_ee  = sigez_eeme + sigez_eeeme*dme;
    sigez_eee = ( sigez_ee - sigez_eeme )/dme;

end
 
if str2double(bnd(2))==1       
    sigez_z   = ( sigez_pz*dmz   - sigez_*dmz   + sigez_*dpz   - sigez_mz*dpz  )/(2*dpz*dmz);
    sigez_zz  = ( sigez_zpz*dmz  - sigez_z*dmz  + sigez_z*dpz  - sigez_zmz*dpz )/(2*dpz*dmz);
    sigez_zzz = ( sigez_zzpz*dmz - sigez_zz*dmz + sigez_zz*dpz - sigez_zzmz*dpz )/(2*dpz*dmz);
     
elseif str2double(bnd(2))==0
    sigez_z   = sigez_zpz - sigez_zzpz*dpz;
    sigez_zz  = sigez_zzpz - sigez_zzzpz*dpz;
    sigez_zzz = ( sigez_zzpz - sigez_zz )/dpz;
    
elseif str2double(bnd(2))==2
    sigez_z   = sigez_zmz  + sigez_zzmz*dmz;
    sigez_zz  = sigez_zzmz + sigez_zzzmz*dmz;
    sigez_zzz = ( sigez_zz - sigez_zzmz )/dmz;
 
end
 
if str2double(bnd(1))==1    
    sigez_ze_ = ( sigez_zpe*dme - sigez_z*dme + sigez_z*dpe - sigez_zme*dpe)/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    sigez_ze_ = sigez_ezpe - sigez_ezepe*dpe;
    
elseif str2double(bnd(1))==2
    sigez_ze_ = sigez_ezme + sigez_ezeme*dme;
    
end
 
if str2double(bnd(2))==1        
    sigez_ez_ = ( sigez_epz*dmz - sigez_e*dmz + sigez_e*dpz - sigez_emz*dpz)/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    sigez_ez_  = sigez_ezpz - sigez_ezzpz*dpz;
 
elseif str2double(bnd(2))==2
    sigez_ez_  = sigez_ezmz + sigez_ezzmz*dmz;
 
end
 
if and(str2double(bnd(1))==3,str2double(bnd(2))==3)
    sigez_e   = 0;
    sigez_ee  = 0;
    sigez_eee = 0;
    
    sigez_z   = 0;
    sigez_zz  = 0;
    sigez_zzz = 0;
    
    sigez_ez  = 0;
    sigez_eze = 0;
    sigez_ezz = 0;
 
else
    sigez_ez = 1/2*sigez_ez_ + 1/2*sigez_ze_;
end
 
if str2double(bnd(1))==1
    sigez_eze  = ( sigez_ezpe*dme  - sigez_ez*dme  + sigez_ez*dpe  - sigez_ezme*dpe  )/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    sigez_eze  = ( sigez_ezpe  - sigez_ez  )/dpe;
 
elseif str2double(bnd(1))==2
    sigez_eze  = ( sigez_ez  - sigez_ezme  )/dme;
end
 
if str2double(bnd(2))==1
    sigez_ezz  = ( sigez_ezpz*dmz  - sigez_ez*dmz  + sigez_ez*dpz  - sigez_ezmz*dpz  )/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    sigez_ezz  = ( sigez_ezpz  - sigez_ez  )/dpz;
 
elseif str2double(bnd(2))==2
    sigez_ezz  = ( sigez_ez  - sigez_ezmz  )/dmz;
end


%% System of Equations for mue
if str2double(bnd(1))==1 
    mue_e   = ( mue_pe*dme   - mue_*dme   + mue_*dpe   - mue_me*dpe   )/(2*dpe*dme);
    mue_ee  = ( mue_epe*dme  - mue_e*dme  + mue_e*dpe  - mue_eme*dpe  )/(2*dpe*dme);
    mue_eee = ( mue_eepe*dme - mue_ee*dme + mue_ee*dpe - mue_eeme*dpe )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    mue_e   = mue_epe - mue_eepe*dpe;
    mue_ee  = mue_eepe - mue_eeepe*dpe;
    mue_eee = ( mue_eepe - mue_ee )/dpe;
 
elseif str2double(bnd(1))==2
    mue_e   = mue_eme + mue_eeme*dme;
    mue_ee  = mue_eeme + mue_eeeme*dme;
    mue_eee = ( mue_ee - mue_eeme )/dme;

end
 
if str2double(bnd(2))==1       
    mue_z   = ( mue_pz*dmz   - mue_*dmz   + mue_*dpz   - mue_mz*dpz  )/(2*dpz*dmz);
    mue_zz  = ( mue_zpz*dmz  - mue_z*dmz  + mue_z*dpz  - mue_zmz*dpz )/(2*dpz*dmz);
    mue_zzz = ( mue_zzpz*dmz - mue_zz*dmz + mue_zz*dpz - mue_zzmz*dpz )/(2*dpz*dmz);
     
elseif str2double(bnd(2))==0
    mue_z   = mue_zpz - mue_zzpz*dpz;
    mue_zz  = mue_zzpz - mue_zzzpz*dpz;
    mue_zzz = ( mue_zzpz - mue_zz )/dpz;
    
elseif str2double(bnd(2))==2
    mue_z   = mue_zmz  + mue_zzmz*dmz;
    mue_zz  = mue_zzmz + mue_zzzmz*dmz;
    mue_zzz = ( mue_zz - mue_zzmz )/dmz;
 
end
 
if str2double(bnd(1))==1    
    mue_ze_ = ( mue_zpe*dme - mue_z*dme + mue_z*dpe - mue_zme*dpe)/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    mue_ze_ = mue_ezpe - mue_ezepe*dpe;
    
elseif str2double(bnd(1))==2
    mue_ze_ = mue_ezme + mue_ezeme*dme;
    
end
 
if str2double(bnd(2))==1        
    mue_ez_ = ( mue_epz*dmz - mue_e*dmz + mue_e*dpz - mue_emz*dpz)/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    mue_ez_  = mue_ezpz - mue_ezzpz*dpz;
 
elseif str2double(bnd(2))==2
    mue_ez_  = mue_ezmz + mue_ezzmz*dmz;
 
end
 
if and(str2double(bnd(1))==3,str2double(bnd(2))==3)
    mue_e   = 0;
    mue_ee  = 0;
    mue_eee = 0;
    
    mue_z   = 0;
    mue_zz  = 0;
    mue_zzz = 0;
    
    mue_ez  = 0;
    mue_eze = 0;
    mue_ezz = 0;
 
else
    mue_ez = 1/2*mue_ez_ + 1/2*mue_ze_;
end
 
if str2double(bnd(1))==1
    mue_eze  = ( mue_ezpe*dme  - mue_ez*dme  + mue_ez*dpe  - mue_ezme*dpe  )/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    mue_eze  = ( mue_ezpe  - mue_ez  )/dpe;
 
elseif str2double(bnd(1))==2
    mue_eze  = ( mue_ez  - mue_ezme  )/dme;
end
 
if str2double(bnd(2))==1
    mue_ezz  = ( mue_ezpz*dmz  - mue_ez*dmz  + mue_ez*dpz  - mue_ezmz*dpz  )/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    mue_ezz  = ( mue_ezpz  - mue_ez  )/dpz;
 
elseif str2double(bnd(2))==2
    mue_ezz  = ( mue_ez  - mue_ezmz  )/dmz;
end

%% System of Equations for mues_z
if str2double(bnd(1))==1 
    muz_e   = ( muz_pe*dme   - muz_*dme   + muz_*dpe   - muz_me*dpe   )/(2*dpe*dme);
    muz_ee  = ( muz_epe*dme  - muz_e*dme  + muz_e*dpe  - muz_eme*dpe  )/(2*dpe*dme);
    muz_eee = ( muz_eepe*dme - muz_ee*dme + muz_ee*dpe - muz_eeme*dpe )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    muz_e   = muz_epe - muz_eepe*dpe;
    muz_ee  = muz_eepe - muz_eeepe*dpe;
    muz_eee = ( muz_eepe - muz_ee )/dpe;
 
elseif str2double(bnd(1))==2
    muz_e   = muz_eme + muz_eeme*dme;
    muz_ee  = muz_eeme + muz_eeeme*dme;
    muz_eee = ( muz_ee - muz_eeme )/dme;

end
 
if str2double(bnd(2))==1       
    muz_z   = ( muz_pz*dmz   - muz_*dmz   + muz_*dpz   - muz_mz*dpz  )/(2*dpz*dmz);
    muz_zz  = ( muz_zpz*dmz  - muz_z*dmz  + muz_z*dpz  - muz_zmz*dpz )/(2*dpz*dmz);
    muz_zzz = ( muz_zzpz*dmz - muz_zz*dmz + muz_zz*dpz - muz_zzmz*dpz )/(2*dpz*dmz);
     
elseif str2double(bnd(2))==0
    muz_z   = muz_zpz - muz_zzpz*dpz;
    muz_zz  = muz_zzpz - muz_zzzpz*dpz;
    muz_zzz = ( muz_zzpz - muz_zz )/dpz;
    
elseif str2double(bnd(2))==2
    muz_z   = muz_zmz  + muz_zzmz*dmz;
    muz_zz  = muz_zzmz + muz_zzzmz*dmz;
    muz_zzz = ( muz_zz - muz_zzmz )/dmz;
 
end
 
if str2double(bnd(1))==1    
    muz_ze_ = ( muz_zpe*dme - muz_z*dme + muz_z*dpe - muz_zme*dpe)/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    muz_ze_ = muz_ezpe - muz_ezepe*dpe;
    
elseif str2double(bnd(1))==2
    muz_ze_ = muz_ezme + muz_ezeme*dme;
    
end
 
if str2double(bnd(2))==1        
    muz_ez_ = ( muz_epz*dmz - muz_e*dmz + muz_e*dpz - muz_emz*dpz)/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    muz_ez_  = muz_ezpz - muz_ezzpz*dpz;
 
elseif str2double(bnd(2))==2
    muz_ez_  = muz_ezmz + muz_ezzmz*dmz;
 
end
 
if and(str2double(bnd(1))==3,str2double(bnd(2))==3)
    muz_e   = 0;
    muz_ee  = 0;
    muz_eee = 0;
    
    muz_z   = 0;
    muz_zz  = 0;
    muz_zzz = 0;
    
    muz_ez  = 0;
    muz_eze = 0;
    muz_ezz = 0;
 
else
    muz_ez = 1/2*muz_ez_ + 1/2*muz_ze_;
end
 
if str2double(bnd(1))==1
    muz_eze  = ( muz_ezpe*dme  - muz_ez*dme  + muz_ez*dpe  - muz_ezme*dpe  )/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    muz_eze  = ( muz_ezpe  - muz_ez  )/dpe;
 
elseif str2double(bnd(1))==2
    muz_eze  = ( muz_ez  - muz_ezme  )/dme;
end
 
if str2double(bnd(2))==1
    muz_ezz  = ( muz_ezpz*dmz  - muz_ez*dmz  + muz_ez*dpz  - muz_ezmz*dpz  )/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    muz_ezz  = ( muz_ezpz  - muz_ez  )/dpz;
 
elseif str2double(bnd(2))==2
    muz_ezz  = ( muz_ez  - muz_ezmz  )/dmz;
end