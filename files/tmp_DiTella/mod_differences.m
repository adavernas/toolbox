% differences for price q
if str2double(bnd(1))==1
    qe   = ( qpe*dme   - q*dme   + q*dpe   - qme*dpe   )/(2*dpe*dme);
    qee  = ( qepe*dme  - qe*dme  + qe*dpe  - qeme*dpe  )/(2*dpe*dme);
    qeee = ( qeepe*dme - qee*dme + qee*dpe - qeeme*dpe )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    qe   = qepe    - qeepe*dpe;
    qee  = qeepe   - qeeepe*dpe;
    qeee = ( qeepe - qee )/dpe;

elseif str2double(bnd(1))==2
    qe   = qeme  + qeeme*dme;
    qee  = qeeme + qeeeme*dme;
    qeee = ( qee - qeeme )/dme;

end

if str2double(bnd(2))==1
    qz   = ( qpz*dmz   - q*dmz   + q*dpz   - qmz*dpz  )/(2*dpz*dmz);
    qzz  = ( qzpz*dmz  - qz*dmz  + qz*dpz  - qzmz*dpz )/(2*dpz*dmz);
    qzzz = ( qzzpz*dmz - qzz*dmz + qzz*dpz - qzzmz*dpz )/(2*dpz*dmz);

elseif str2double(bnd(2))==0
    qz   = qzpz    - qzzpz*dpz;
    qzz  = qzzpz   - qzzzpz*dpz;
    qzzz = ( qzzpz - qzz )/dpz;

elseif str2double(bnd(2))==2
    qz   = qzmz  + qzzmz*dmz;
    qzz  = qzzmz + qzzzmz*dmz;
    qzzz = ( qzz - qzzmz )/dmz;

end

if str2double(bnd(1))==1
    qze_ = ( qzpe*dme - qz*dme + qz*dpe - qzme*dpe)/(2*dpe*dme);

elseif str2double(bnd(1))==0
    qze_ = qezpe - qezepe*dpe;

elseif str2double(bnd(1))==2
    qze_ = qezme + qezeme*dme;

end

if str2double(bnd(2))==1
    qez_ = ( qepz*dmz - qe*dmz + qe*dpz - qemz*dpz)/(2*dpz*dmz);

elseif str2double(bnd(2))==0
    qez_  = qezpz - qezzpz*dpz;

elseif str2double(bnd(2))==2
    qez_  = qezmz + qezzmz*dmz;

end

if and(str2double(bnd(1))==3,str2double(bnd(2))==3)
    qe   = 0;
    qee  = 0;
    qeee = 0;

    qz   = 0;
    qzz  = 0;
    qzzz = 0;

    qez  = 0;
    qeze = 0;
    qezz = 0;

else
    qez = 1/2*qez_ + 1/2*qze_;
end

if str2double(bnd(1))==1
    qeze  = ( qezpe*dme - qez*dme  + qez*dpe  - qezme*dpe  )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    qeze  = ( qezpe - qez  )/dpe;

elseif str2double(bnd(1))==2
    qeze  = ( qez - qezme  )/dme;
end

if str2double(bnd(2))==1
    qezz  = ( qezpz*dmz - qez*dmz  + qez*dpz  - qezmz*dpz  )/(2*dpz*dmz);

elseif str2double(bnd(2))==0
    qezz  = ( qezpz - qez  )/dpz;

elseif str2double(bnd(2))==2
    qezz  = ( qez - qezmz  )/dmz;
end

