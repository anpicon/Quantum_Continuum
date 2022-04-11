
vec2x    b0(NEi,vec1x(NR,complexd(0,0)));   // Arrangement pattern: [NEi][NR]
vec2x    b1(NEi,vec1x(NR,complexd(0,0)));
vec2x    b2(NEi,vec1x(NR,complexd(0,0)));
vec2x    bv(NEi,vec1x(NR,complexd(0,0)));
vec2d    PES(NEi,vec1d(NR,0.));             // Potential energy surfaces: PES[Ei][R]
vec3d    Dip1(NEi,vec2d(NEi,vec1d(NR,0.))); // Dipole transitions pump: Dip[Ei][Ej][R]
vec3d    Dip2(NEi,vec2d(NEi,vec1d(NR,0.))); // Dipole transitions Stokes: Dip[Ei][Ej][R]
vec1d    Gamma(NEi,0.);                     // Core-hole decay: Gamma[Ei]
vec1x    Vt(NR,complexd(0,0));              // Vector to be used in the Runge-Kutta

R.resize(NR);                               // R grid (interpolated) R[iR]
for(int iR=0;iR<R.size();iR++) R[iR]=0.;
