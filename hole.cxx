#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;
#define DHI 20  // initial payout depth, m
#define DSI 40  // initial hole depth, m
#define DSH 10  // height of the depth slice, m
#define DSN 300 // number of Depth Slices
#define IRN 200 // Number of Radii in Ice

#define HSL DSH // length of surface hose segment
#define HSN DSN // number of surface hose segments

const double hour=3600;
const double L=DSH*0;  // extension of nozzle beyond tabulated depth
const bool only=false; // drill-only, no bulk ice, no hose simulation
const bool dfrm=true;  // calculate deformations
const bool strh=true;  // update coordinates after stretch
const double deo=0.03; // m, compression buffer
double dex=0; // max deo needed

struct segment{
  double hose_T;
  double hflo_H; // heat flow through hose (into the water)
  double hole_T;
  double hole_R;
  double ice_y[IRN]; // T/sqrt(r)
  double ice_t[IRN], ice_u[IRN], ice_r[IRN], ice_h[IRN], ice_z[IRN]; // temperature (C), strain, stress (radial, hoop, z, Pa)
  double ice_du[IRN], ice_ds[IRN]; // deviatoric stress
  int ice_jr; // current hole wall slice
};

ostream& operator<<(ostream& out, const segment& a){
  double T0=a.ice_y[0];
  int j0=a.ice_jr;
  out<<a.hose_T<<" "<<(T0<0?T0:a.hole_T)<<" "<<a.hole_R<<" "<<a.ice_u[j0]<<" "<<a.ice_r[j0];
  return out;
}

class hole{
public:
  const double pi=3.141592654;
  const double xx=1.e-5;

  const double rho_i=917;  // kg/m^3 density of ice
  const double rho_w=982;  // kg/m^3 density of water at 60 C
  const double delta=rho_i/rho_w; // specific volume water over ice
  const double c_i=1950;   // J/(kg C) specific heat of ice
  const double c_w=4170;   // J/(kg C) specific heat of water
  const double c_f=335000; // J/kg heat of fusion of ice

  const double k_i=2.2;    // W/(m C) thermal conductivity of ice
  const double k_w=0.655;  // W/(m C) thermal conductivity of water at 60C (average in 0...100C)
  const double alf=k_i/(rho_i*c_i); // m^2/s alfa in the thermal flow equation
  const double h_w=50;     // W/(m^2 C) convective heat transfer coefficient for still water

  const double L_hose=2900; // total hose length, m
  const double L_surf=200;  // length of open-air hose, m
  double dT_surf=20;   // templerature difference of the outer-inner hose on the surface
  double dT_spool=10;  // templerature difference of the outer-inner hose on the spool

  const double k_hose=only?0:0.26; // W/(m C) thermal conductivity of hose wall
  const double r_inner=0.03175; // m, inner radius of the hose (2.5 in)
  const double r_core=0.047625; // m, outer radius of the hose and drill (3.75 in)
  const double r_step=0.001;    // r step for drill head calculation
  const double Rz=1.e-5; // smallest hole radius
  const double Ry=10.; // radius of the simulated ice around frozen hole
  const double dx=1./(IRN-1);
  const double dr=Ry/(IRN-1);

  double Tinf=-50;    // far-field temperature of ice (-50C)
  double v=2.25/60;   // m/s drill advance speed
  double Vr=1.262e-2; // m^3/s (200 gal/min) water flow rate
  double Ur=Vr;  // average water flow rate
  double dd=DHI; // depth of drill head
  double tt=0;   // time since start

  const double dtt=hour/2; // time-delta (s) between printouts
  const int tti=150;   // depth slice for bulk ice printout
  double ttl=0;  // time of the last printout

  double A;  // m^2 cross-sectional area of the flow region between the hose and hole wall
  double Rh, Dh; // m hydraulic radius and diameter (4A/p)
  double mu; // kg/(m s) dynamic viscosity of water
  double h;  // W/(m^2 C) convective heat transfer coefficient
  double Tw; // C bulk water temperature in the hole
  double Re; // Reynolds number
  double Pr; // Prandtl number

  const double Yo=9.25e9; // Young's modulus, Pa
  const double si=0.33;   // Poisson ratio
  const double al=5.5e-6; // Linear expansion coefficient, K^-1
  const double be=al*(1+si)/(1-si);
  const double ga=(1-2*si)/(2*(1-si));
  const double Yp=Yo/(1+si), Ym=Yo/(1-si);
  const double dl=si/(1-si);

  double Ttip=80;  // C temperature of water exiting the drill tip
  double Tini=88;  // C temperature of water entering the hose
  int drill=-1; // m depth slice of drill
  int dhole=-1; // m depth slice of hole

  double Ri=0;  // m radius of the hole
  double Ka=0, Ya=0, Ra=r_core+r_step;
  double dD=0, Rs=0, Ks=0;

  double RL;

  segment s[DSN];
  double surf_T[HSN];

  double sq(double x){ // square
    return x*x;
  }

  double ca(double x){ // circle area
    return pi*sq(x);
  }

  double neu(){
    return sq(r_inner)/(2*alf);
  }

  double depth(int i){ // depth of depth slice mid-slice
    return DSI+(i+0.5)*DSH;
  }

  double depth_bottom(int i){ // depth of depth slice bottom
    return DSI+(i+1)*DSH;
  }

  double temp(double x){
    if(only) return Tinf;
    const double data[2][6]={
      // fit to https://icecube.wisc.edu/data-releases/2020/05/south-pole-ice-temperature/
      {0.11022E-05, 0.19325E-02, -0.20837E-01, 0.35285, 2201.3, -63.830},
      // fit to table 3 of doi:10.3189/2014AoG68A033
      {-0.17758E-04, 0.13180E-01, -0.29533E-01, 1.1386, 1361.7, -50.754}
    };
    const double * par = data[0];
    double dx=x-par[4];
    double a=par[1]*par[3];
    double b=(par[1]*par[2]+par[0]*par[3])*dx;
    double c=par[0]*par[2]*dx*dx-1;
    double D=sqrt(abs(b*b-4*a*c));
    double fit=(-b+D)/(2*a)+par[5];
    return fit;
  }

  void chtc(double U, double Tb = 0){ // convective heat transfer coefficient
    mu=1/(27*Tw+500); Dh=2*Rh;
    Pr=1/(0.00493*Tw+0.055);
    Re=fabs(U/A)*rho_w*Dh/mu;

    double nu;
    double mub=1/(27*Tb+500);
    double Ta=(Tw+Tb)/2;
    double mua=1/(27*Ta+500);
    double Pra=1/(0.00493*Ta+0.055);
    double Rea=fabs(U/A)*rho_w*Dh/mua;
    double f8=(1/sq(1.82*log10(Rea)-1.64))/8;
    double DL=Dh/(Dh+RL);

    switch(1){
      // fully developed turbulent flow
    case 1: // Dittus and Boelter
      nu=0.023*pow(Re, 0.8)*pow(Pr, Tw>Tb?0.3:0.4);
      break;
    case 2: // Gnielinski (Pr: 0.5 ... 500; Re: 1.e4 ... 5.e6 and 3000 ... 1.e6)
      nu=(Pr<1.5 ? 0.0214*(pow(Re, 0.8)-100) : 0.012*(pow(Re, 0.87)-280))*pow(Pr, 0.4);
      break;
      // property variations across flow cross-section
    case 3: // Sieder and Tate
      nu=0.027*pow(Re, 0.8)*pow(Pr, 1./3)*pow(mu/mub, 0.14);
      break;
      // near entrance region, flow is not yet developed
    case 4: // Nusselt L/d: 10 ... 400
      nu=0.036*pow(Re, 0.8)*pow(Pr, 1./3)*pow(DL, 0.055);
      break;
    case 5: // Petukhnov (Pr: 0.5 ... 2000; Re: 1.e4 ... 5.e6; mub/mu: 0.8 ... 40)
      nu=f8*Rea*Pra/(1.07+12.7*sqrt(f8)*(pow(Pra, 2./3)-1))*pow(mub/mu, Tw>Tb?0.11:0.25);
      break;
    case 6: // (roughness: 1.e-6 ... 1.e-3; Re: 5000 ...1.e8)
      f8=(1.325/sq(log(1.e-3/3.7)+5.74/pow(Re, 0.9)))/8;
      nu=Re*pow(Pr, 1./3)*f8;
      break;
      // fully developed laminar flow
    case 7: // Hausen
      nu=Re*Pr*DL;
      nu=3.66+0.0668*nu/(1+0.04*pow(nu, 2./3));
      break;
    case 8: // Sieder and Tate (Re*Pr*DL>10)
      nu=1.86*pow(Re*Pr*DL, 1./3)*pow(mub/mu, 0.14);
      break;
    }
    nu*=5.0; // value of 5.0 was fitted to Upgrade strings (or case 6 with 0.1 roughness)
    h=max(h_w, nu*k_w/Dh);
  }

  hole(){
    cerr<<"Neumann stability criterion: dt<"<<neu()<<endl;

    for(int i=0; i<DSN; i++){
      s[i].hose_T=0;
      s[i].hole_T=0;
      s[i].hole_R=0;
    }

    for(int i=0; i<HSN; i++){
      surf_T[i]=Tini;
    }

    {
      char * tmp = getenv("surf");
      if(tmp!=NULL) dT_surf=atof(tmp);

      tmp = getenv("spool");
      if(tmp!=NULL) dT_spool=atof(tmp);

      cerr<<"Surface losses: "<<dT_surf<<" "<<dT_spool<<endl;
    }
  }

  void dri(double dd){
    double c=c_f-c_i*temp(dd);
    for(double dr=r_step, R=r_core+dr, Y=0; Y<dd; R+=dr){
      Tw=(Ttip*Vr-ca(R)*v*delta*c/c_w)/(Vr+delta*ca(R)*v);
      Rh=R-r_core, A=ca(R)-ca(r_core);
      Ur=Vr-v*(ca(R)*(1-delta)-ca(r_core));
      RL=Y;
      chtc(Ur);

      cout<<R<<" "<<Y<<" "<<Tw<<endl;
      double dY=v*rho_i*c/(Tw*h)*dr;
      Y+=dY;
    }
  }

  double rJ(int j, double u[], int k){
    return Ry*sq(j/(double) (IRN-1))-deo+u[max(j, k)];
  }

  int JR(double r, double u[], int k = 0){
    int j=min(IRN-1, (int) floor(sqrt(max(0., r+deo)/Ry)*(IRN-1)));
    if(strh){
      double uj=u[max(j, k)];
      if(uj<0){
	for(; j<IRN-1;) if(rJ(++j, u, k)>r) break;
	j--;
      }
      else if(uj>0){
	for(; j>0; j--) if(rJ(j, u, k)<=r) break;
      }
    }
    return j;
  }

  double RJ(int j, double R, double u[], int k){
    return max(R, Ry*sq(j/(double) (IRN-1))-deo+(strh?u[max(j, k)]:0));
  }

  double fRJ(double r[], int k, double R, double R0, double u[]){
    if(k>=IRN-1) return r[IRN-1];
    double rl=RJ(k, R0, u, k);
    double rh=RJ(k+1, R0, u, k);
    return rh>rl?(r[k]*(rh-R)+r[k+1]*(R-rl))/(rh-rl):r[k];
  }

  void postFreeze(int i){
    double & R = s[i].hole_R;
    double * y = s[i].ice_y;
    double T[IRN], To=0, ro=0;

    for(int i=0, j=0; i<IRN; i++){
      double ri=i*dr;
      for(; j<IRN; j++){
	double rj=j<IRN-1?tf==1?R/(1-j*dx):R*pow(Ry/R, j*dx):Ry;
	double Tj=j<IRN-1?y[j]*(tf==1?sqrt(rj):1):Tinf;

	if(rj>ri || j==IRN-1){
	  T[i]=(To*(rj-ri)+Tj*(ri-ro))/(rj-ro);
	  break;
	}
	else To=Tj, ro=rj;
      }
    }

    for(int i=0; i<IRN; i++) y[i]=T[i];
  }

  void run(double dt, double vv, double fl, double rt){
    v=vv;
    Vr=fl*6.30902e-5;
    Tini=rt;
    int n=1;
    while(dt*v>DSH*n) n*=2;
    if(n>1){
      cerr<<"Drill moved "<<(v*dt)<<" m per time step, which is more than simulated depth slice height of "<<DSH<<endl;
      cerr<<"splitting the time step "<<n<<"-fold"<<endl;
    }
    for(int i=0; i<n; i++) run(dt/n);
  }

  const int tf=2; // parameterization: 1: x=1-R/r; y=T/sqrt(r); 2: x=log(r/R)/log(Ry/R); y=T

  void run(double dt){ // run simulation for time duration dt
    if(DSH<v*dt){
      cerr<<"Drill moved "<<(v*dt)<<" m per time step, which is more than simulated depth slice height of "<<DSH<<endl;
      cerr<<"this is currently not supported"<<endl;
      exit(1);
    }
    Ttip=Tini;
    drill=-1;
    for(int i=0; i<HSN; ++i*HSL<L_hose-dd){
      double hflo=2*pi*k_hose*(i*HSL<L_surf?dT_surf:dT_spool)*HSL/(log(r_core/r_inner)*rho_w*c_w);
      double V=ca(r_inner)*HSL;
      Ttip-=(dt*hflo+(Ttip-surf_T[i])*V)/(dt*Vr+V);
      surf_T[i]=Ttip;
    }

    for(int i=0; i<DSN; i++){
      double & hflo = s[i].hflo_H;
      if(depth_bottom(i)<dd){
	hflo=2*pi*k_hose*(s[i].hose_T-s[i].hole_T)*DSH/(log(r_core/r_inner)*rho_w*c_w);
	if(s[i].hose_T==0) s[i].hose_T=Ttip;
	double V=ca(r_inner)*DSH;
	Ttip-=(dt*hflo+(Ttip-s[i].hose_T)*V)/(dt*Vr+V);
	s[i].hose_T=Ttip;
	drill=i;
      }
      else{
	hflo=0;
	s[i].hose_T=0;
      }
    }

    if(drill>dhole){
      dhole=drill;
      Ya=0, Ra=r_core+r_step; dD=0, Rs=0, Ks=0;
      int i=dhole;
      if(s[i].hole_R==0){
	s[i].hole_R=Ri;
	s[i].hole_T=Ka; if(Ka<=0) cerr<<"Ka="<<Ka<<" Ri="<<Ri<<" tt="<<tt<<" dd="<<dd<<endl;
	Tinf=temp(depth(i));
	for(int j=0; j<IRN; j++){
	  s[i].ice_y[j]=j>0?Tinf*(tf==1?sqrt((1-j*dx)/Ri):1):0;
	  s[i].ice_u[j]=0, s[i].ice_r[j]=0, s[i].ice_h[j]=0, s[i].ice_z[j]=0, s[i].ice_du[i]=0, s[i].ice_ds[i]=0;
	  s[i].ice_t[j]=RJ(j, Ri, s[i].ice_u, 0)>Ri?Tinf:0;
	}
	s[i].ice_jr=JR(Ri, s[i].ice_u);
      }
    }

    tt+=dt; dd+=v*dt; Ur=Vr;
    double ttt=floor(tt/dtt)*dtt;

    if(dhole==drill){ // simulation above the drill head until the next depth bin
      double d=dd-depth_bottom(drill)+L, c=c_f-c_i*temp(dd);

      if(v>0 && d>Ya){
	Ri=0; Ka=0; Tw=Ttip;
	double & R = Ra, & Y = Ya;
	for(double dr=r_step; Y<d; R+=dr){
	  Tw=(Ttip*Vr-ca(R)*v*delta*c/c_w)/(Vr+delta*ca(R)*v);
	  Rh=R-r_core, A=ca(R)-ca(r_core);
	  Ur=Vr-v*(ca(R)*(1-delta)-ca(r_core));
	  RL=Y;
	  chtc(Ur);

	  double dY, Yf;
	  if(Tw>0){
	    const double rate=v*rho_i*c/(Tw*h);
	    dY=rate*dr, Yf=Y+dY;
	    if(Yf>d){
	      Yf=d, dY=d-Y;
	      dr=dY/rate;
	    }
	  }
	  else{
	    Yf=d, dY=d-Y;
	    dr=0;
	  }
	  if(Yf>d-DSH){
	    double dYb=Yf-max(d-DSH, Y); // min(dY, d-Y)
	    Ks+=A*dYb*Tw;
	    Rs+=A*dYb;
	    dD+=dYb;
	  }
	  Y=Yf;
	}
	if(dD>0){
	  Ka=Ks/Rs;
	  Ri=sqrt(Rs/(pi*dD)+sq(r_core));
	  Ttip=Ka;
	}
	if(Ka<0) cerr<<"Ka="<<Ka<<" v="<<v<<" tt="<<tt<<" dd="<<dd<<endl;
      }
    }

    double res[DSN] = {0};
    double tos[DSN] = {0};

    double U=0; // water flow rate between cells: up(positive) or down(negative)
    for(int i=dhole; i>=0; i--){
      double & R = s[i].hole_R;
      Rh=R, A=ca(R); if(i<=drill) Rh-=r_core, A-=ca(r_core);

      double Ua=U*s[U>0?i+1:i].hole_T;
      if(i==drill) U+=Ur, Ua+=Ur*Ttip;

      double Vol=A*DSH+U*dt;
      if(Vol<0) cerr<<"Vol="<<Vol<<" is negative tt="<<tt<<" dd="<<dd<<endl;
      Tw=(s[i].hole_T*A*DSH+Ua*dt)/Vol;
      RL=(dhole+0.5-i)*DSH;
      chtc(U);

      double Qice=h*Tw, Rt=0;
      Tinf=temp(depth(i));

      if(only){
	Rt=Qice/(rho_i*(c_f+(Qice>0?-c_i*Tinf:c_w*Tw)));
	R+=Rt*dt;
	if(R<Rz) R=Rz;
      }
      else{
	double * y = s[i].ice_y;
	double R0=R, dt1=0, dt2=dt;

	if(R>Rz){
	  dt1=dt, dt2=0;

	  double dydx=(4*y[1]-3*y[0]-y[2])/(2*dx);
	  double dTdr=tf==1?dydx/sqrt(R)+y[0]/(2*sqrt(R)):dydx/(R*log(Ry/R));
	  double Qnet=Qice+k_i*dTdr;

	  if(Qnet>0){
	    double a=(rho_i*c_i/2)*dTdr, b=rho_i*c_f, c=Qnet*dt;
	    double D=b*b-4*a*c;
	    double dr=D<0?c/b:(b-sqrt(D))/(2*a);
	    Rt=dr/dt;
	  }
	  else{
	    Rt=Qnet/(rho_i*(c_f+c_w*Tw));
	    double dtx=(sq(Rz)-sq(R))/(2*R*Rt);
	    double aux=Rt*dt/R;
	    if(dt<dtx){
	      if(fabs(aux)>xx) Rt=R*(sqrt(1+2*aux)-1)/dt;
	    }
	    else{
	      Rt=(Rz-R)/dtx;
	      dt1=dtx, dt2=dt-dtx;
	    }
	  }

	  if(dt1>0){
	    // substitution [R ... inf] --> [0, 1], trade dT/dx for y term): x=1-R/r; y=T/sqrt(r)
	    // Crank-Nicolson method (average of both time end points), solution by Thomas algorithm

	    double RR=R+Rt*dt1/2;
	    double A[IRN-1], B[IRN-1];
	    A[0]=0; B[0]=y[0];

	    if(tf==1){
	      for(int j=1; j<IRN-1; j++){
		double xi=1-j*dx, xR=xi/RR, ci=alf*sq(xR);
		double bj=dt1*ci*sq(xi/dx)/2;
		double wj=dt1*Rt*xR/(4*dx);
		double gj=dt1*ci/8-2*bj;

		double bb=bj*(y[j-1]+y[j+1])+(1+gj)*y[j]+wj*(y[j+1]-y[j-1]);
		double D=1-gj-(bj-wj)*A[j-1];
		A[j]=(bj+wj)/D;
		B[j]=(bb+(bj-wj)*B[j-1])/D;
	      }
	    }
	    else if(tf==2){
	      double LL=log(Ry/RR);
	      double q=alf*dt1/sq(LL*dx), p=dt1*Rt/(4*RR*LL*dx);
	      for(int j=1; j<IRN-1; j++){
		double xi=j*dx, ri=RR*exp(LL*xi);
		double a=q/(2*sq(ri))-p*(1-xi); // -A[j]
		double b=1+q/sq(ri); // B[j]
		double c=q/(2*sq(ri))+p*(1-xi); // -C[j]
		double d=b-a*A[j-1];
		A[j]=c/d;
		B[j]=(a*y[j-1]+(2-b)*y[j]+c*y[j+1]+a*B[j-1])/d;
	      }
	    }

	    for(int j=IRN-2; j>0; j--) y[j]=B[j]+A[j]*y[j+1];

	    R+=Rt*dt1;
	  }
	  if(R<Rz || dt2>0) R=Rz;

	  if(R==Rz){
	    bool verbose=false;

	    if(verbose){
	      cerr<<"A "<<R<<endl;
	      for(int i=0; i<IRN-1; i++){
		double ri, Ti;
		if(tf==1) ri=R/(1-i*dx), Ti=y[i]*sqrt(ri);
		else if(tf==2) ri=R*pow(Ry/R, i*dx), Ti=y[i];
		cerr<<ri<<" "<<Ti<<endl;
	      }
	    }

	    postFreeze(i);

	    if(verbose){
	      cerr<<"B "<<R<<endl;
	      for(int i=0; i<IRN; i++){
		double ri=i*dr, Ti=y[i];
		cerr<<ri<<" "<<Ti<<endl;
	      }
	    }
	  }
	}

	if(R==Rz){
	  // y=T
	  double A[IRN-1], B[IRN-1];
	  double q=alf*dt2/sq(dr);
	  // y[IRN-1]=Tinf;

	  {
	    double b=1+2*q;
	    A[0]=2*q/b;
	    B[0]=((1-2*q)*y[0]+2*q*y[1])/b;
	  }

	  for(int j=1; j<IRN-1; j++){
	    double a=q*(1-1./(2*j))/2;
	    double b=1+q;
	    double c=q*(1+1./(2*j))/2;
	    double d=b-a*A[j-1];
	    A[j]=c/d;
	    B[j]=(a*y[j-1]+(2-b)*y[j]+c*y[j+1]+a*B[j-1])/d;
	  }

	  for(int j=IRN-2; j>=0; j--) y[j]=B[j]+A[j]*y[j+1];
	}

	double * u = s[i].ice_u;
	double * r = s[i].ice_r;
	double * h = s[i].ice_h;
	double * z = s[i].ice_z;
	double * du = s[i].ice_du;
	double * ds = s[i].ice_ds;

	if(dfrm){
	  double T[IRN], To=0, ro=0;
	  double * t = s[i].ice_t;

	  int j0=s[i].ice_jr, j1=JR(R, u, j0); s[i].ice_jr=j1;
	  int jx=max(j0, j1);

	  for(int i=j1, j=0; i<IRN; i++){
	    double ri=RJ(i, R, u, j0);
	    for(; j<IRN; j++){
	      double rj, Tj;
	      if(j<IRN-1){
		if(R>Rz){
		  if(tf==1) rj=R/(1-j*dx), Tj=y[j]*sqrt(rj);
		  else if(tf==2) rj=R*pow(Ry/R, j*dx), Tj=y[j];
		}
		else rj=j*dx*Ry, Tj=y[j];
	      }
	      else rj=Ry, Tj=Tinf;

	      if(rj>ri || j==IRN-1){
		T[i]=(To*(rj-ri)+Tj*(ri-ro))/(rj-ro);
		break;
	      }
	      else To=Tj, ro=rj;
	    }
	  }

	  // u,r,z,h are given at exactly j except boundary bin,
	  // where the values at Rj<R are given at R

	  double P0=r[j0], dP=0;

	  if(R>R0){
	    r[j1]=fRJ(r, j1, R, R0, u);
	    z[j1]=fRJ(z, j1, R, R0, u);
	    h[j1]=fRJ(h, j1, R, R0, u);
	    t[j1]=fRJ(t, j1, R, R0, u);
	    u[j1]=fRJ(u, j1, R, R0, u);
	    dP=(r[j1]-P0)*sq(R)/Yp;
	    for(int j=j0; j<j1; j++) u[j]=0, r[j]=0, z[j]=0, h[j]=0, t[j]=0, du[j]=0, ds[j]=0;
	  }
	  else if(R<R0){
	    double u0=u[j0];
	    for(int j=j1; j<j0; j++) u[j]=u0, r[j]=P0, z[j]=P0, h[j]=P0; // t[j]=0;
	    // z[j0]=P0, h[j0]=P0;
	  }

	  for(int j=j1; j<IRN; j++) du[j]=0;

	  double fr=0;
	  while(fr<1){
	    double qv=0;
	    double sum=0, sum_1=0, sum_2=0;
	    double Ri, Ti, ti, pz, prh;
	    double vp[IRN] = {0};
	    double vq[IRN] = {0};
	    double dw[IRN] = {0};

	    if(true){ // Visco-Plastic Bulk Relaxation
	      const double mc=1.e-32; // Mobility coefficient, m^4/(Pa^3*s)

	      double vr[IRN], vm[IRN];
	      for(int j=j1; j<IRN; j++){
		vr[j]=RJ(j, max(R, R0), u, j1);
		vm[j]=(r[j]+z[j]+h[j])/3;
	      }
	      // if(R==Rz) vr[j1]=0;

	      if(true){ // A[j]*P[j-1] + B[j]*P[j] + C[j]*P[j+1] = D[j]
		for(int j=j1; j<IRN; j++){
		  double h1, h2, P0, P1, P2, dPdr;

		  if(j<j0) dPdr=0;
		  else if(j==jx){
		    double R0=vr[j], R1=vr[j+1], R2=vr[j+2]; h1=R1-R0, h2=R2-R1;
		    P0=vm[j], P1=vm[j+1], P2=vm[j+2];
		    if(R>Rz){
		      if(h1<=0) dPdr=(P2-P1)/h2;
		      else if(h2<=0) dPdr=(P1-P0)/h1;
		      else dPdr=(-h2*(2*h1+h2)*P0+sq(h1+h2)*P1-sq(h1)*P2)/(h1*h2*(h1+h2));
		    }
		    else{ // gradient spans first cell, preventing stagnation-lock isolation
		      if(h1>0) dPdr=(P1-P0)/h1;
		      else dPdr=(P2-P0)/(h1+h2);
		    }
		  }
		  else if(j<IRN-1){
		    double R0=vr[j-1], R1=vr[j], R2=vr[j+1]; h1=R1-R0, h2=R2-R1;
		    P0=vm[j-1], P1=vm[j], P2=vm[j+1];
		    if(h1<=0) dPdr=(P2-P1)/h2;
		    else if(h2<=0) dPdr=(P1-P0)/h1;
		    else dPdr=(sq(h1)*(P2-P1)+sq(h2)*(P1-P0))/(h1*h2*(h1+h2));
		  }
		  else{
		    double R0=vr[j], R1=vr[j-1], R2=vr[j-2]; h1=R0-R1, h2=R1-R2;
		    P0=vm[j], P1=vm[j-1], P2=vm[j-2];
		    if(h1<=0) dPdr=(P2-P1)/h2;
		    else if(h2<=0) dPdr=(P1-P0)/h1;
		    else dPdr=(h2*(2*h1+h2)*P0-sq(h1+h2)*P1+sq(h1)*P2)/(h1*h2*(h1+h2));
		  }

		  if(j>=j0 && (h1<=0 || h2<=0)) cerr<<"H Error "<<i<<" "<<j<<" "<<j0<<" "<<j1<<" "<<h1<<" "<<h2<<" "<<dPdr<<endl;
		  double VT=mc*sq(dPdr)*dt;
		  vp[j]=3*Yp*VT;
		}

		double A[IRN]={0};
		double B[IRN]={0};
		double C[IRN]={0};
		double D[IRN]={0};

		for(int j=jx; j<IRN; j++){
		  if(j==jx){
		    if(R>Rz){ // Open Hole Dirichlet
		      A[j]=0, B[j]=1, C[j]=0, D[j]=vm[j];
		    }
		    else{ // Centerline Symmetry Neumann
		      double h = vr[j+1]-vr[j];
		      double fc = 2*vp[j]/sq(h);
		      A[j]=0, B[j]=1+fc, C[j]=-fc, D[j]=vm[j];
		    }
		  }
		  else if(j==IRN-1){ // Outer Boundary Neumann
		    double h = vr[j]-vr[j-1];
		    double fo = 2*vp[j]/sq(h);
		    A[j]=-fo, B[j]=1+fo, C[j]=0, D[j]=vm[j];
		  }
		  else { // Interior Cylindrical Nodes
		    double h1 = vr[j]-vr[j-1];
		    double h2 = vr[j+1]-vr[j];
		    double dn = h1*h2*(h1+h2);
		    double fi = -2*vp[j]/dn;

		    A[j] = fi*h2;
		    B[j] = 1-fi*(h1+h2);
		    C[j] = fi*h1;
		    D[j] = vm[j];
		  }
		}

		double cp[IRN]={0};
		double dp[IRN]={0};

		cp[jx] = C[jx]/B[jx];
		dp[jx] = D[jx]/B[jx];

		for(int j=jx+1; j<IRN; j++){
		  double dn=B[j]-A[j]*cp[j-1];

		  if(fabs(dn) > 0){
		    cp[j] = C[j]/dn;
		    dp[j] = (D[j]-A[j]*dp[j-1])/dn;
		  }
		  else{
		    cp[j] = 0;
		    dp[j] = dp[j-1];
		  }
		}

		vm[IRN-1] = dp[IRN-1];
		for(int j=IRN-2; j>=jx; j--) vm[j]=dp[j]-cp[j]*vm[j+1];
	      }

	      for(int j=j1; j<IRN; j++){
		double h1, h2, P0, P1, P2, dPdr, dPd2;

		if(j<j0) dPdr=0, dPd2=0;
		else if(j==jx){
		  if(R>Rz){
		    double R0=vr[j], R1=vr[j+1], R2=vr[j+2]; h1=R1-R0, h2=R2-R1;
		    P0=vm[j], P1=vm[j+1], P2=vm[j+2];
		    if(h1<=0) dPdr=(P2-P1)/h2;
		    else if(h2<=0) dPdr=(P1-P0)/h1;
		    else dPdr=(-h2*(2*h1+h2)*P0+sq(h1+h2)*P1-sq(h1)*P2)/(h1*h2*(h1+h2));
		    dPd2=2*(h2*P0-(h1+h2)*P1+h1*P2)/(h1*h2*(h1+h2));
		  }
		  else dPdr=0, dPd2=0;
		}
		else if(j<IRN-1){
		  double R0=vr[j-1], R1=vr[j], R2=vr[j+1]; h1=R1-R0, h2=R2-R1;
		  P0=vm[j-1], P1=vm[j], P2=vm[j+1];
		  if(h1<=0) dPdr=(P2-P1)/h2;
		  else if(h2<=0) dPdr=(P1-P0)/h1;
		  else dPdr=(sq(h1)*(P2-P1)+sq(h2)*(P1-P0))/(h1*h2*(h1+h2));
		  dPd2=2*(h1*P2-(h1+h2)*P1+h2*P0)/(h1*h2*(h1+h2));
		}
		else{
		  double R0=vr[j], R1=vr[j-1], R2=vr[j-2]; h1=R0-R1, h2=R1-R2;
		  P0=vm[j], P1=vm[j-1], P2=vm[j-2];
		  if(h1<=0) dPdr=(P2-P1)/h2;
		  else if(h2<=0) dPdr=(P1-P0)/h1;
		  else dPdr=(h2*(2*h1+h2)*P0-sq(h1+h2)*P1+sq(h1)*P2)/(h1*h2*(h1+h2));
		  dPd2=2*(h2*P0-(h1+h2)*P1+h1*P2)/(h1*h2*(h1+h2));
		}

		double VT=mc*sq(dPdr)*dt;
		vp[j]=3*VT*dPd2;
		vq[j]=VT*dPdr/vr[j];
	      }
	    }

	    double qr[IRN], qz[IRN], qh[IRN], fx=0;

	    for(int j=j1; j<IRN; j++){
	      double se=sq(r[j]-z[j])+sq(z[j]-h[j])+sq(h[j]-r[j]); // 2*sigma_e^2
	      double sm=(r[j]+z[j]+h[j])/3;

	      double sigma_e=sqrt(se/2); // True von Mises effective stress
	      ds[j]=sigma_e;

	      if(true){
		const double iys=20e6; // ice yield strength limit for shear/crushing (MPa)
		if(sigma_e>iys){
		  cerr<<"The ice has crushed or reached its plastic flow limit."<<endl;
		  double scale=iys/sigma_e;
		  r[j]=sm+(r[j]-sm)*scale;
		  z[j]=sm+(z[j]-sm)*scale;
		  h[j]=sm+(h[j]-sm)*scale;
		  sigma_e=iys;
		  se=2*sq(sigma_e);
		}
	      }

	      double BT=2.4e-24; // Pa^-3 s^-1, near 0 C
	      { // my fit to temperature dependence
		const double bk=8.617333262e-5, Tz=273.15;
		double kt=1/(bk*(Tz+T[j]));
		BT=exp(T[j]>-10?-56.152+0.46397/pow(kt-41.907, 2.4276):-20.5468-0.805838*kt);
	      }

	      double gl=BT*0.75*se*dt; // Glen Law
	      gl/=(1+Yp*gl); // dt-relaxation

	      qr[j]=gl*(r[j]-sm);
	      qz[j]=gl*(z[j]-sm);
	      qh[j]=gl*(h[j]-sm);

	      double pr=vp[j];
	      double pz=-(vp[j]+vq[j]);
	      double ph=vq[j];

	      double fse=(sq(pr-pz)+sq(pz-ph)+sq(ph-pr))/3;
	      fx=max(fx, dt*BT*Yp*sq(Yp)*fse);
	    }

	    const double mx=0.1;
	    double mr=mx/(mx+fx);

	    for(int j=j1; j<IRN; j++){
	      double pr=vp[j];
	      double pz=-(vp[j]+vq[j]);
	      double ph=vq[j];

	      double dr=qr[j]+mr*pr;
	      double dz=qz[j]+mr*pz;
	      double dh=qh[j]+mr*ph;

	      qr[j]=dr, qz[j]=dz, qh[j]=dh;
	      qv=max(qv, fabs(dr));
	      qv=max(qv, fabs(dz));
	      qv=max(qv, fabs(dh));
	    }

	    {
	      const double qx=1.e-4;
	      qv=qv>qx?qx/qv:1;
	      if(fr+qv>1) qv=1-fr, fr=1;
	      else fr+=qv;
	      if(qv<1) cerr<<"Time step: "<<qv<<" "<<fr<<endl;
	    }

	    double ar[IRN], az[IRN], ah[IRN], at[IRN];

	    for(int j=j1; j<IRN; j++){
	      double Rj=RJ(j, R, u, j1), Tj=T[j], tj=t[j];

	      double dr=qv*qr[j];
	      double dz=qv*qz[j];
	      double dh=qv*qh[j];
	      double drh=dr-dh;

	      if(j>j1){ // new ice sum_1,2 - between ri and rc is stress-free
		double ri=Ri, rc=min(max(R0, Ri), Rj), rj=Rj;
		sum += (Ti*(rj+2*ri)+Tj*(2*rj+ri))*(rj-ri)/6 - ti*(rc+ri)*(rc-ri)/2 - (ti*(rj+2*rc)+tj*(2*rj+rc))*(rj-rc)/6;
		sum_1 += (pz*(rj+2*rc)+dz*(2*rj+rc))*(rj-rc)/6; // +pz*(rc-ri)*(ri+rc)/2 (need zero pz, prh)
		sum_2 += (rj>rc?drh-prh+log(rj/rc)*(prh*rj-drh*rc)/(rj-rc):0); // +prh*log(rc/ri)
	      }

	      Ri=Rj, Ti=Tj, ti=tj, pz=dz, prh=drh;
	      t[j]=T[j];
	      if(fr==qv) if(j==j0 && Rj<R0) z[j]=P0, h[j]=P0, dr=0, dz=0, dh=0;

	      dw[j]=(be*sum+qv*dP)/Rj-ga*(sum_1/Rj+sum_2*Rj);

	      double ez=-Ym*al*(Tj-tj);
	      ar[j]=sum_2;
	      az[j]=dr-dz/dl;
	      ah[j]=dr-dh/dl;
	      at[j]=ez;
	    }
	    for(int j=j1; j<IRN; j++){
	      double Rj=RJ(j, R, u, j1);

	      dw[j]+=sum_2*(ga*(sq(Rj)-sq(R))-(R>Rz?sq(R):0))/Rj;
	      double er=Yp*dw[j]/Rj;

	      r[j]+=(ar[j]-sum_2)*Yp-er;
	      z[j]+=(ar[j]-sum_2+az[j])*Yp*dl+at[j];
	      h[j]+=(ar[j]-sum_2+ah[j])*Yp*dl+at[j]+er;
	      du[j]+=dw[j];
	    }
	    // dP=0;
	  }

	  {
	    double dr=Rt*dt1, dP=0;

	    if(R>Rz){
	      double delta=916.7/999.84; // at 0 C
	      double eta=1/delta-1; // Volumetric expansion factor of freezing water (~9%)
	      double Kw=1.96e9-7.15*min(0., r[j1]); // Bulk modulus of liquid water, Pa (Tait's Equation)
	      dP=-(eta*(R0+dr/2)*dr+R0*du[j1])/(1+Yp/(2*Kw));
	    }

	    res[i]=dP, tos[i]=R0;
	  }
	}
      }

      if(R>Rz){
	double WA=2*pi*R*DSH;
	U-=(1-delta)*WA*Rt;
	double dU=U<0&&i>0?-U*dt:0;
	Tw=Tw*Vol+dt*(s[i].hflo_H-WA*Qice/(rho_w*c_w))-(U<0&&i>0?U*s[i-1].hole_T*dt:0);
	Tw/=Vol+dt*WA*Rt*delta-(U<0&&i>0?U*dt:0);
	if(Tw<0) cerr<<"Tw="<<Tw<<" U="<<U<<" tt="<<tt<<" dd="<<dd<<endl;
	if(Tw<xx) Tw=0;
      }
      else Tw=0, U=0;
      s[i].hole_T=Tw;
    }

    double dP=0, dR=0;
    for(int i=dhole, k=i; i>=0; i--){
      double & R = s[i].hole_R;

      double * y = s[i].ice_y;
      double * u = s[i].ice_u;
      double * r = s[i].ice_r;
      double * h = s[i].ice_h;
      double * z = s[i].ice_z;
      double * du = s[i].ice_du;
      double * ds = s[i].ice_ds;

      int j1=s[i].ice_jr;

      if(dfrm){
	if(i<k) k=i;
	if(i==k) dP=0, dR=0;
	for(; k>=0 && s[k].hole_R>Rz; k--) dP+=res[k], dR+=sq(tos[k]);
	if(k<0) dP=0, dR=0;

	double pp=dR>0?sq(R)*dP/dR:0;

	for(int j=j1; j<IRN; j++){
	  double Rj=RJ(j, R, u, j1);
	  double pu=pp/Rj, er=Yp*pu/Rj;
	  du[j]+=pu;
	  r[j]-=er;
	  h[j]+=er;
	}

	{
	  double vr[IRN];
	  for(int j=j1; j<IRN; j++){
	    vr[j]=RJ(j, 0, u, j1);
	    if(!isfinite(u[j])){
	      cerr<<"Bad error: "<<i<<" "<<k<<" "<<j<<" "<<j1<<" "<<R<<" "<<vr[j]<<" "<<u[j]<<" "<<r[j]<<endl;
	      exit(1);
	    }
	  }
	  for(int j=j1; j<IRN-1; j++)
	    if(vr[j]>=vr[j+1]) cerr<<"Bin order error: "<<i<<" "<<k<<" "<<j<<" "<<j1<<" "<<R<<" "<<vr[j]<<" "<<vr[j+1]<<endl;

	  { // check bin boundaries
	    for(int j=j1; j<IRN; j++){
	      const double q=0.1;
	      double v=du[j], w;
	      int err=0;
	      if(v<0){
		w=j>j1?vr[j-1]-vr[j]:vr[j]-vr[j+1];
		if(v<q*w) err=1;
	      }
	      else{
		w=j<IRN-1?vr[j+1]-vr[j]:vr[j]-vr[j-1];
		if(v>q*w) err=2;
	      }
	      if(err>0){
		cerr<<"Bin boundary error: "<<i<<" "<<j<<" "<<j1<<" "<<err<<" "<<v<<" "<<w<<endl;
		v=q*w; // comment out here to use du[j] directly
	      }
	      u[j]+=v;
	    }
	    if(u[j1]>dex){
	      dex=u[j1];
	      if(dex>deo) cerr<<"Max grid compression buffer exceeded: "<<dex<<" "<<deo<<endl;
	    }
	  }
	}

	if(R>Rz){
	  R+=du[j1];
	  if(R<Rz){
	    R=Rz;
	    postFreeze(i);
	  }
	}
      }

      Tinf=temp(depth(i));

      if(i==tti && ttt>ttl){
	ttl=ttt;
	for(int j=0; j<IRN; j++){
	  double rj, Tj;
	  if(j<IRN-1){
	    if(R>Rz){
	      if(tf==1) rj=R/(1-j*dx), Tj=y[j]*sqrt(rj);
	      else if(tf==2) rj=R*pow(Ry/R, j*dx), Tj=y[j];
	    }
	    else rj=j*dx*Ry, Tj=y[j];
	  }
	  else rj=Ry, Tj=Tinf;
	  if(dfrm){
	    double rr=RJ(j, R, u, s[i].ice_jr);
	    cout<<"S "<<(ttt/hour)<<" "<<rj<<" "<<Tj<<" "<<rr<<" "<<u[j]<<" "<<r[j]<<" "<<h[j]<<" "<<z[j]<<" "<<ds[j]<<endl;
	  }
	  else{
	    cout<<"S "<<(ttt/hour)<<" "<<rj<<" "<<Tj<<endl;
	  }
	}
      }
    }
  }

  bool take(int i){
    return s[i].hole_R>=Rz || i>0 && s[i-1].hole_R>=Rz || i<DSN-1 && s[i+1].hole_R>=Rz;
  }

  double PP(int i, double rr){ // rr = reporting radius
    segment& a = s[i];
    double R=a.hole_R;
    double * r = a.ice_r;
    int j0=a.ice_jr;
    if(R<rr){
      double * u = a.ice_u;
      int j=JR(rr, u, j0);
      return fRJ(r, j, rr, R, u);
    }
    else return r[j0];
  }

  void output(double tl, int i){
    cout<<(tl/hour)<<" "<<depth(i)<<" "<<s[i];
    for(int k=1; k<=10; k++) cout<<" "<<PP(i, 0.01*k);
    cout<<endl;
  }
} a;

int main(int arg_c, char *arg_a[]){
  if(arg_c>1){
    cerr<<"Simultating drill: "<<arg_a[0]<<" [dd] [drill speed] [flow rate] [temperature]"<<endl;
    if(arg_c>2) a.v=atof(arg_a[2])/60;
    if(arg_c>3) a.Vr=atof(arg_a[3])*6.30902e-5;
    if(arg_c>4) a.Ttip=atof(arg_a[4]);
    a.dri(atof(arg_a[1]));
    exit(0);
  }

  string in;
  double t=0, tl=0, del=hour/2;
  double dt=1, v=2.25/60, fl=200, rt=88;
  while(getline(cin, in)){
    if(sscanf(in.c_str(), "%lf %lf %lf %lf", &dt, &v, &fl, &rt)==4){
      a.run(dt, v, fl, rt);
      if(t>tl){
	for(int i=0; i<DSN; i++) if(a.take(i)) a.output(tl, i);
	tl+=del;
      }
      t+=dt;
    }
  }

  if(false){
    del=hour/2; dt=1, v=-2.25/60;
    for(; t<10*24*hour; t+=dt){
      if(a.dd<0) v=0, fl=0, rt=0;
      a.run(dt, v, fl, rt); // retrieve the drill
      if(t>tl){
	for(int i=0; i<DSN; i++) if(a.take(i)) a.output(tl, i);
	tl+=del;
      }
    }
  }

  {
    del=6*hour; dt=10, v=-2.25/60;
    for(; t<100*24*hour; t+=dt){
      if(a.dd<0) v=0, fl=0, rt=0;
      a.run(dt, v, fl, rt); // retrieve the drill
      if(t>tl){
	for(int i=0; i<DSN; i++) if(a.take(i)) a.output(tl, i);
	tl+=del;
      }
    }
  }
}
