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

const double L=DSH*0; // extension of nozzle beyond tabulated depth
const bool only=false; // drill-only, no bulk ice, no hose simulation
const double hour=3600;

struct segment{
  double hose_T;
  double hflo_H; // heat flow through hose (into the water)
  double hole_T;
  double hole_R;
  double ice_y[IRN]; // T/sqrt(r)
};

ostream& operator<<(ostream& out, const segment& a){
  double T0=a.ice_y[0];
  out<<a.hose_T<<" "<<(T0<0?T0:a.hole_T)<<" "<<a.hole_R;
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
  const int tti=-1;   // depth slice for bulk ice printout
  double ttl=0;  // time of the last printout

  double A;  // m^2 cross-sectional area of the flow region between the hose and hole wall
  double Rh, Dh; // m hydraulic radius and diameter (4A/p)
  double mu; // kg/(m s) dynamic viscosity of water
  double h;  // W/(m^2 C) convective heat transfer coefficient
  double Tw; // C bulk water temperature in the hole
  double Re; // Reynolds number
  double Pr; // Prandtl number

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
      // near entrance region, flow it not yet developed
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
    nu*=1.0; // value of 5.0 was fitted to Upgrade strings (or case 6) with 0.1 roughness
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

  const int pz=2; // parameterization: 1: x=1-R/r; y=T/sqrt(r); 2: x=log(r/R)/log(Ry/R); y=T

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
	for(int j=0; j<IRN; j++) s[i].ice_y[j]=j>0?Tinf*(pz==1?sqrt((1-j*dx)/Ri):1):0;
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

	if(R>Rz){
	  double dt1=dt, dt2=0;
	  double dydx=(4*y[1]-3*y[0]-y[2])/(2*dx);
	  double dTdr=pz==1?dydx/sqrt(R)+y[0]/(2*sqrt(R)):dydx/(R*log(Ry/R));
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

	    if(pz==1){
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
	    else if(pz==2){
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

	    if(i==tti && ttt>ttl){
	      ttl=ttt;
	      for(int i=0; i<IRN-1; i++){
		double ri, Ti;
		if(pz==1) ri=RR/(1-i*dx), Ti=y[i]*sqrt(ri);
		else if(pz==2) ri=RR*pow(Ry/R, i*dx), Ti=y[i];
		cout<<(ttt/hour)<<" "<<ri<<" "<<Ti<<endl;
	      }
	    }
	  }
	  if(R<Rz || dt2>0) R=Rz;

	  if(R==Rz){
	    bool verbose=false;

	    if(verbose){
	      cerr<<"A "<<R<<endl;
	      for(int i=0; i<IRN-1; i++){
		double ri, Ti;
		if(pz==1) ri=R/(1-i*dx), Ti=y[i]*sqrt(ri);
		else if(pz==2) ri=R*pow(Ry/R, i*dx), Ti=y[i];
		cerr<<ri<<" "<<Ti<<endl;
	      }
	    }

	    double T[IRN], To=0, ro=0;

	    for(int i=0, j=0; i<IRN; i++){
	      double ri=i*dr;
	      for(; j<IRN; j++){
		double rj=j<IRN-1?pz==1?R/(1-j*dx):R*pow(Ry/R, j*dx):Ry;
		double Tj=j<IRN-1?y[j]*(pz==1?sqrt(rj):1):Tinf;

		if(rj>ri || j==IRN-1){
		  T[i]=(To*(rj-ri)+Tj*(ri-ro))/(rj-ro);
		  break;
		}
		else To=Tj, ro=rj;
	      }
	    }

	    for(int i=0; i<IRN; i++) y[i]=T[i];
	    dt=dt2;

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
	  double q=alf*dt/sq(dr);
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

	  if(i==tti && ttt>ttl){
	    ttl=ttt;
	    for(int i=0; i<IRN-1; i++){
	      double ri=i*dx*Ry, Ti=y[i];
	      cout<<(ttt/hour)<<" "<<ri<<" "<<Ti<<endl;
	    }
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
  }

  bool take(int i){
    return s[i].hole_R>=Rz || i>0 && s[i-1].hole_R>=Rz || i<DSN-1 && s[i+1].hole_R>=Rz;
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
	for(int i=0; i<DSN; i++) if(a.take(i)) cout<<(tl/hour)<<" "<<a.depth(i)<<" "<<a.s[i]<<endl;
	tl+=del;
      }
      t+=dt;
    }
  }

  {
    del=6*hour; dt=10, v=-2.25/60;
    for(; t<100*24*hour; t+=dt){
      if(a.dd<0) v=0, fl=0, rt=0;
      a.run(dt, v, fl, rt); // retrieve the drill
      if(t>tl){
	for(int i=0; i<DSN; i++) if(a.take(i)) cout<<(tl/hour)<<" "<<a.depth(i)<<" "<<a.s[i]<<endl;
	tl+=del;
      }
    }
  }
}
