#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;
#define DSI 0   // initial depth, m
#define DSH 10  // height of the depth slice, m
#define DSN 265 // number of Depth Slices
#define IRN 200 // Number of Radii in Ice

const double L=DSH*0; // extension of nozzle beyond tabulated depth
const bool only=false; // drill-only, no bulk ice, no hose simulation

struct segment{
  double hose_T;
  double hflo_H; // heat flow through hose (into the water)
  double hole_T;
  double hole_R;
  double ice_y[IRN]; // T/sqrt(r)
};

ostream& operator<<(ostream& out, const segment& a){
  out<<a.hose_T<<" "<<a.hole_T<<" "<<a.hole_R;
  return out;
}

class hole{
public:
  const double pi=3.141592654;
  const double dx=1./(IRN-1);
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

  const double k_hose=only?0:0.26; // W/(m C) thermal conductivity of hose wall
  const double r_inner=0.03175; // m, inner radius of the hose (2.5 in)
  const double r_core=0.047625; // m, outer radius of the hose and drill (3.75 in)
  const double r_step=0.001;    // r step for drill head calculation
  const double Rz=1.e-5; // smallest hole radius

  double Tinf=-50;    // far-field temperature of ice (-50C)
  double v=2.25/60;   // m/s drill advance speed
  double Vr=1.262e-2; // m^3/s (200 gal/min) water flow rate
  double Ur=Vr; // average water flow rate
  double dd=0;  // depth of drill head

  double A;  // m^2 cross-sectional area of the flow region between the hose and hole wall
  double Dh; // m hydrolic diameter (4A/p)
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

  segment s[DSN];

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

  void chtc(double U){ // convective heat transfer coefficient
    mu=1/(27*Tw+500);
    Pr=1/(0.00493*Tw+0.055);
    Re=fabs(U/A)*rho_w*Dh/mu;
    h=max(h_w, (k_w/Dh)*0.023*pow(Re, 0.8)*pow(Pr, 0.3));
  }

  hole(){
    cerr<<"Neumann stability criterion: dt<"<<neu()<<endl;

    for(int i=0; i<DSN; i++){
      s[i].hose_T=0;
      s[i].hole_T=0;
      s[i].hole_R=0;
    }
  }

  void dri(double dd){
    double c=c_f-c_i*temp(dd);
    for(double dr=r_step, R=r_core+dr, Y=0; Y<dd; R+=dr){
      Tw=(Ttip*Vr-ca(R)*v*delta*c/c_w)/(Vr+delta*ca(R)*v);
      Dh=R-r_core, A=ca(R)-ca(r_core);
      Ur=Vr-v*(ca(R)*(1-delta)-ca(r_core));
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
    run(dt);
  }

  void run(double dt){ // run simulation for time duration dt
    if(DSH<v*dt){
      cerr<<"Drill moved "<<(v*dt)<<" m per time step, which is more than simulated depth slice height of "<<DSH<<endl;
      cerr<<"This is currently not supported"<<endl;
      exit(1);
    }
    Ttip=Tini;
    drill=-1;
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
	s[i].hole_T=Ka; if(Ka<=0) cerr<<"Ka="<<Ka<<" Ri="<<Ri<<endl;
	Tinf=temp(depth(i));
	for(int j=0; j<IRN; j++) s[i].ice_y[j]=Tinf*sqrt((1-j*dx)/Ri);
      }
    }

    dd+=v*dt; Ur=Vr;

    if(dhole==drill){ // simulation above the drill head until the next depth bin
      double d=dd-depth_bottom(drill)+L, c=c_f-c_i*temp(dd);

      if(v>0 && d>Ya){
	Ri=0; Ka=0; Tw=Ttip;
	double & R = Ra, & Y = Ya;
	for(double dr=r_step; Y<d; R+=dr){
	  Tw=(Ttip*Vr-ca(R)*v*delta*c/c_w)/(Vr+delta*ca(R)*v);
	  Dh=R-r_core, A=ca(R)-ca(r_core);
	  Ur=Vr-v*(ca(R)*(1-delta)-ca(r_core));
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
	if(Ka<0) cerr<<"Ka="<<Ka<<" v="<<v<<endl;
      }
    }

    double U=0; // water flow rate between cells: up(positive) or down(negative)
    for(int i=dhole; i>=0; i--){
      double & R = s[i].hole_R;
      Dh=R, A=ca(R); if(i<=drill) Dh-=r_core, A-=ca(r_core);

      double Ua=U*s[U>0?i+1:i].hole_T;
      if(i==drill) U+=Ur, Ua+=Ur*Ttip;

      double Vol=A*DSH+U*dt;
      if(Vol<0) cerr<<"Vol="<<Vol<<" is negative"<<endl;
      Tw=(s[i].hole_T*A*DSH+Ua*dt)/Vol;
      chtc(U);

      double Qice=h*Tw, Rt=0;

      if(only){
	Rt=Qice/(rho_i*(c_f+(Qice>0?-c_i*temp(dd):c_w*Tw)));
	R+=Rt*dt;
	if(R<Rz) R=Rz;
      }
      else if(R>Rz){
	// substitution [R ... inf] --> [0, 1], trade dT/dx for y term): x=1-R/r; y=T/sqrt(r)
	// Crank-Nicolson method (average of both time end points), solution by Thomas algorithm
	double A[IRN-1], B[IRN], b[IRN], w[IRN], g[IRN];
	double * y = s[i].ice_y;
	y[0]=0;

	double dydx=(4*y[1]-3*y[0]-y[2])/(2*dx);
	double dTdr=dydx/sqrt(R)+y[0]/(2*sqrt(R));
	double Qnet=Qice+k_i*dTdr;

	if(Qnet>0){
	  double a=(rho_i*c_i/2)*dTdr, b=rho_i*c_f, c=Qnet*dt;
	  double D=b*b-4*a*c;
	  double dr=D<0?c/b:(b-sqrt(D))/(2*a);
	  Rt=dr/dt;
	}
	else{
	  Rt=Qnet/(rho_i*(c_f+c_w*Tw));
	}

	for(int j=0; j<IRN; j++){
	  double xi=1-j*dx, xR=xi/R, ci=alf*sq(xR);
	  b[j]=dt*ci*sq(xi/dx)/2;
	  w[j]=dt*Rt*xR/(4*dx);
	  g[j]=dt*ci/8-2*b[j];
	}
	for(int j=1; j<IRN-1; j++) B[j]=b[j]*(y[j-1]+y[j+1])+(1+g[j])*y[j]+w[j]*(y[j+1]-y[j-1]);
	A[0]=0; B[0]=y[0], B[IRN-1]=y[IRN-1];
	for(int j=1; j<IRN-1; j++){
	  double D=1-g[j]-(b[j]-w[j])*A[j-1];
	  A[j]=(b[j]+w[j])/D;
	  B[j]=(B[j]+(b[j]-w[j])*B[j-1])/D;
	}
	for(int j=IRN-2; j>0; j--) y[j]=B[j]+A[j]*y[j+1];
	R+=Rt*dt;
	if(R<Rz) R=Rz;

	if(false) for(int j=0; j<IRN-1; j++){
	    double xi=1-i*dx, Rx=R/xi;
	    cout<<Rx<<" "<<(y[j]*sqrt(Rx))<<endl;
	  }
      }

      if(R>Rz){
	double WA=2*pi*R*DSH;
	U-=(1-delta)*WA*Rt;
	double dU=U<0&&i>0?-U*dt:0;
	Tw=Tw*Vol+dt*(s[i].hflo_H-WA*Qice/(rho_w*c_w))-(U<0&&i>0?U*s[i-1].hole_T*dt:0);
	Tw/=Vol+dt*WA*Rt*delta-(U<0&&i>0?U*dt:0);
	if(Tw<0) cerr<<"Tw="<<Tw<<" U="<<U<<endl;
	if(Tw<xx) Tw=0;
      }
      else Tw=0, U=0;
      s[i].hole_T=Tw;
    }
  }

  bool take(int i){
    return s[i].hole_R>Rz || i>0 && s[i-1].hole_R>Rz || i<DSN-1 && s[i+1].hole_R>Rz;
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
  const double hour=3600;
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
