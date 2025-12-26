#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits>

using namespace std;

static constexpr double pi = 3.14159265358979323846;
inline double deg2rad(double x){ return x*pi/180.0; }
inline double rad2deg(double x){ return x*180.0/pi; }
template<class T> inline T clampv(T v, T lo, T hi){ return std::min(hi, std::max(lo, v)); }

struct AtmState{
    double T, p, rho, a, g, mu, nu, k;
};

struct AtmosphereGOST {
    AtmState at(double H) const {
        static const double Hb[] = {
            0, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000,
            9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000
        };
        static const double Tb[] = {
            288.150, 284.900, 281.651, 275.154, 268.659, 262.166, 255.676, 249.187, 242.700, 236.215,
            229.733, 223.252, 216.774, 216.650, 216.650, 216.650, 216.650, 216.650, 216.650, 216.650, 216.650, 216.650
        };
        static const double pb[] = {
            1.01325e5, 9.54613e4, 8.98763e4, 7.95014e4, 7.01212e4, 6.16604e4, 5.40483e4, 4.72167e4, 4.11051e4, 3.56516e4,
            3.08007e4, 2.64999e4, 2.26999e4, 1.93994e4, 1.65796e4, 1.41703e4, 1.21118e4, 1.03528e4, 8.84970e3, 7.56521e3, 6.46747e3, 5.52929e3
        };
        static const double rhob[] = {
            1.22500, 1.16727, 1.11166, 1.00655, 0.909254, 0.819347, 0.736429, 0.660111, 0.590018, 0.526783,
            0.467063, 0.413510, 0.364801, 0.311937, 0.266595, 0.227855, 0.194755, 0.166470, 0.142301, 0.121647, 0.103995, 0.0889097
        };
        static const double gb[] = {
            9.8066, 9.8051, 9.8036, 9.8005, 9.7974, 9.7943, 9.7912, 9.7882, 9.7851, 9.7820,
            9.7789, 9.7759, 9.7728, 9.7697, 9.7667, 9.7636, 9.7605, 9.7575, 9.7544, 9.7513, 9.7483, 9.7452
        };
        static const double He[] = {
            0, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
            10000, 11000, 12000, 13000, 14000, 15000, 16000
        };
        static const double ae[] = {
            340.294, 338.370, 336.435, 332.532, 328.584, 324.589, 320.545, 316.452, 312.306, 308.105, 303.848,
            299.532, 295.154, 295.069, 295.069, 295.069, 295.069, 295.069
        };
        static const double mue[] = {
            1.7894e-5, 1.7737e-5, 1.7579e-5, 1.7260e-5, 1.6938e-5, 1.6612e-5, 1.6282e-5, 1.5949e-5, 1.5612e-5, 1.5271e-5, 1.4926e-5,
            1.4577e-5, 1.4223e-5, 1.4216e-5, 1.4216e-5, 1.4216e-5, 1.4216e-5, 1.4216e-5
        };
        static const double nue[] = {
            1.4607e-5, 1.5195e-5, 1.5813e-5, 1.7147e-5, 1.8628e-5, 2.0275e-5, 2.2110e-5, 2.4162e-5, 2.6461e-5, 2.9044e-5, 3.1957e-5,
            3.5251e-5, 3.8988e-5, 4.5574e-5, 5.333e-5, 6.2391e-5, 7.2995e-5, 8.5397e-5
        };
        static const double ke[] = {
            2.5343e-2, 2.5087e-2, 2.4830e-2, 2.4314e-2, 2.3795e-2, 2.3273e-2, 2.2747e-2, 2.2218e-2, 2.1687e-2, 2.1152e-2, 2.0614e-2,
            2.0072e-2, 1.9528e-2, 1.9518e-2, 1.9518e-2, 1.9518e-2, 1.9518e-2, 1.9518e-2
        };

        auto lerp = [](double x1,double x2,double t){ return (1.0-t)*x1 + t*x2; };
        const int Nb = sizeof(Hb)/sizeof(Hb[0]);
        AtmState s;
        
        if (H <= Hb[0]) s = {Tb[0], pb[0], rhob[0], 0, gb[0], 0, 0, 0};
        else if (H >= Hb[Nb-1]) s = {Tb[Nb-1], pb[Nb-1], rhob[Nb-1], 0, gb[Nb-1], 0, 0, 0};
        else {
            int k=1; while(k<Nb && Hb[k]<H) ++k;
            double t=(H-Hb[k-1])/(Hb[k]-Hb[k-1]);
            s = {lerp(Tb[k-1],Tb[k],t), lerp(pb[k-1],pb[k],t), lerp(rhob[k-1],rhob[k],t), 0, lerp(gb[k-1],gb[k],t), 0, 0, 0};
        }

        const int Ne = sizeof(He)/sizeof(He[0]);
        if (H <= He[0]) { s.a=ae[0]; s.mu=mue[0]; s.nu=nue[0]; s.k=ke[0]; }
        else if (H >= He[Ne-1]) {
            const double gamma=1.4, R=287.05287; s.a = sqrt(gamma*R*s.T);
            s.mu=mue[Ne-1]; s.nu=nue[Ne-1]; s.k=ke[Ne-1];
        } else {
            int k=1; while(k<Ne && He[k]<H) ++k;
            double t=(H-He[k-1])/(He[k]-He[k-1]);
            s.a = lerp(ae[k-1],ae[k],t);
            s.mu= lerp(mue[k-1],mue[k],t);
            s.nu= lerp(nue[k-1],nue[k],t);
            s.k = lerp(ke[k-1],ke[k],t);
        }
        if (s.a<=0){ const double gamma=1.4, R=287.05287; s.a=sqrt(gamma*R*s.T); }
        return s;
    }
};

struct AeroConfig{
    double S=127.0, CL0=0.2, CLa=5.5, CD0=0.022, k=0.045;
    double CLmax=1.6, q_max=1e12;
};

struct Aerodynamics{
    AeroConfig cfg;
    inline double CL(double a) const { return cfg.CL0 + cfg.CLa*a; }
    inline double CD(double a) const { double cl=CL(a); return cfg.CD0 + cfg.k*cl*cl; }
};

struct EngineConfig{
    double Tmax_SL=133000.0, throttle_max=1.10, throttle_min=0.30;
    double sfc0=0.070, sfc_kH=0.0, phi_p_deg=0.0;
};

struct Engine{
    EngineConfig cfg; 
    AtmosphereGOST atm;
    
    double Tmax(double H, double M) const{
        auto s=atm.at(H); 
        double sigma=s.rho/1.225;
        double mach_loss=clampv(1.0 - 0.20*max(0.0, M-0.5), 0.60, 1.0);
        return cfg.Tmax_SL*(0.70+0.30*sigma)*mach_loss;
    }
    
    double sfc(double H,double M,double thr) const{
        return cfg.sfc0 * (1.0 + cfg.sfc_kH*(H/10000.0)) * (1.05 - 0.05*thr);
    }
};

enum class Objective{ MinFuel, MinTime, Weighted };

struct Config{
    double V0 = 330.0/3.6, H0 = 500.0, m0 = 47000.0;
    double Vf = 840.0/3.6, Hf = 7000.0;
    int NV=25, NH=22;
    double Hmin=0.0, Hmax=8000.0, Vmin=70.0, Vmax=240.0;
    double alpha_min_deg=-2.0, alpha_max_deg=12.0, alpha_step_deg=0.1;
    double throttle_min=0.30, throttle_max=1.10, throttle_step=0.05;
    double theta_abs_max_deg=15.0;
    double Mmin=0.15, Mmax=0.82;
    Objective objective=Objective::MinFuel;
    double lambda_mix=1.0, time_as_fuel=0.0;
    AtmosphereGOST atmosphere;
    Aerodynamics aero;
    Engine engine;
};

struct StepRes{
    bool ok=false;
    double dt=0.0, fuel=0.0, Theta=0.0, alpha=0.0, throttle=0.0, P=0.0;
    double CL=0,CD=0,q=0,rho=0,a=0,g=0;
};

static double J_inc(const Config& cfg, double fuel, double dt){
    if(cfg.objective==Objective::MinFuel) return fuel;
    if(cfg.objective==Objective::MinTime) return dt;
    return cfg.lambda_mix*fuel + (1.0-cfg.lambda_mix)*(cfg.time_as_fuel*dt);
}

static StepRes step_A(const Config& cfg, double V, double H, double m, double dV){
    StepRes r; 
    if(fabs(dV)<=1e-9) return r;
    auto atm=cfg.atmosphere.at(H); 
    double rho=atm.rho, a=atm.a, g=atm.g;
    double M=V/max(1e-6,a); 
    if(M<cfg.Mmin||M>cfg.Mmax) return r;
    double S=cfg.aero.cfg.S, phi=deg2rad(cfg.engine.cfg.phi_p_deg);
    double best=1e300;

    for(double ad=cfg.alpha_min_deg; ad<=cfg.alpha_max_deg+1e-9; ad+=cfg.alpha_step_deg){
        double aR=deg2rad(ad);
        double CL=cfg.aero.CL(aR), CD=cfg.aero.CD(aR);
        if(fabs(CL)>cfg.aero.cfg.CLmax) continue;
        double q=0.5*rho*V*V; 
        if(q>cfg.aero.cfg.q_max) continue;
        double Y=q*S*CL, X=q*S*CD;
        
        for(double thr=cfg.throttle_min; thr<=cfg.throttle_max+1e-12; thr+=cfg.throttle_step){
            double Pmax=cfg.engine.Tmax(H,M); 
            double P=thr*Pmax;
            double Fv=P*sin(aR+phi)+Y - m*g; 
            if(fabs(Fv) > 0.02*m*g) continue;
            double dVdt=(P*cos(aR+phi)-X)/m; 
            if((dV>0&&dVdt<=1e-6)||(dV<0&&dVdt>=-1e-6)) continue;
            double dt=fabs(dV)/fabs(dVdt);
            double fuel=P*cfg.engine.sfc(H,M,thr)/3600.0 * dt;
            double cost=J_inc(cfg,fuel,dt);
            if(cost<best){ best=cost; r={true,dt,fuel,0.0,aR,thr,P,CL,CD,q,rho,a,g}; }
        }
    }
    return r;
}

static StepRes step_B(const Config& cfg, double V, double H, double m, double dH){
    StepRes r; 
    if(fabs(dH)<=1e-9) return r;
    auto atm=cfg.atmosphere.at(H); 
    double rho=atm.rho, a=atm.a, g=atm.g;
    double M=V/max(1e-6,a); 
    if(M<cfg.Mmin||M>cfg.Mmax) return r;
    double S=cfg.aero.cfg.S, phi=deg2rad(cfg.engine.cfg.phi_p_deg);
    double best=1e300;

    for(double ad=cfg.alpha_min_deg; ad<=cfg.alpha_max_deg+1e-9; ad+=cfg.alpha_step_deg){
        double aR=deg2rad(ad);
        double CL=cfg.aero.CL(aR), CD=cfg.aero.CD(aR);
        if(fabs(CL)>cfg.aero.cfg.CLmax) continue;
        double q=0.5*rho*V*V; 
        if(q>cfg.aero.cfg.q_max) continue;
        double Y=q*S*CL, X=q*S*CD;

        for(double thr=cfg.throttle_min; thr<=cfg.throttle_max+1e-12; thr+=cfg.throttle_step){
            double Pmax=cfg.engine.Tmax(H,M); 
            double P=thr*Pmax;
            double sinT=(P*cos(aR+phi)-X)/(m*g);
            double cosT=(P*sin(aR+phi)+Y)/(m*g);
            double norm=sinT*sinT+cosT*cosT; 
            if(norm<0.95||norm>1.05) continue;
            double Theta=atan2(sinT,cosT);
            if(fabs(rad2deg(Theta))>cfg.theta_abs_max_deg) continue;
            if((dH>0 && sinT<=1e-5) || (dH<0 && sinT>=-1e-5)) continue;

            double rateH=V*sinT; 
            double dt=fabs(dH)/fabs(rateH);
            double fuel=P*cfg.engine.sfc(H,M,thr)/3600.0 * dt;
            double cost=J_inc(cfg,fuel,dt);
            if(cost<best){ best=cost; r={true,dt,fuel,Theta,aR,thr,P,CL,CD,q,rho,a,g}; }
        }
    }
    return r;
}

static StepRes step_C(const Config& cfg, double V, double H, double m, double dV, double dH){
    StepRes r;
    auto atm=cfg.atmosphere.at(H); 
    double rho=atm.rho, a=atm.a, g=atm.g;
    double M=V/max(1e-6,a); 
    if(M<cfg.Mmin||M>cfg.Mmax) return r;
    double S=cfg.aero.cfg.S, phi=deg2rad(cfg.engine.cfg.phi_p_deg);
    double targetK = (fabs(dV)<1e-9 ? 1e300 : (dH/dV));
    double best=1e300;

    for(double ad=cfg.alpha_min_deg; ad<=cfg.alpha_max_deg+1e-9; ad+=cfg.alpha_step_deg){
        double aR=deg2rad(ad);
        double CL=cfg.aero.CL(aR), CD=cfg.aero.CD(aR);
        if(fabs(CL)>cfg.aero.cfg.CLmax) continue;
        double q=0.5*rho*V*V; 
        if(q>cfg.aero.cfg.q_max) continue;
        double Y=q*S*CL, X=q*S*CD;

        for(double thr=cfg.throttle_min; thr<=cfg.throttle_max+1e-12; thr+=cfg.throttle_step){
            double Pmax=cfg.engine.Tmax(H,M); 
            double P=thr*Pmax;
            double cosT=(P*sin(aR+phi)+Y)/(m*g);
            if(cosT<-1.0||cosT>1.0) continue;
            double sinT_abs=sqrt(max(0.0,1.0-cosT*cosT));
            double sinT=(dH>=0 ? +sinT_abs : -sinT_abs);
            double Theta=atan2(sinT,cosT);
            if(fabs(rad2deg(Theta))>cfg.theta_abs_max_deg) continue;

            double dVdt=(P*cos(aR+phi)-X)/m - g*sinT;
            double dHdt=V*sinT;
            if((dV>0 && dVdt<=1e-6) || (dV<0 && dVdt>=-1e-6)) continue;
            if((dH>0 && dHdt<=1e-6) || (dH<0 && dHdt>=-1e-6)) continue;

            double K=dHdt/dVdt; 
            if(!std::isfinite(K)) continue;
            if(fabs((K-targetK)/max(1e-6,fabs(targetK)))>0.10) continue;

            double dt=fabs(dV)/fabs(dVdt);
            double fuel=P*cfg.engine.sfc(H,M,thr)/3600.0 * dt;
            double cost=J_inc(cfg,fuel,dt);
            if(cost<best){ best=cost; r={true,dt,fuel,Theta,aR,thr,P,CL,CD,q,rho,a,g}; }
        }
    }
    return r;
}

struct Node{
    double cost = numeric_limits<double>::infinity();
    double time = 0.0, mass = numeric_limits<double>::quiet_NaN();
    int pi=-1, pj=-1;
    char move='?';
    double H=0, V=0, Theta=0, alpha=0, thr=0, P=0, dt=0;
    double CL=0, CD=0, q=0, rho=0, a=0, g=0;
};

struct EdgeCost{ 
    bool ok=false; 
    double cost=numeric_limits<double>::infinity(), dt=0.0, fuel=0.0; 
};

struct EdgeMatrices{ 
    vector<vector<EdgeCost>> A,B,C; 
    void init(int NV,int NH){ 
        A.assign(NV,vector<EdgeCost>(NH)); 
        B.assign(NV,vector<EdgeCost>(NH)); 
        C.assign(NV,vector<EdgeCost>(NH)); 
    } 
};

struct Result{
    vector<double> V,H,t,m; 
    vector<char> move; 
    vector<double> alpha, thr, Theta, dt, P;
    double total_time=0.0, fuel_used=0.0; 
    bool ok=false;
};

static void saveCSV(const string& path, const Result& R){
    ofstream f(path); 
    f.setf(std::ios::fixed); 
    f<<setprecision(6);
    f<<"i,move,V,H,t,m,alpha_deg,throttle,Theta_deg,dt,P\n";
    for(size_t i=0;i<R.V.size();++i){
        f<<i<<","<<(R.move[i]?string(1,R.move[i]):"")<<","<<R.V[i]<<","<<R.H[i]<<","<<R.t[i]<<","<<R.m[i]<<","
         <<rad2deg(R.alpha[i])<<","<<R.thr[i]<<","<<rad2deg(R.Theta[i])<<","<<R.dt[i]<<","<<R.P[i]<<"\n";
    }
}

struct SolverVH{
    Config cfg; 
    bool debug_last_layer=false; 
    EdgeMatrices edges;

    void dump_pred_last_layer(const vector<double>& Vg,const vector<double>& Hg,
                              int iF,int jF,int di,int dj,
                              const vector<vector<Node>>& dp) const
    {
        ofstream f("pred_last_layer.csv"); 
        f.setf(std::ios::fixed); 
        f<<setprecision(6);
        f<<"i,j,move,dV,dH,dt,fuel,dJ\n";
        int s=iF+jF-1;
        auto inside=[&](int i,int j){ return i>=0&&i<(int)Vg.size()&&j>=0&&j<(int)Hg.size(); };
        for(int i=0;i<(int)Vg.size();++i){
            int j=s-i; 
            if(!inside(i,j)) continue; 
            if(!std::isfinite(dp[i][j].cost)) continue;
            double V=Vg[i], H=Hg[j], m=dp[i][j].mass; 
            int in,jn;
            in=i+di; jn=j;
            if(inside(in,jn)){ 
                double dV=Vg[in]-V; 
                if(dV*di>0){ 
                    auto r=step_A(cfg,V,H,m,dV); 
                    if(r.ok) f<<i<<","<<j<<",A,"<<dV<<",0,"<<r.dt<<","<<r.fuel<<","<<J_inc(cfg,r.fuel,r.dt)<<"\n"; 
                } 
            }
            in=i; jn=j+dj;
            if(inside(in,jn)){ 
                double dH=Hg[jn]-H; 
                if(dH*dj>0){ 
                    auto r=step_B(cfg,V,H,m,dH); 
                    if(r.ok) f<<i<<","<<j<<",B,0,"<<dH<<","<<r.dt<<","<<r.fuel<<","<<J_inc(cfg,r.fuel,r.dt)<<"\n"; 
                } 
            }
            in=i+di; jn=j+dj;
            if(inside(in,jn)){ 
                double dV=Vg[in]-V, dH=Hg[jn]-H; 
                if(dV*di>0 && dH*dj>0){ 
                    auto r=step_C(cfg,V,H,m,dV,dH); 
                    if(r.ok) f<<i<<","<<j<<",C,"<<dV<<","<<dH<<","<<r.dt<<","<<r.fuel<<","<<J_inc(cfg,r.fuel,r.dt)<<"\n"; 
                } 
            }
        }
    }
    
    static void dumpEdgeCSV(const string& name,const vector<double>& Vg,const vector<double>& Hg,const EdgeMatrices& E){
        auto save=[&](const string& fname,const vector<vector<EdgeCost>>& M){
            ofstream f(fname); 
            f.setf(std::ios::fixed); 
            f<<setprecision(6);
            f<<"i,j,V,H,ok,cost,dt,fuel\n";
            for(int i=0;i<(int)Vg.size();++i)
                for(int j=0;j<(int)Hg.size();++j){
                    const auto&e=M[i][j];
                    f<<i<<","<<j<<","<<Vg[i]<<","<<Hg[j]<<","<<(e.ok?1:0)<<","<<e.cost<<","<<e.dt<<","<<e.fuel<<"\n";
                }
        }; 
        save(name+"_A.csv",E.A); 
        save(name+"_B.csv",E.B); 
        save(name+"_C.csv",E.C);
    }

    Result solve(){
        double Vmin=min(cfg.V0,cfg.Vf), Vmax=max(cfg.V0,cfg.Vf);
        double Hmin=min(cfg.H0,cfg.Hf), Hmax=max(cfg.H0,cfg.Hf);
        cfg.Vmin=min(cfg.Vmin,Vmin); cfg.Vmax=max(cfg.Vmax,Vmax);
        cfg.Hmin=min(cfg.Hmin,Hmin); cfg.Hmax=max(cfg.Hmax,Hmax);

        vector<double> Vg(cfg.NV), Hg(cfg.NH);
        for(int i=0;i<cfg.NV;i++) Vg[i]=cfg.Vmin + (cfg.Vmax-cfg.Vmin)*i/(cfg.NV-1);
        for(int j=0;j<cfg.NH;j++) Hg[j]=cfg.Hmin + (cfg.Hmax-cfg.Hmin)*j/(cfg.NH-1);

        auto idxN=[&](double x,const vector<double>& a){
            int n=a.size(),b=0; 
            double d=fabs(a[0]-x);
            for(int i=1;i<n;i++){ 
                double di=fabs(a[i]-x); 
                if(di<d){d=di;b=i;} 
            } 
            return b;
        };
        int i0=idxN(cfg.V0,Vg), j0=idxN(cfg.H0,Hg);
        int iF=idxN(cfg.Vf,Vg), jF=idxN(cfg.Hf,Hg);

        int di=(iF>=i0)?+1:-1, dj=(jF>=j0)?+1:-1;
        edges.init(cfg.NV,cfg.NH);

        vector<vector<Node>> dp(cfg.NV, vector<Node>(cfg.NH));
        dp[i0][j0].cost=0.0; 
        dp[i0][j0].time=0.0; 
        dp[i0][j0].mass=cfg.m0; 
        dp[i0][j0].V=Vg[i0]; 
        dp[i0][j0].H=Hg[j0];

        auto inside=[&](int i,int j){ return i>=0 && i<cfg.NV && j>=0 && j<cfg.NH; };

        int si0=i0+j0, siF=iF+jF, sdir=(siF>=si0)?+1:-1;
        for(int s=si0+sdir; sdir>0? s<=siF : s>=siF; s+=sdir){
            for(int i=0;i<cfg.NV;i++){
                int j=s-i; 
                if(!inside(i,j)) continue;
                Node best; 
                double bestCost=numeric_limits<double>::infinity();

                int ip=i-di, jp=j;
                if(inside(ip,jp) && std::isfinite(dp[ip][jp].cost)){
                    double Vprev=Vg[ip], Hprev=Hg[jp];
                    double dV=Vg[i]-Vprev; 
                    if(dV*di>0){
                        double mprev=dp[ip][jp].mass;
                        auto r=step_A(cfg,Vprev,Hprev,mprev,dV);
                        if(r.ok){
                            double inc=J_inc(cfg,r.fuel,r.dt); 
                            edges.A[i][j]={true,inc,r.dt,r.fuel};
                            double newCost=dp[ip][jp].cost+inc;
                            if(newCost<bestCost){
                                bestCost=newCost; best=dp[ip][jp];
                                best.cost=newCost; best.time+=r.dt; best.mass-=r.fuel;
                                best.pi=ip; best.pj=jp; best.move='A';
                                best.V=Vg[i]; best.H=Hg[j]; best.Theta=r.Theta;
                                best.alpha=r.alpha; best.thr=r.throttle; best.P=r.P; best.dt=r.dt;
                                best.CL=r.CL; best.CD=r.CD; best.q=r.q; best.rho=r.rho; best.a=r.a; best.g=r.g;
                            }
                        }
                    }
                }
                
                ip=i; jp=j-dj;
                if(inside(ip,jp) && std::isfinite(dp[ip][jp].cost)){
                    double Vprev=Vg[ip], Hprev=Hg[jp];
                    double dH=Hg[j]-Hprev; 
                    if(dH*dj>0){
                        double mprev=dp[ip][jp].mass;
                        auto r=step_B(cfg,Vprev,Hprev,mprev,dH);
                        if(r.ok){
                            double inc=J_inc(cfg,r.fuel,r.dt); 
                            edges.B[i][j]={true,inc,r.dt,r.fuel};
                            double newCost=dp[ip][jp].cost+inc;
                            if(newCost<bestCost){
                                bestCost=newCost; best=dp[ip][jp];
                                best.cost=newCost; best.time+=r.dt; best.mass-=r.fuel;
                                best.pi=ip; best.pj=jp; best.move='B';
                                best.V=Vg[i]; best.H=Hg[j]; best.Theta=r.Theta;
                                best.alpha=r.alpha; best.thr=r.throttle; best.P=r.P; best.dt=r.dt;
                                best.CL=r.CL; best.CD=r.CD; best.q=r.q; best.rho=r.rho; best.a=r.a; best.g=r.g;
                            }
                        }
                    }
                }
                
                ip=i-di; jp=j-dj;
                if(inside(ip,jp) && std::isfinite(dp[ip][jp].cost)){
                    double Vprev=Vg[ip], Hprev=Hg[jp];
                    double dV=Vg[i]-Vprev, dH=Hg[j]-Hprev;
                    if(dV*di>0 && dH*dj>0){
                        double mprev=dp[ip][jp].mass;
                        auto r=step_C(cfg,Vprev,Hprev,mprev,dV,dH);
                        if(r.ok){
                            double inc=J_inc(cfg,r.fuel,r.dt); 
                            edges.C[i][j]={true,inc,r.dt,r.fuel};
                            double newCost=dp[ip][jp].cost+inc;
                            if(newCost<bestCost){
                                bestCost=newCost; best=dp[ip][jp];
                                best.cost=newCost; best.time+=r.dt; best.mass-=r.fuel;
                                best.pi=ip; best.pj=jp; best.move='C';
                                best.V=Vg[i]; best.H=Hg[j]; best.Theta=r.Theta;
                                best.alpha=r.alpha; best.thr=r.throttle; best.P=r.P; best.dt=r.dt;
                                best.CL=r.CL; best.CD=r.CD; best.q=r.q; best.rho=r.rho; best.a=r.a; best.g=r.g;
                            }
                        }
                    }
                }
                if(std::isfinite(best.cost)) dp[i][j]=best;
            }
            if(debug_last_layer){
                if((sdir>0 && s==siF-1) || (sdir<0 && s==siF+1)){
                    dump_pred_last_layer(Vg,Hg,iF,jF,di,dj,dp);
                }
            }
        }

        Result R;
        if(!std::isfinite(dp[iF][jF].cost)){ R.ok=false; return R; }

        vector<pair<int,int>> path; 
        int ci=iF,cj=jF;
        while(!(ci==i0 && cj==j0)){
            path.push_back({ci,cj});
            int ni=dp[ci][cj].pi, nj=dp[ci][cj].pj;
            if(ni<0||nj<0){ R.ok=false; return R; }
            ci=ni; cj=nj;
        }
        path.push_back({i0,j0}); 
        reverse(path.begin(),path.end());

        R.ok=true; 
        int N=path.size();
        R.V.resize(N); R.H.resize(N); R.t.resize(N); R.m.resize(N);
        R.move.resize(N); R.alpha.resize(N); R.thr.resize(N); R.Theta.resize(N); R.dt.resize(N); R.P.resize(N);
        for(int k=0;k<N;k++){
            int i=path[k].first, j=path[k].second; 
            const auto& n=dp[i][j];
            R.V[k]=Vg[i]; R.H[k]=Hg[j]; R.t[k]=n.time; R.m[k]=n.mass; R.move[k]=n.move;
            R.alpha[k]=n.alpha; R.thr[k]=n.thr; R.Theta[k]=n.Theta; R.dt[k]=n.dt; R.P[k]=n.P;
        }
        R.total_time=R.t.back();
        R.fuel_used=cfg.m0 - R.m.back();

        dumpEdgeCSV("edges", Vg, Hg, edges);
        return R;
    }
};

int main(int argc, char** argv){
    ios::sync_with_stdio(false); 
    cin.tie(nullptr);

    Config cfg;
    cfg.time_as_fuel = cfg.engine.cfg.Tmax_SL * cfg.engine.cfg.sfc0 / 3600.0;

    int a=1; 
    string mode="vh"; 
    if(argc>a) mode=argv[a++];
    if(argc>a){
        string obj=argv[a++];
        if(obj=="fuel") cfg.objective=Objective::MinFuel;
        else if(obj=="time") cfg.objective=Objective::MinTime;
        else if(obj=="mix"){ 
            cfg.objective=Objective::Weighted; 
            if(argc>a) cfg.lambda_mix=atof(argv[a++]); 
        }
    }
    bool debug=false; 
    if(argc>a){ 
        string last=argv[a++]; 
        if(last=="debug") debug=true; 
    }

    SolverVH solver{cfg}; 
    solver.debug_last_layer = debug;
    auto R = solver.solve();
    if(!R.ok){ 
        cerr<<"Не удалось построить траекторию (проверьте сетки/ограничения).\n"; 
        return 1; 
    }

    cerr<<"Итог: Время="<<R.total_time<<" с, Расход="<<R.fuel_used<<" кг\n";
    cerr<<"Финал: V="<<R.V.back()<<" м/с, H="<<R.H.back()<<" м, m="<<R.m.back()<<" кг\n";
    saveCSV("trajectory_vh.csv", R);
    cerr<<"CSV: trajectory_vh.csv; матрицы весов: edges_A/B/C.csv";
    if (debug) cerr<<"; отладка: pred_last_layer.csv";
    cerr<<"\n";
    return 0;
}
