#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
using namespace Rcpp;
struct Ranking {
 std::vector<double>s,t;std::vector<int>id,y,o,ot;bool good;
 Ranking(List z){NumericVector a=z["score"];IntegerVector ii=z["id"],yy=z["label"];s=as<std::vector<double>>(a);id=as<std::vector<int>>(ii);y=as<std::vector<int>>(yy);for(int &i:id)--i;good=s.size()>0;for(double v:s)if(!R_finite(v))good=false;o.resize(s.size());std::iota(o.begin(),o.end(),0);if(good)std::stable_sort(o.begin(),o.end(),[&](int i,int j){return s[i]<s[j];});if(z.containsElementNamed("second")){t=as<std::vector<double>>(z["second"]);for(double v:t)if(!R_finite(v))good=false;ot=o;if(good)std::stable_sort(ot.begin(),ot.end(),[&](int i,int j){return t[i]<t[j];});}}
 std::vector<double> calc(const std::vector<double>&w,int metric){std::vector<double>v(s.size(),NA_REAL);if(!good)return v;
 if(metric==4){double pos=0,neg=0;for(size_t j=0;j<s.size();j++)(y[j]?pos:neg)+=w[id[j]];if(pos<=0||neg<=0)return v;double bp=0,bn=0;
 for(size_t j=0;j<o.size();){size_t k=j;double pp=0,nn=0;while(k<o.size()&&s[o[k]]==s[o[j]]){int a=o[k++];(y[a]?pp:nn)+=w[id[a]];}for(size_t l=j;l<k;l++){int a=o[l];v[a]=y[a]?(bn+.5*nn)/neg:(pos-bp-.5*pp)/pos;}bp+=pp;bn+=nn;j=k;}
 }else{double N=0;for(int i:id)N+=w[i];if(N<3)return v;std::vector<double>rx(s.size()),ry(s.size());auto ranks=[&](std::vector<double>&r,const std::vector<double>&a,const std::vector<int>&ord){double below=0;for(size_t j=0;j<ord.size();){size_t k=j;double tie=0;while(k<ord.size()&&a[ord[k]]==a[ord[j]])tie+=w[id[ord[k++]]];for(size_t l=j;l<k;l++)r[ord[l]]=below+(tie+1)/2;below+=tie;j=k;}};ranks(rx,s,o);ranks(ry,t,ot);double mean=(N+1)/2,xx=0,yy=0;for(size_t j=0;j<s.size();j++){xx+=w[id[j]]*(rx[j]-mean)*(rx[j]-mean);yy+=w[id[j]]*(ry[j]-mean)*(ry[j]-mean);}if(xx<=0||yy<=0)return v;for(size_t j=0;j<s.size();j++)v[j]=N*(rx[j]-mean)*(ry[j]-mean)/sqrt(xx*yy);
 }return v;}
};
// [[Rcpp::export]]
NumericVector reassay_draws(List unit_segments,IntegerVector repetitions,NumericMatrix weights,int metric){
 int G=weights.ncol(),D=weights.nrow(),R=max(repetitions);std::vector<std::vector<Ranking>>units;for(int u=0;u<unit_segments.size();u++){List ss=unit_segments[u];std::vector<Ranking>z;for(int j=0;j<ss.size();j++)z.emplace_back(as<List>(ss[j]));units.push_back(z);}NumericVector out(D,NA_REAL);
 for(int d=0;d<D;d++){std::vector<double>w(G);for(int i=0;i<G;i++)w[i]=weights(d,i);std::vector<std::vector<double>>rs(R,std::vector<double>(G)),rc(R,std::vector<double>(G));
 for(size_t u=0;u<units.size();u++){std::vector<double>su(G),nu(G);std::vector<bool>bad(G,false);for(Ranking &z:units[u]){auto v=z.calc(w,metric);for(size_t j=0;j<v.size();j++){int i=z.id[j];if(R_finite(v[j])){su[i]+=v[j];nu[i]++;}else bad[i]=true;}}int r=repetitions[u]-1;for(int i=0;i<G;i++)if(nu[i]>0&&!bad[i]){rs[r][i]+=su[i]/nu[i];rc[r][i]++;}}
 double num=0,den=0;for(int i=0;i<G;i++){double s=0,n=0;for(int r=0;r<R;r++)if(rc[r][i]>0){s+=rs[r][i]/rc[r][i];n++;}if(n>0){num+=w[i]*s/n;den+=w[i];}}if(den>0)out[d]=num/den;
 }return out;
}
