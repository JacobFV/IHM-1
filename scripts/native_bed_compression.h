#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
namespace ihm_bed {
struct Response {double pressure=0,skin_indentation=0,bed_indentation=0,skin_energy=0,bed_energy=0;};
class Curve {
public:
    std::string material;double thickness=0;std::vector<double> strain,pressure,integral;
    void read(const std::string& path){
        std::ifstream input(path);std::string version;int count;
        if(!(input>>version>>material>>thickness>>count)||version!="IHM_BED_COMPRESSION_V1"||count<2||count>100||!std::isfinite(thickness)||thickness<=0)
            throw std::runtime_error("invalid measured bed curve header");
        for(int i=0;i<count;i++){double e,p;if(!(input>>e>>p)||!std::isfinite(e)||!std::isfinite(p)||e<0||e>=1||p<0||(i&&(e<=strain.back()||p<pressure.back())))throw std::runtime_error("invalid monotone finite bed curve");strain.push_back(e);pressure.push_back(p);}
        if(strain[0]!=0||pressure[0]!=0)throw std::runtime_error("bed curve needs explicit zero-load origin");
        std::string extra;if(input>>extra)throw std::runtime_error("trailing bed curve input");
        integral.push_back(0);for(int i=1;i<count;i++)integral.push_back(integral.back()+.5*(pressure[i-1]+pressure[i])*(strain[i]-strain[i-1]));
    }
    int interval(double e) const{return std::max(0,std::min(int(strain.size())-2,int(std::upper_bound(strain.begin(),strain.end(),e)-strain.begin())-1));}
    double stress(double e) const {const int i=interval(e);return pressure[i]+(e-strain[i])*(pressure[i+1]-pressure[i])/(strain[i+1]-strain[i]);}
    double energy(double e) const {const int i=interval(e);const double d=e-strain[i],slope=(pressure[i+1]-pressure[i])/(strain[i+1]-strain[i]);return thickness*(integral[i]+pressure[i]*d+.5*slope*d*d);}
    Response solve(double approach,double h,double mu,double lambda,double minimum_ratio) const {
        Response result;if(approach<=0)return result;
        const double skin_max=h*(1-minimum_ratio),bed_max=thickness*strain.back();
        double low=std::max(0.,approach-bed_max),high=std::min(skin_max,approach);
        if(low>high+1e-12)throw std::runtime_error("combined skin/bed compression outside retained domains");
        const auto skin_pressure=[&](double d){const double stretch=1-d/h;return -mu*(stretch-1/stretch)-lambda*std::log(stretch)/stretch;};
        const auto residual=[&](double d){return skin_pressure(d)-stress((approach-d)/thickness);};
        if(residual(low)>1e-7||residual(high)<-1e-7)throw std::runtime_error("equal-pressure skin/bed solution exceeds retained domains");
        // Safeguarded Newton solves the unchanged monotone pressure equation.
        // Retain the original bracket/bisection for nonconvergence and endpoints.
        const double original_low=low,original_high=high;
        double indentation=.5*(low+high);bool converged=false;
        for(int iteration=0;iteration<50;iteration++){
            const double value=residual(indentation);
            if(std::abs(value)<=1e-10){converged=true;break;}
            if(value>0)high=indentation;else low=indentation;
            const double stretch=1-indentation/h;
            const int i=interval((approach-indentation)/thickness);
            const double derivative=(mu*(1+1/(stretch*stretch))+
                lambda*(1-std::log(stretch))/(stretch*stretch))/h+
                (pressure[i+1]-pressure[i])/(thickness*(strain[i+1]-strain[i]));
            const double proposed=indentation-value/derivative;
            const double next=std::isfinite(proposed)&&proposed>low&&proposed<high?
                proposed:.5*(low+high);
            if(next==indentation)break;
            indentation=next;
        }
        if(!converged){
            low=original_low;high=original_high;
            for(int i=0;i<50;i++){const double middle=.5*(low+high);if(residual(middle)>0)high=middle;else low=middle;}
            indentation=.5*(low+high);
        }
        result.skin_indentation=indentation;result.bed_indentation=approach-result.skin_indentation;
        result.pressure=skin_pressure(result.skin_indentation);const double stretch=1-result.skin_indentation/h,log=std::log(stretch);
        result.skin_energy=h*(.5*mu*(stretch*stretch-1)-mu*log+.5*lambda*log*log);
        result.bed_energy=energy(result.bed_indentation/thickness);return result;
    }
};
} // namespace ihm_bed
