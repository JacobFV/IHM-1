#pragma once
#include <OpenSim/OpenSim.h>
#include "native_bed_compression.h"
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <cmath>

namespace ihm_surface {
struct Point {SimTK::Vec3 station;double area;int index=-1;bool selected=false;};
struct Observation {int index;std::string body;SimTK::Vec3 location,force;double area,indentation,bed_indentation=0,total_approach=0;};
struct Group {std::string body;std::vector<Point> points;bool has_sensors=false;SimTK::Vec3 low{SimTK::Infinity},high{-SimTK::Infinity};};
struct Wrench {SimTK::Vec3 force{0},moment{0};double maximum_penetration=0;int contacting_points=0;};
struct Sample {std::vector<Observation> observations;std::map<std::string,Wrench> bodies;SimTK::Vec3 force{0},bed_moment{0};double energy=0,skin_energy=0,bed_energy=0,maximum_bed_indentation=0,maximum_total_approach=0,power=0,dissipative_power=0,maximum_penetration=0;int contacting_points=0;};
class Foundation final:public OpenSim::Force {
    OpenSim_DECLARE_CONCRETE_OBJECT(Foundation,OpenSim::Force);
public:
    double plane=0,h=0,mu=0,lambda=0,minimum_ratio=0,dissipation=0,friction=0,viscous=0,transition=0;
    std::vector<Group> groups;ihm_bed::Curve bed;bool bed_enabled=false;
    void readBed(const std::string& path){bed.read(path);bed_enabled=true;}
    void read(const std::string& path,const OpenSim::Model& model) {
        std::ifstream input(path);std::string version;int count;
        if (!(input>>version>>plane>>h>>mu>>lambda>>minimum_ratio>>dissipation>>friction>>viscous>>transition>>count)
            ||version!="IHM_SURFACE_FOUNDATION_V1"||count<1||count>50000)
            throw std::runtime_error("invalid surface foundation input header");
        for(double value:{plane,h,mu,lambda,minimum_ratio,dissipation,friction,viscous,transition})
            if(!std::isfinite(value))throw std::runtime_error("nonfinite surface foundation parameter");
        if(h<=0||mu<=0||lambda<0||minimum_ratio<=0||minimum_ratio>=1||dissipation<0||friction<0||viscous<0||transition<=0)
            throw std::runtime_error("invalid surface foundation material/domain");
        std::map<std::string,int> indices;
        for(int i=0;i<count;++i){
            std::string body;Point point;
            if(!(input>>body>>point.station[0]>>point.station[1]>>point.station[2]>>point.area)
                ||!point.station.isFinite()||!std::isfinite(point.area)||point.area<=0)
                throw std::runtime_error("invalid surface foundation quadrature");
            model.getBodySet().get(body);
            if(!indices.count(body)){indices[body]=(int)groups.size();groups.push_back(Group{});groups.back().body=body;}
            auto& group=groups[indices.at(body)];point.index=i;group.points.push_back(point);
            for(int k=0;k<3;k++){group.low[k]=std::min(group.low[k],point.station[k]);group.high[k]=std::max(group.high[k],point.station[k]);}
        }
        std::string extra;if(input>>extra)throw std::runtime_error("trailing surface foundation input");
    }
    void select(const std::string& path){
        std::ifstream input(path);int count;
        if(!(input>>count)||count<0||count>128)throw std::runtime_error("invalid selected skin sensor count");
        std::set<int> selected;for(int i=0;i<count;i++){int index;if(!(input>>index)||index<0||!selected.insert(index).second)throw std::runtime_error("invalid selected skin sensor index");}
        std::string extra;if(input>>extra)throw std::runtime_error("trailing skin sensor selection");
        for(auto& group:groups)for(auto& point:group.points)if(selected.erase(point.index)){point.selected=true;group.has_sensors=true;}
        if(!selected.empty())throw std::runtime_error("selected skin sensor outside quadrature");
    }
    Sample sample(const SimTK::State& state,bool with_velocity=true) const {
        Sample result;
        for(const auto& group:groups){
            const auto& body=getModel().getBodySet().get(group.body);const auto transform=body.getTransformInGround(state);
            const auto rotation=transform.R().asMat33();const auto origin=transform.p();
            double minimum_x=origin[0];
            for(int k=0;k<3;k++)minimum_x+=rotation(0,k)*(rotation(0,k)>=0?group.low[k]:group.high[k]);
            auto& wrench=result.bodies[group.body];
            if(minimum_x>=plane&&!group.has_sensors)continue;
            SimTK::SpatialVec velocity(SimTK::Vec3(0),SimTK::Vec3(0));
            if(with_velocity)velocity=body.getMobilizedBody().getBodyVelocity(state);
            for(const auto& point:group.points){
                if(minimum_x>=plane&&!point.selected)continue;
                const auto offset=transform.R()*point.station;const auto location=origin+offset;
                const double penetration=plane-location[0];if(penetration<=0){if(point.selected)result.observations.push_back({point.index,group.body,location,SimTK::Vec3(0),point.area,0.});continue;}
                double skin_indentation=penetration,bed_indentation=0,skin_energy=0,bed_energy=0,pressure=0;
                if(bed_enabled){const auto response=bed.solve(penetration,h,mu,lambda,minimum_ratio);
                    skin_indentation=response.skin_indentation;bed_indentation=response.bed_indentation;
                    skin_energy=response.skin_energy;bed_energy=response.bed_energy;pressure=response.pressure;
                }else{const double stretch=1-penetration/h;
                    if(stretch<minimum_ratio-1e-12)throw std::runtime_error("surface foundation compression exceeds declared domain");
                    const double log=std::log(stretch);pressure=-mu*(stretch-1/stretch)-lambda*log/stretch;
                    skin_energy=h*(.5*mu*(stretch*stretch-1)-mu*log+.5*lambda*log*log);}
                const double elastic=point.area*pressure;
                const auto speed=velocity[1]+cross(velocity[0],offset);
                const double normal=std::max(0.,elastic*(1-(bed_enabled?0.:dissipation)*speed[0]));
                const double tangent=std::hypot(speed[1],speed[2]);
                const double coefficient=friction*std::tanh(tangent/transition)+viscous*tangent;
                SimTK::Vec3 force(normal,0,0);
                if(tangent>0){force[1]=-normal*coefficient*speed[1]/tangent;force[2]=-normal*coefficient*speed[2]/tangent;}
                if(point.selected)result.observations.push_back({point.index,group.body,location,force,point.area,skin_indentation,bed_indentation,penetration});
                wrench.force+=force;wrench.moment+=cross(offset,force);wrench.contacting_points++;
                wrench.maximum_penetration=std::max(wrench.maximum_penetration,skin_indentation);
                result.force+=force;result.bed_moment-=cross(location,force);result.contacting_points++;
                result.maximum_penetration=std::max(result.maximum_penetration,skin_indentation);
                result.skin_energy+=point.area*skin_energy;result.bed_energy+=point.area*bed_energy;result.energy+=point.area*(skin_energy+bed_energy);
                result.maximum_bed_indentation=std::max(result.maximum_bed_indentation,bed_indentation);result.maximum_total_approach=std::max(result.maximum_total_approach,penetration);
                result.power+=dot(force,speed);result.dissipative_power+=dot(force-SimTK::Vec3(elastic,0,0),speed);
            }
        }
        if(result.dissipative_power>1e-8)throw std::runtime_error("surface foundation created dissipative energy");
        return result;
    }
    void computeForce(const SimTK::State& state,SimTK::Vector_<SimTK::SpatialVec>& forces,SimTK::Vector&) const override {
        const auto result=sample(state);
        for(const auto& value:result.bodies){const auto& body=getModel().getBodySet().get(value.first);
            forces[body.getMobilizedBodyIndex()][0]+=value.second.moment;
            forces[body.getMobilizedBodyIndex()][1]+=value.second.force;}
        forces[0][0]+=result.bed_moment;forces[0][1]-=result.force;
    }
    double computePotentialEnergy(const SimTK::State& state) const override {return sample(state,false).energy;}
};
} // namespace ihm_surface
