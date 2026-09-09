 // Finite differences include activation/fiber states. Muscle equilibration is
 // performed ONLY at the base point, never after perturbing a dynamic state.
 auto names=model.getStateVariableNames();std::vector<int> selected;
 if(fs::exists(source/"reference_state.txt")){
  std::ifstream input(source/"reference_state.txt");std::string schema;int count;
  input>>schema>>count;if(!input||schema!="IHM_REFERENCE_STATE_V1"||count!=names.size())throw std::runtime_error("Reference state requires every native state variable");
  std::map<std::string,double> reference;
  for(int i=0;i<count;i++){std::string name;double value;input>>name>>value;if(!input||!std::isfinite(value)||reference.count(name))throw std::runtime_error("Invalid or duplicate reference state variable");reference[name]=value;}
  std::string extra;if(input>>extra)throw std::runtime_error("Trailing reference state data");
  auto requested=model.getStateVariableValues(state);
  for(int i=0;i<names.size();i++){if(!reference.count(names[i]))throw std::runtime_error("Missing reference state variable "+names[i]);requested[i]=reference.at(names[i]);}
  model.setStateVariableValues(state,requested);
  model.getMultibodySystem().projectQ(state,1e-12);model.getMultibodySystem().projectU(state,1e-12);
  if((model.getStateVariableValues(state)-requested).norm()>1e-8)throw std::runtime_error("Reference state violates native position/velocity constraints");
  // Do not equilibrate: the provided copied-state fiber solution is the point
  // to linearize, and replacing it would recreate the initialization mismatch.
 }
 for(int i=0;i<names.size();i++){
  bool skip=false;
  for(const auto& c:model.getComponentList<Coordinate>()){
   const auto prefix=c.getAbsolutePathString()+"/";
   if(names[i].find(prefix)==0&&(c.isDependent(state)||c.getLocked(state)||c.isPrescribed(state)))skip=true;
  }
  if(!skip)selected.push_back(i);
 }
 std::vector<std::string> muscles;for(const auto& m:model.getComponentList<Muscle>())muscles.push_back(m.getName());
 auto x0=model.getStateVariableValues(state);
 auto values=[&](const SimTK::State& s){auto all=model.getStateVariableValues(s);SimTK::Vector v(selected.size());for(int j=0;j<(int)selected.size();j++)v[j]=all[selected[j]];return v;};
 auto derivative=[&](SimTK::State& s){model.realizeAcceleration(s);SimTK::Vector f(selected.size());for(int j=0;j<(int)selected.size();j++)f[j]=model.getStateVariableDerivativeValue(s,names[selected[j]]);return f;};
 auto jsonv=[&](std::ostream& o,const SimTK::Vector& v){o<<'[';for(int i=0;i<v.size();i++){if(i)o<<',';num(o,v[i]);}o<<']';};
 auto f0=derivative(state);
 std::ofstream result(out/"linearization.json");result<<"{\"state_names\":[";
 for(int j=0;j<(int)selected.size();j++){if(j)result<<',';str(result,names[selected[j]]);}result<<"],\"muscle_names\":[";
 for(int j=0;j<(int)muscles.size();j++){if(j)result<<',';str(result,muscles[j]);}result<<"],\"x0\":";jsonv(result,values(state));result<<",\"f0\":";jsonv(result,f0);
 result<<",\"u0\":[";for(int j=0;j<(int)muscles.size();j++){if(j)result<<',';num(result,excitation->values.at(muscles[j]));}
 result<<"],\"minimum_activation\":[";
 for(int j=0;j<(int)muscles.size();j++){
  if(j)result<<',';const auto& m=model.getMuscles().get(muscles[j]);
  if(const auto* millard=dynamic_cast<const Millard2012EquilibriumMuscle*>(&m))num(result,millard->getMinimumActivation());
  else if(const auto* thelen=dynamic_cast<const Thelen2003Muscle*>(&m))num(result,thelen->getMinimumActivation());
  else throw std::runtime_error("Minimum activation unavailable for native muscle type "+m.getConcreteClassName());
 }
 result<<"],\"state_columns\":[";
 for(int j=0;j<(int)selected.size();j++){
  double h=1e-6;SimTK::State plus=state,minus=state;auto xp=x0,xm=x0;xp[selected[j]]+=h;xm[selected[j]]-=h;
  model.setStateVariableValues(plus,xp);model.setStateVariableValues(minus,xm);
  model.getMultibodySystem().projectQ(plus,1e-12);model.getMultibodySystem().projectU(plus,1e-12);
  model.getMultibodySystem().projectQ(minus,1e-12);model.getMultibodySystem().projectU(minus,1e-12);
  auto dx=(values(plus)-values(minus))/(2*h);auto df=(derivative(plus)-derivative(minus))/(2*h);
  if(j)result<<',';result<<"{\"dx\":";jsonv(result,dx);result<<",\"df\":";jsonv(result,df);result<<'}';
 }
 result<<"],\"input_columns\":[";
 for(int j=0;j<(int)muscles.size();j++){
  double h=1e-6,base=excitation->values.at(muscles[j]);SimTK::State plus=state,minus=state;
  excitation->values[muscles[j]]=base+h;plus.invalidateAllCacheAtOrAbove(SimTK::Stage::Dynamics);auto fp=derivative(plus);
  excitation->values[muscles[j]]=base-h;minus.invalidateAllCacheAtOrAbove(SimTK::Stage::Dynamics);auto fm=derivative(minus);excitation->values[muscles[j]]=base;
  if(j)result<<',';jsonv(result,(fp-fm)/(2*h));
 }
 result<<"]}";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
