// Exact CPU implementation of execution_graph.py's finite event rules.
// Build without fast-math or contraction: Python's ordered binary64 operations
// are part of the reference semantics, including event ties and wait totals.
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
struct Unsupported : std::runtime_error {using std::runtime_error::runtime_error;};
struct PythonError {};
struct Ref {
 PyObject* p;
 explicit Ref(PyObject* value):p(value){if(!p)throw PythonError();}
 ~Ref(){Py_DECREF(p);}
 Ref(const Ref&)=delete;
 operator PyObject*() const{return p;}
};
std::string text(PyObject* value){
 if(!PyUnicode_Check(value))throw Unsupported("Native scheduler requires string keys");
 Py_ssize_t length;const char* p=PyUnicode_AsUTF8AndSize(value,&length);if(!p)throw PythonError();std::string result(p,length);if(result.find(char(0))!=std::string::npos)throw Unsupported("Embedded NUL key");return result;
}
PyObject* mapping(PyObject* value){if(!PyDict_Check(value))throw Unsupported("Non-dictionary graph field");return value;}
PyObject* field(PyObject* value,const char* key){
 PyObject* p=PyDict_GetItemString(mapping(value),key);if(!p)throw Unsupported("Missing graph field");return p;
}
double number(PyObject* p){double value=PyFloat_AsDouble(p);if(PyErr_Occurred())throw PythonError();return value;}
int64_t positive_integer(PyObject* p,const char* error="Token capacities must be positive integers"){
 if(PyBool_Check(p)||!PyLong_Check(p))throw std::invalid_argument(error);
 int64_t value=PyLong_AsLongLong(p);if(PyErr_Occurred())throw Unsupported("Large Python integer");
 if(value<1)throw std::invalid_argument(error);return value;
}
struct Demand{int resource;double amount;};
struct Action{int token;int64_t amount;};
struct Node{
 std::string name;double seconds=0.;std::vector<int> deps,children;
 std::vector<Demand> demands;std::array<std::vector<Action>,3> actions;
 std::vector<int> waits_after,waits_until;
 int dequeue=-1,enqueue=-1,consumer=-1,attempt=-1,ready=-1;
};
struct Wait{int after,until;std::vector<Demand> demands;};
struct Token{std::string name;int64_t capacity;int shared=-1;};
struct Graph{
 std::vector<Node> nodes;std::vector<Wait> waits;std::vector<Token> tokens;
 std::unordered_map<std::string,int> node_ids,token_ids,queue_ids;
 int complete=-1;
};
struct Item{int64_t order;std::shared_ptr<Graph> graph;};
struct NodeRef{int slot,node;};
struct WaitRef{int slot,wait;};
struct Trace{int item,node;double time;};
struct Conditional{int item,node;bool blocked;double waited,duration;};
struct Active{NodeRef ref;double remaining,rate=1.;};
struct Part{
 int item;std::vector<int> degree;std::vector<double> ends;std::vector<unsigned char> ended;
 std::vector<int64_t> free;std::vector<std::deque<int>> queues;size_t finished=0;
};
struct Solver{
 bool streamed,trace;std::string prefix;
 std::vector<std::string> resource_names;std::vector<double> capacities,usage;
 std::vector<unsigned char> usage_seen,capacity_given;std::vector<int> usage_order;
 std::unordered_map<std::string,int> resource_ids,shared_ids;
 std::unordered_set<std::string> shared_names;
 std::vector<int64_t> shared_capacity,shared_free;
 std::unordered_map<PyObject*,std::shared_ptr<Graph>> templates;
 std::vector<Item> items;std::vector<std::vector<int>> chains;std::vector<size_t> next;
 std::vector<std::unique_ptr<Part>> parts;
 std::deque<NodeRef> ready;std::vector<Active> active;
 std::map<std::pair<int64_t,int>,WaitRef> waits;
 std::vector<Trace> starts,ends;std::vector<Conditional> conditional;
 double now=0.,handoff_seconds=0.;int64_t steps=0,added=0,finished=0,live=0,peak=0,possible=0,blocked=0;
 int resource(const std::string& name){
  auto found=resource_ids.find(name);if(found!=resource_ids.end())return found->second;
  int id=resource_names.size();resource_ids[name]=id;resource_names.push_back(name);capacities.push_back(0.);capacity_given.push_back(0);return id;
 }
 void read_capacities(PyObject* graph){
  Ref values(PyObject_GetAttrString(graph,"capacities"));mapping(values);
  PyObject *key,*value;Py_ssize_t pos=0;
  while(PyDict_Next(values,&pos,&key,&value)){
   int id=resource(text(key));double cap=number(value);
   if(!std::isfinite(cap))throw Unsupported("Nonfinite capacity uses Python reference");
   if(capacity_given[id] && capacities[id]!=cap)throw std::invalid_argument("Conflicting shared resource capacity");
   capacities[id]=cap;capacity_given[id]=1;
  }
 }
 std::vector<Demand> demands(PyObject* row,bool waiting=false){
  std::vector<Demand> result;mapping(row);PyObject *key,*value;Py_ssize_t pos=0;
  while(PyDict_Next(row,&pos,&key,&value)){
   int id=resource(text(key));double amount=number(value);
   if(!std::isfinite(amount)||amount<0.)throw std::invalid_argument(waiting?"Invalid resource-wait demand or missing capacity":"Invalid resource demand");
   if(amount && capacities[id]<=0.)throw std::invalid_argument(waiting?"Invalid resource-wait demand or missing capacity":"Missing or zero resource "+resource_names[id]);
   result.push_back({id,amount});
  }
  return result;
 }
 int node_id(Graph& g,PyObject* name,const char* error){
  auto it=g.node_ids.find(text(name));if(it==g.node_ids.end())throw std::invalid_argument(error);return it->second;
 }
 int queue_id(Graph& g,PyObject* name){
  std::string key=text(name);auto it=g.queue_ids.find(key);if(it!=g.queue_ids.end())return it->second;
  int id=g.queue_ids.size();g.queue_ids[key]=id;return id;
 }
 std::shared_ptr<Graph> compile(PyObject* object){
  auto cached=templates.find(object);if(cached!=templates.end())return cached->second;
  auto g=std::make_shared<Graph>();
  Ref raw_nodes(PyObject_GetAttrString(object,"nodes"));mapping(raw_nodes);
  PyObject *key,*value;Py_ssize_t pos=0;
  while(PyDict_Next(raw_nodes,&pos,&key,&value)){
   std::string name=text(key);if(streamed && name=="complete")throw std::invalid_argument("Graph completion name collision");
   Ref row(PySequence_Fast(value,"Expected node tuple"));if(PySequence_Fast_GET_SIZE(row.p)!=2)throw Unsupported("Node tuple size");
   double seconds=number(PySequence_Fast_GET_ITEM(row.p,0));
   if(!std::isfinite(seconds)||seconds<0.)throw Unsupported("Modified node duration");
   int id=g->nodes.size();g->node_ids[name]=id;Node node;node.name=name;node.seconds=seconds;g->nodes.push_back(std::move(node));
  }
  pos=0;int index=0;
  while(PyDict_Next(raw_nodes,&pos,&key,&value)){
   Ref row(PySequence_Fast(value,"Expected node tuple"));Ref deps(PySequence_Fast(PySequence_Fast_GET_ITEM(row.p,1),"Expected dependencies"));
   std::unordered_set<int> seen;
   for(Py_ssize_t j=0;j<PySequence_Fast_GET_SIZE(deps.p);++j){
    int dep=node_id(*g,PySequence_Fast_GET_ITEM(deps.p,j),"Unknown predecessor");
    if(seen.insert(dep).second){g->nodes[index].deps.push_back(dep);g->nodes[dep].children.push_back(index);}
   }
   ++index;
  }
  Ref raw_demands(PyObject_GetAttrString(object,"demands"));mapping(raw_demands);pos=0;
  while(PyDict_Next(raw_demands,&pos,&key,&value)){
   auto it=g->node_ids.find(text(key));if(it==g->node_ids.end())continue;
   g->nodes[it->second].demands=demands(value);
  }
  Ref raw_tokens(PyObject_GetAttrString(object,"token_capacities"));mapping(raw_tokens);pos=0;
  while(PyDict_Next(raw_tokens,&pos,&key,&value)){
   std::string name=text(key);int64_t cap=positive_integer(value);int shared=-1;
   if(streamed && shared_names.count(name)){
    auto it=shared_ids.find(name);
    if(it==shared_ids.end()){shared=shared_ids.size();shared_ids[name]=shared;shared_capacity.push_back(cap);shared_free.push_back(cap);}
    else{shared=it->second;if(shared_capacity[shared]!=cap)throw std::invalid_argument("Token name collision or conflicting shared capacity");}
   }
   g->token_ids[name]=g->tokens.size();g->tokens.push_back({name,cap,shared});
  }
  Ref raw_actions(PyObject_GetAttrString(object,"token_actions"));mapping(raw_actions);pos=0;
  while(PyDict_Next(raw_actions,&pos,&key,&value)){
   int node=node_id(*g,key,"Unknown token node");mapping(value);PyObject *kind,*requests;Py_ssize_t ap=0;
   while(PyDict_Next(value,&ap,&kind,&requests)){
    std::string action=text(kind);int phase=action=="acquire"?0:action=="release_start"?1:action=="release_finish"?2:-1;
    if(phase<0)throw std::invalid_argument("Unknown token action");
    mapping(requests);PyObject *token,*amount;Py_ssize_t tp=0;
    while(PyDict_Next(requests,&tp,&token,&amount)){
     auto it=g->token_ids.find(text(token));if(it==g->token_ids.end())throw Unsupported("Token may be declared by another graph");
     int64_t count=positive_integer(amount,"Invalid token request");if(count>g->tokens[it->second].capacity)throw std::invalid_argument("Invalid token request");
     g->nodes[node].actions[phase].push_back({it->second,count});
    }
   }
  }
  Ref raw_enqueues(PyObject_GetAttrString(object,"fifo_enqueues"));mapping(raw_enqueues);pos=0;
  while(PyDict_Next(raw_enqueues,&pos,&key,&value)){
   auto it=g->node_ids.find(text(key));if(it==g->node_ids.end())throw Unsupported("Unused FIFO enqueue key");
   Ref row(PySequence_Fast(value,"Expected FIFO tuple"));if(PySequence_Fast_GET_SIZE(row.p)!=2)throw Unsupported("FIFO tuple size");
   auto consumer=g->node_ids.find(text(PySequence_Fast_GET_ITEM(row.p,1)));
   if(consumer==g->node_ids.end())throw Unsupported("Unused FIFO consumer");
   g->nodes[it->second].enqueue=queue_id(*g,PySequence_Fast_GET_ITEM(row.p,0));g->nodes[it->second].consumer=consumer->second;
  }
  Ref raw_dequeues(PyObject_GetAttrString(object,"fifo_dequeues"));mapping(raw_dequeues);pos=0;
  while(PyDict_Next(raw_dequeues,&pos,&key,&value)){
   auto it=g->node_ids.find(text(key));if(it==g->node_ids.end())throw Unsupported("Unused FIFO dequeue key");
   if(value!=Py_None)g->nodes[it->second].dequeue=queue_id(*g,value);
  }
  Ref raw_waits(PyObject_GetAttrString(object,"resource_waits"));mapping(raw_waits);pos=0;
  while(PyDict_Next(raw_waits,&pos,&key,&value)){
   if(streamed)text(key); // Prefix composition requires string wait names.
   int after=node_id(*g,field(value,"after"),"Unknown resource-wait endpoint");
   int until=node_id(*g,field(value,"until"),"Unknown resource-wait endpoint");
   int id=g->waits.size();g->waits.push_back({after,until,demands(field(value,"resources"),true)});
   g->nodes[after].waits_after.push_back(id);g->nodes[until].waits_until.push_back(id);
  }
  Ref raw_conditions(PyObject_GetAttrString(object,"conditional_delays"));mapping(raw_conditions);pos=0;
  while(PyDict_Next(raw_conditions,&pos,&key,&value)){
   int id=node_id(*g,key,"Unknown conditional-delay node");Node& node=g->nodes[id];
   node.attempt=node_id(*g,field(value,"attempt"),"Unknown conditional-delay endpoint");
   node.ready=node_id(*g,field(value,"ready"),"Unknown conditional-delay endpoint");
   for(int endpoint:{node.attempt,node.ready})if(std::find(node.deps.begin(),node.deps.end(),endpoint)==node.deps.end())
    throw std::invalid_argument("Conditional-delay endpoints must be direct dependencies");
  }
  if(streamed){
   Node complete;complete.name="complete";g->complete=g->nodes.size();
   for(size_t i=0;i<g->nodes.size();++i)if(g->nodes[i].children.empty()){
    complete.deps.push_back(i);g->nodes[i].children.push_back(g->complete);
   }
   g->nodes.push_back(std::move(complete));
  }
  templates[object]=g;return g;
 }
 Solver(PyObject* graph,PyObject* raw_chains,PyObject* shared,bool do_trace,const std::string& p):streamed(raw_chains!=Py_None),trace(do_trace),prefix(p){
  Ref shared_seq(PySequence_List(shared));
  for(Py_ssize_t i=0;i<PyList_GET_SIZE(shared_seq.p);++i)shared_names.insert(text(PyList_GET_ITEM(shared_seq.p,i)));
  read_capacities(graph);
  if(!streamed){items.push_back({0,compile(graph)});chains={{0}};}
  else{
   Ref list(PySequence_Fast(raw_chains,"Expected chains"));std::unordered_set<int64_t> orders;
   // Capacities are fixed across all future graphs, just as solve_chains checks
   // before admitting any work. Templates retain their input insertion order.
   for(Py_ssize_t i=0;i<PySequence_Fast_GET_SIZE(list.p);++i){
    Ref chain(PySequence_Fast(PySequence_Fast_GET_ITEM(list.p,i),"Expected chain"));
    for(Py_ssize_t j=0;j<PySequence_Fast_GET_SIZE(chain.p);++j){
     Ref pair(PySequence_Fast(PySequence_Fast_GET_ITEM(chain.p,j),"Expected ordered graph"));
     if(PySequence_Fast_GET_SIZE(pair.p)!=2)throw Unsupported("Chain tuple size");
     read_capacities(PySequence_Fast_GET_ITEM(pair.p,1));
    }
   }
   for(Py_ssize_t i=0;i<PySequence_Fast_GET_SIZE(list.p);++i){
    Ref chain(PySequence_Fast(PySequence_Fast_GET_ITEM(list.p,i),"Expected chain"));
    std::vector<int> row;int64_t previous=-1;
    for(Py_ssize_t j=0;j<PySequence_Fast_GET_SIZE(chain.p);++j){
     Ref pair(PySequence_Fast(PySequence_Fast_GET_ITEM(chain.p,j),"Expected ordered graph"));
     PyObject* raw_order=PySequence_Fast_GET_ITEM(pair.p,0);
     if(PyBool_Check(raw_order)||!PyLong_Check(raw_order))throw std::invalid_argument("Unique increasing graph order required");
     int64_t order=PyLong_AsLongLong(raw_order);if(PyErr_Occurred())throw Unsupported("Large graph order");
     if(order<0 || order<=previous || !orders.insert(order).second)throw std::invalid_argument("Unique increasing graph order required");
     previous=order;int id=items.size();items.push_back({order,compile(PySequence_Fast_GET_ITEM(pair.p,1))});row.push_back(id);
    }
    if(!row.empty())chains.push_back(std::move(row));
   }
   std::sort(chains.begin(),chains.end(),[&](const auto& a,const auto& b){return items[a.front()].order<items[b.front()].order;});
  }
  next.assign(chains.size(),0);parts.resize(chains.size());usage.assign(capacities.size(),0.);usage_seen.assign(capacities.size(),0);
 }
 Graph& graph(Part& part){return *items[part.item].graph;}
 int64_t& token(Part& part,int id){int shared=graph(part).tokens[id].shared;return shared<0?part.free[id]:shared_free[shared];}
 void admit(int slot){
  if(next[slot]>=chains[slot].size())return;
  auto part=std::make_unique<Part>();part->item=chains[slot][next[slot]++];Graph& g=graph(*part);
  part->degree.reserve(g.nodes.size());part->ends.assign(g.nodes.size(),0.);part->ended.assign(g.nodes.size(),0);
  for(const auto& t:g.tokens)part->free.push_back(t.capacity);
  part->queues.resize(g.queue_ids.size());
  for(size_t i=0;i<g.nodes.size();++i){part->degree.push_back(g.nodes[i].deps.size());if(g.nodes[i].deps.empty())ready.push_back({slot,(int)i});}
  live+=g.nodes.size();added+=g.nodes.size();peak=std::max(peak,live);parts[slot]=std::move(part);
 }
 void release(Part& part,const Node& node,int phase){
  Graph& g=graph(part);
  for(const auto& action:node.actions[phase]){
   int64_t& free=token(part,action.token);
   if(free>g.tokens[action.token].capacity-action.amount)throw std::invalid_argument("Unbalanced token release");
   free+=action.amount;
  }
 }
 void finish(NodeRef ref){
  Part& part=*parts[ref.slot];Graph& g=graph(part);const Node& node=g.nodes[ref.node];
  part.ends[ref.node]=now;part.ended[ref.node]=1;
  int64_t order=items[part.item].order;
  for(int w:node.waits_until)waits.erase({order,w});
  for(int w:node.waits_after)if(!part.ended[g.waits[w].until])waits[{order,w}]={ref.slot,w};
  if(node.enqueue>=0)part.queues[node.enqueue].push_back(node.consumer);
  release(part,node,2);
  for(int child:node.children){if(--part.degree[child]==0)ready.push_back({ref.slot,child});}
  ++finished;++part.finished;
  if(trace || !streamed)ends.push_back({part.item,ref.node,now});
  if(streamed && ref.node==g.complete){
   if(part.finished!=g.nodes.size())throw std::invalid_argument("Cyclic or incomplete graph at chain completion");
   for(size_t w=0;w<g.waits.size();++w)if(waits.count({order,(int)w}))throw std::invalid_argument("Graph completed with an active resource wait");
   live-=g.nodes.size();parts[ref.slot].reset();admit(ref.slot);
  }
 }
 void run(){
  for(size_t i=0;i<chains.size();++i)admit(i);
  std::vector<double> totals(capacities.size()),resource_rates(capacities.size());
  while(!ready.empty() || !active.empty()){
   while(!ready.empty()){
    bool progressed=false;size_t count=ready.size();
    for(size_t i=0;i<count;++i){
     NodeRef ref=ready.front();ready.pop_front();Part& part=*parts[ref.slot];const Node& node=graph(part).nodes[ref.node];
     if(node.dequeue>=0 && (part.queues[node.dequeue].empty() || part.queues[node.dequeue].front()!=ref.node)){
      ready.push_back(ref);continue;
     }
     bool available=true;for(const auto& action:node.actions[0])if(token(part,action.token)<action.amount){available=false;break;}
     if(!available){ready.push_back(ref);continue;}
     for(const auto& action:node.actions[0])token(part,action.token)-=action.amount;
     if(node.dequeue>=0)part.queues[node.dequeue].pop_front();
     release(part,node,1);progressed=true;
     if(trace || !streamed)starts.push_back({part.item,ref.node,now});
     double duration=node.seconds;
     if(node.attempt>=0){
      double waited=part.ends[node.ready]-part.ends[node.attempt];bool did_block=waited>1e-12;
      if(!did_block)duration=0.;
      if(trace || !streamed)conditional.push_back({part.item,ref.node,did_block,std::max(0.,waited),duration});
      if(streamed){++possible;blocked+=did_block;handoff_seconds+=duration;}
     }
     if(duration)active.push_back({ref,duration,1.});else finish(ref);
    }
    if(!progressed)break;
   }
   if(active.empty()){
    if(!ready.empty())throw std::invalid_argument("Token dependency deadlock");
    continue;
   }
   std::fill(totals.begin(),totals.end(),0.);
   // Active nodes retain insertion order, including simultaneous completions.
   for(const auto& a:active){const Node& node=graph(*parts[a.ref.slot]).nodes[a.ref.node];for(const auto& d:node.demands)totals[d.resource]+=d.amount;}
   // std::map uses the same global (graph order, wait order) summation order.
   for(const auto& entry:waits){auto ref=entry.second;const Wait& w=graph(*parts[ref.slot]).waits[ref.wait];for(const auto& d:w.demands)totals[d.resource]+=d.amount;}
   for(size_t r=0;r<totals.size();++r)resource_rates[r]=totals[r]?capacities[r]/totals[r]:1.;
   double interval=std::numeric_limits<double>::infinity();
   for(auto& a:active){
    a.rate=1.;const Node& node=graph(*parts[a.ref.slot]).nodes[a.ref.node];
    for(const auto& d:node.demands)if(d.amount && totals[d.resource])a.rate=std::min(a.rate,resource_rates[d.resource]);
    if(a.rate==0.)throw std::underflow_error("float division by zero");
    interval=std::min(interval,a.remaining/a.rate);
   }
   now+=interval;++steps;
   for(const auto& entry:waits){
    auto ref=entry.second;const Wait& w=graph(*parts[ref.slot]).waits[ref.wait];double rate=1.;
    for(const auto& d:w.demands)if(d.amount && totals[d.resource])rate=std::min(rate,resource_rates[d.resource]);
    for(const auto& d:w.demands){
     if(!usage_seen[d.resource]){usage_seen[d.resource]=1;usage_order.push_back(d.resource);}
     usage[d.resource]+=interval*rate*d.amount;
    }
   }
   std::vector<NodeRef> done;size_t keep=0;
   for(size_t i=0;i<active.size();++i){
    auto a=active[i];a.remaining-=interval*a.rate;
    const Node& node=graph(*parts[a.ref.slot]).nodes[a.ref.node];
    if(a.remaining<=std::max(1e-14,1e-12*node.seconds))done.push_back(a.ref);
    else active[keep++]=a;
   }
   active.resize(keep);
   for(auto ref:done)finish(ref);
  }
  if(finished!=added)throw std::invalid_argument("Cyclic execution dependencies");
 }
 std::string label(int item,int node){
  const auto& entry=items[item];const auto& name=entry.graph->nodes[node].name;
  return streamed?prefix+std::to_string(entry.order)+":"+name:name;
 }
 void put(PyObject* target,const char* key,PyObject* value){
  Ref owned(value);if(PyDict_SetItemString(target,key,owned)<0)throw PythonError();
 }
 PyObject* trace_dict(const std::vector<Trace>& events){
  Ref result(PyDict_New());
  for(const auto& event:events){Ref value(PyFloat_FromDouble(event.time));if(PyDict_SetItemString(result,label(event.item,event.node).c_str(),value)<0)throw PythonError();}
  return Py_NewRef(result.p);
 }
 PyObject* result(){
  Ref out(PyDict_New());put(out,"seconds",PyFloat_FromDouble(now));
  put(out,"start",trace_dict(starts));put(out,"end",trace_dict(ends));
  put(out,"resource_event_steps",PyLong_FromLongLong(steps));
  Ref wait_result(PyDict_New());
  for(int r:usage_order){Ref value(PyFloat_FromDouble(usage[r]));if(PyDict_SetItemString(wait_result,resource_names[r].c_str(),value)<0)throw PythonError();}
  put(out,"wait_resource_seconds",Py_NewRef(wait_result.p));
  Ref delays(PyDict_New());
  for(const auto& value:conditional){
   Ref row(PyDict_New());put(row,"blocked",PyBool_FromLong(value.blocked));put(row,"dependency_wait_seconds",PyFloat_FromDouble(value.waited));
   put(row,"extra_elapsed_service_seconds",PyFloat_FromDouble(value.duration));
   if(PyDict_SetItemString(delays,label(value.item,value.node).c_str(),row)<0)throw PythonError();
  }
  put(out,"conditional_delays",Py_NewRef(delays.p));
  put(out,"resource_policy",PyUnicode_FromString("Proportional active-demand fluid sharing; capacity released at every completion."));
  if(streamed){
   Ref summary(PyDict_New());put(summary,"possible_waits",PyLong_FromLongLong(possible));put(summary,"blocked_waits",PyLong_FromLongLong(blocked));
   put(summary,"extra_elapsed_service_seconds",PyFloat_FromDouble(handoff_seconds));put(out,"conditional_summary",Py_NewRef(summary.p));
   put(out,"peak_active_nodes",PyLong_FromLongLong(peak));put(out,"scheduled_nodes",PyLong_FromLongLong(added));
  }
  return Py_NewRef(out.p);
 }
};

PyObject* solve(PyObject*,PyObject* args){
 PyObject *graph,*chains,*shared;int trace;const char* prefix;
 if(!PyArg_ParseTuple(args,"OOOps",&graph,&chains,&shared,&trace,&prefix))return nullptr;
 try{
  Solver solver(graph,chains,shared,trace,prefix);
  PyThreadState* state=PyEval_SaveThread();
  try{solver.run();}catch(...){PyEval_RestoreThread(state);throw;}
  PyEval_RestoreThread(state);return solver.result();
 }catch(const Unsupported& error){PyErr_Clear();PyErr_SetString(PyExc_NotImplementedError,error.what());}
 catch(const PythonError&){if(!PyErr_Occurred())PyErr_SetString(PyExc_RuntimeError,"Python graph conversion failed");}
 catch(const std::underflow_error& error){PyErr_SetString(PyExc_ZeroDivisionError,error.what());}
 catch(const std::bad_alloc&){PyErr_NoMemory();}
 catch(const std::exception& error){PyErr_SetString(PyExc_ValueError,error.what());}
 return nullptr;
}
PyMethodDef methods[]={{"solve",solve,METH_VARARGS,"Run the reference finite event rules on CPU."},{nullptr,nullptr,0,nullptr}};
PyModuleDef module={PyModuleDef_HEAD_INIT,"_event_solver",nullptr,-1,methods};
}
PyMODINIT_FUNC PyInit__event_solver(){
 PyObject* result=PyModule_Create(&module);if(!result)return nullptr;
 if(PyModule_AddStringConstant(result,"BUILD_KEY",TORCHGWAS_SOLVER_BUILD_KEY)<0){Py_DECREF(result);return nullptr;}return result;
}
