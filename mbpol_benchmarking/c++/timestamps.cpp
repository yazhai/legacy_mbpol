#include "timestamps.h"
#include <chrono>
#include <string>
#include <iostream>
#include <time.h>
#include <climits>
#include <iomanip>

using namespace std;

timestamp::timestamp(int _threadid, std::string _label):timespan(0), threadid(_threadid),label(_label){};
timestamp::~timestamp(){};
timestamp::timestamp(const timestamp& that){
    start=that.start;
    end=that.end;
    timespan = that.timespan;
    threadid = that.threadid;
    label = that.label;
}



void timestamp::setthread(int _id){
     threadid=_id;
}


void timestamp::setlabel(string _label){
     label = _label;
};

void timestamp::stampstart(){
     start = chrono::high_resolution_clock::now();
};

void timestamp::stampend(){
     end = chrono::high_resolution_clock::now();  
     timespan = (long long int)chrono::duration_cast<chrono::microseconds>(end-start).count();
};
     


timers_t::timers_t(){
     timers_list.clear();
     timecollections.clear();
}

timers_t::~timers_t(){
     timers_list.clear();
     timecollections.clear();
}

timers_t::timers_t(const timers_t& that){
     timers_list.clear();
     for(auto item : that.timers_list){
          timers_list.insert(item);
     }
     timecollections.clear();
     for(auto item : that.timecollections){
          timecollections.insert(item);
     }
}




// insert a timer into the list with a given unique ID
bool timers_t::insert_timer(id _id, int _threadid, string _label){
     if (timers_list.find(_id)==timers_list.end()) {          
          timestamp stamp(_threadid,_label);
          timers_list.insert({_id, stamp});        
          return true;
     }
     return false;
}  

// insert a timer into the list with a random ID, and return with this id.
void timers_t::insert_random_timer(id & id, int _threadid, std::string _label){
     srand(time(NULL));
     id=0;
     bool ifinsert;
     do {
          id = rand() % ULLONG_MAX;  
          ifinsert = this->insert_timer(id, _threadid, _label);         
     } while (!ifinsert);
     return;
}


bool timers_t::timer_start(id _id){
     if (timers_list.find(_id)==timers_list.end()) {
          return false;
     }
     timers_list[_id].stampstart();
     return true;
}


bool timers_t::timer_end(id _id, bool ifadd, bool ifsave){
     if (timers_list.find(_id)==timers_list.end()) {
          return false;
     }
     timers_list[_id].stampend();
     
     if(ifadd){     
          add_time(_id);
     }     
     if(!ifsave){
          timers_list.erase(_id);
     }     
     return true;          
}


long long int timers_t::get_time_span(id _id){
     if (timers_list.find(_id)==timers_list.end()) {
          return 0;
     }
     return timers_list[_id].timespan;
}

int timers_t::get_thread_id(id _id){
     if (timers_list.find(_id)==timers_list.end()) {
          return 0;
     }
     return timers_list[_id].threadid;
}

string timers_t::get_label(id _id){
     if (timers_list.find(_id)==timers_list.end()) {
          return 0;
     }
     return timers_list[_id].label;
}

bool timers_t::get_all_timers_info(){     
     for (auto item : timers_list) {
          cout << setw(4)  << item.second.threadid   
               << setw(10) << item.second.label
               << setw(16) << item.second.timespan
               << endl;
     };    
     return true;
}


int timers_t::get_num_timers(){
     return timers_list.size();
}


void timers_t::get_all_timers_id(){
     for (auto item : timers_list) {
          cout << setw(20) << item.first << endl;
     };
     return;
};

void timers_t::get_time_collections(){
     for(auto itr = timecollections.begin(); itr != timecollections.end(); itr++) {
          cout << " Thread= " << setw(3) << itr->threadid 
               << "    Label= " << setw(10) << itr->label
               << "    Time[ms]= " << setw(15) << itr->timespan
               << endl;
     }
};

bool timers_t::add_time(id _id){
     timestamp& target = timers_list[_id];
     auto findtarget = timecollections.find(target);
     if ( findtarget == timecollections.end()){
          //If not find, initialize a new record in timecollections
          timecollections.insert(target);
     }else {
          //If find, add time to timecollections
          findtarget->timespan += target.timespan;
     }
     return true;
};


