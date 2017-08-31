#ifndef TIMESTAMPS_H
#define TIMESTAMPS_H


#include <chrono>
#include <map>
#include <set>
#include <string>


struct timestamp {
     std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
     std::string label;
     mutable long long int timespan;
     int threadid;
     
     timestamp(int id=0, std::string label="");
     ~timestamp();
     timestamp(const timestamp& that);
     //timestamp operator=(const timestamp& that);
     
     void setthread(int id);
     void setlabel(std::string label);
     void stampstart();
     void stampend();
};


struct compare_by_label {
     bool operator () (const timestamp & lhs, const timestamp & rhs){
          if (lhs.threadid < rhs.threadid) {
               return true;
          } else if ((lhs.threadid == rhs.threadid) && (lhs.label < rhs.label)){
               return true;
          }
          return false;
     };
};




typedef unsigned long long int id;

class timers_t {
private:
     std::map<id, timestamp> timers_list;
     std::set<timestamp, compare_by_label> timecollections;     
     
public:     
     
     timers_t();
     ~timers_t();
     timers_t(const timers_t& that);

// Following functions are used to deal with a timer according to its unique id:     
     bool insert_timer(id _id, int threadid=0, std::string _label="");
     void insert_random_timer(id & _id, int threadid=0, std::string _label="");
     bool timer_start(id _id);
     bool timer_end(id _id, bool ifadd=true, bool ifsave=false);
     
// Get timer information based on the id:
     long long int get_time_span(id _id);
     int get_thread_id(id _id);
     std::string get_label(id _id);
     

// Get info:     
     int get_num_timers();
     bool get_all_timers_info();
     void get_all_timers_id();
     void get_time_collections();
     
// Add time to collection according to thread and label
     bool add_time(id _id);     
          
};

#endif


