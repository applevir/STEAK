// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#ifndef BUFFER_HPP
#define BUFFER_HPP


#include <boost/circular_buffer.hpp>
#include "data.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <vector>
#include <iostream>

#define WAIT_TIME 1000

template<class T> class buffer{
private:
    bool keep_waiting;
    boost::circular_buffer<T>* data;
    boost::mutex m;
public:
    buffer(int i_size){
      keep_waiting=true;
      data=new boost::circular_buffer<T>(i_size);
    }
    ~buffer(){
      delete data;
    }
    inline int size() const {
        return data->size();
    }
    int put(const T& t){
      while(data->size()==data->capacity()){
        wait(WAIT_TIME);
      }
      m.lock();
      if(data->size()<data->capacity()){
        data->push_back(t);
        m.unlock();
        return 0;
      }else{
        m.unlock();
        return -1;
      }
    }


    int put(const std::vector<T>& v){
      if(v.empty()){
	return 0;
      }
      while(data->size()>=data->capacity()-v.size()){
        wait(WAIT_TIME);
      }
      m.lock();
    if(data->size()<data->capacity()){
      for(int i=0; i<v.size(); ++i){
        data->push_back(v[i]);
      }
      m.unlock();
      return 0;
    }else{
      m.unlock();
      return -1;
    }
}
    int take(T& t){
    while(keep_waiting && data->empty()){
        wait(WAIT_TIME);
    }
    m.lock();
    if(!data->empty()){
      t=data->front();
      data->pop_front();
      m.unlock();
      return 0;
    }else{
        m.unlock();
        return -1;
    }
}

    int take(std::vector<T>& v) {
      v.clear();
  while(keep_waiting && data->empty()){
        wait(WAIT_TIME);
    }
    m.lock();
    if(!data->empty()){
      while(!data->empty()){
        v.push_back(data->front());
        data->pop_front();
      }
      m.unlock();
      return 0;
    }else{
        m.unlock();
        return -1;
    }
} 
    void stop_waiting(){
    keep_waiting=false;
}

    void wait(int n) const{
  boost::this_thread::sleep(boost::posix_time::microseconds(n));
}

};


#endif
