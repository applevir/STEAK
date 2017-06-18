#include "buffer.hpp"


buffer::buffer(int i_size){
    keep_waiting=true;
    data=new boost::circular_buffer<Read>(i_size);
}

buffer::~buffer(){
    delete data;
}

int buffer::put(const Read& read){
    while(data->size()==data->capacity()){
        wait(WAIT_TIME);
    }
    m.lock();
    if(data->size()<data->capacity()){
      data->push_back(read);
      m.unlock();
      return 0;
    }else{
      m.unlock();
      return -1;
    }
}

int buffer::take(Read& read){
    while(keep_waiting && data->empty()){
        wait(WAIT_TIME);
    }
    m.lock();
    if(!data->empty()){
      read=data->front();
      data->pop_front();
      m.unlock();
      return 0;
    }else{
        m.unlock();
        return -1;
    }
}

void buffer::stop_waiting(){
    keep_waiting=false;
}

void buffer::wait(int n) const {
  boost::this_thread::sleep(boost::posix_time::microseconds(n));
}
