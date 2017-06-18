#ifndef BUFFER_HPP
#define BUFFER_HPP


#include <boost/circular_buffer.hpp>
#include "data.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#define WAIT_TIME 1000

class buffer{
private:
    bool keep_waiting;
    boost::circular_buffer<Read>* data;
    boost::mutex m;
public:
    buffer(int i_size);
    ~buffer();
    inline int size() const {
        return data->size();
    }
    int put(const Read& read);
    int take(Read& read);
    void stop_waiting();
    void wait(int n) const;
};


#endif
