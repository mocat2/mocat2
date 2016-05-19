// ============================================================================
// popenstream, C++ iostream classes wrapping the popen()
// Copyright (C) 2013  Luis Pedro Coelho
//
// based on gzstream
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library".
// ============================================================================

#ifndef POPENSTREAM_H_INCLUDE_GUARD_
#define POPENSTREAM_H_INCLUDE_GUARD_

// standard C++ with new header file names and std:: namespace
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef POPENSTREAM_NAMESPACE
namespace popen_ns {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement popenstream. See below for user classes.
// ----------------------------------------------------------------------------

class popenstreambuf : public std::streambuf {
private:
    static const int bufferSize = 512;

    FILE*            file;               // file handle for compressed file
    char             buffer_[bufferSize]; // data buffer_
    bool             is_open_;             // open/close state of stream
    int              mode;               // I/O mode

    int flush_buffer();
public:
    popenstreambuf() : is_open_(false) {
        setp(buffer_, buffer_ + (bufferSize-1));
        setg(buffer_ + 4,     // beginning of putback area
             buffer_ + 4,     // read position
             buffer_ + 4);    // end position      
        // ASSERT: both input & output capabilities will not be used together
    }
    bool is_open() { return is_open_; }
    popenstreambuf* open(const char* name, int open_mode);
    popenstreambuf* close();
    ~popenstreambuf() { close(); }
    
    virtual int     overflow( int c = EOF);
    virtual int     underflow();
    virtual int     sync();
};

class popenstreambase : virtual public std::ios {
protected:
    popenstreambuf buf;
public:
    popenstreambase() { init(&buf); }
    popenstreambase(const char* name, int open_mode);
    ~popenstreambase();
    void open(const char* name, int open_mode);
    void close();
    popenstreambuf* rdbuf() { return &buf; }
};


class ipopenstream : public popenstreambase, public std::istream {
public:
    ipopenstream() : std::istream( &buf) {} 
    ipopenstream(const char* name, int open_mode = std::ios::in)
        : popenstreambase( name, open_mode), std::istream( &buf) {}  
    popenstreambuf* rdbuf() { return popenstreambase::rdbuf(); }
    void open( const char* name, int open_mode = std::ios::in) {
        popenstreambase::open(name, open_mode);
    }
};

class opopenstream : public popenstreambase, public std::ostream {
public:
    opopenstream() : std::ostream(&buf) {}
    opopenstream(const char* name, int mode = std::ios::out)
        : popenstreambase(name, mode), std::ostream( &buf) {}  
    popenstreambuf* rdbuf() { return popenstreambase::rdbuf(); }
    void open(const char* name, int open_mode = std::ios::out) {
        popenstreambase::open(name, open_mode);
    }
};

class opopengzstream : public opopenstream {
public:
    opopengzstream(const std::string& name, const std::string& other_options = std::string()) {
        std::ostringstream ss;
        ss << "gzip " << other_options << " > \"" << name << '\"';
        this->open(ss.str().c_str());
    }
};

class ipopengzstream : public ipopenstream {
public:
    ipopengzstream(const std::string& name, const std::string& other_options = std::string()) {
        std::ostringstream ss;
        ss << "gzip -dc " << other_options << " \"" << name << '\"';
        this->open(ss.str().c_str());
    }
};

class opopenbzstream : public opopenstream {
public:
    opopenbzstream(const std::string& name, const std::string& other_options = std::string()) {
        std::ostringstream ss;
        ss << "bzip2 " << other_options << " > \"" << name << '\"';
        this->open(ss.str().c_str());
    }
};

class ipopenbzstream : public ipopenstream {
public:
    ipopenbzstream(const std::string& name, const std::string& other_options = std::string()) {
        std::ostringstream ss;
        ss << "bzip2 -dc " << other_options << " \"" << name << '\"';
        this->open(ss.str().c_str());
    }
};

#ifdef POPENSTREAM_NAMESPACE
} // namespace popen_ns
#endif

#endif // POPENSTREAM_H_INCLUDE_GUARD_
