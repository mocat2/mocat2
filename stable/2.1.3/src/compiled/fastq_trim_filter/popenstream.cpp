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
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library".
// ============================================================================

#include "popenstream.h"
#include <iostream>
#include <cassert>
#include <exception>
#include <stdexcept>
#include <string.h>  // for memcpy
#include <errno.h>

#ifdef POPENSTREAM_NAMESPACE
namespace popen_ns {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement gzstream. See header file for user classes.
// ----------------------------------------------------------------------------

// --------------------------------------
// class popenstreambuf:
// --------------------------------------

popenstreambuf* popenstreambuf::open( const char* name, int open_mode) {
    if (is_open()) return 0;
    mode = open_mode;
    // no append nor read/write mode
    if ((mode & std::ios::ate) || (mode & std::ios::app)
        || ((mode & std::ios::in) && (mode & std::ios::out)))
        return 0;
    const char* fmode;
    if (mode & std::ios::in) fmode = "r";
    else if (mode & std::ios::out) fmode = "w";
    else fmode = "r"; // alternatively, we could throw an error here
    file = popen(name, fmode);
    if (!file) {
        throw std::runtime_error("popen call failed");
    }
    is_open_ = true;
    return this;
}

popenstreambuf * popenstreambuf::close() {
    if (is_open()) {
        sync();
        is_open_ = false;
        if (pclose(file) == -1) return 0;
        return this;
    }
    return 0;
    return (popenstreambuf*)0;
}

int popenstreambuf::underflow() { // used for input buffer only
    if (gptr() && ( gptr() < egptr()))
        return * reinterpret_cast<unsigned char *>(gptr());

    if (! (mode & std::ios::in) || !is_open_)
        return EOF;
    // Josuttis' implementation of inbuf
    int n_putback = gptr() - eback();
    if (n_putback > 4) n_putback = 4;
    memcpy(buffer_ + (4 - n_putback), gptr() - n_putback, n_putback);

    int num = fread(buffer_ + 4, 1, bufferSize - 4, file);
    if (num <= 0) // ERROR or EOF
        return EOF;

    // reset buffer pointers
    setg(buffer_ + (4 - n_putback),   // beginning of putback area
         buffer_ + 4,                 // read position
         buffer_ + 4 + num);          // end of buffer

    // return next character
    return * reinterpret_cast<unsigned char *>( gptr());    
}

int popenstreambuf::flush_buffer() {
    // Separate the writing of the buffer from overflow() and
    // sync() operation.
    int w = pptr() - pbase();
    assert(w >= 0);
    if (fwrite(pbase(), 1, w, file) != size_t(w))
        return EOF;
    pbump(-w);
    return w;
}

int popenstreambuf::overflow( int c) { // used for output buffer_ only
    if ( ! ( mode & std::ios::out) || !is_open_)
        return EOF;
    if (c != EOF) {
        *pptr() = c;
        pbump(1);
    }
    if ( flush_buffer() == EOF)
        return EOF;
    return c;
}

int popenstreambuf::sync() {
    // Changed to use flush_buffer() instead of overflow( EOF)
    // which caused improper behavior with std::endl and flush(),
    // bug reported by Vincent Ricard.
    if ( pptr() && pptr() > pbase()) {
        if ( flush_buffer() == EOF)
            return -1;
    }
    return 0;
}

// --------------------------------------
// class popenstreambase:
// --------------------------------------

popenstreambase::popenstreambase(const char* name, int mode) {
    init(&buf);
    open(name, mode);
}

popenstreambase::~popenstreambase() {
    buf.close();
}

void popenstreambase::open( const char* name, int open_mode) {
    if (!buf.open( name, open_mode))
        clear( rdstate() | std::ios::badbit);
}

void popenstreambase::close() {
    if (buf.is_open() && !buf.close()) {
        clear(rdstate() | std::ios::badbit);
    }
}

#ifdef POPENSTREAM_NAMESPACE
} // namespace popen_ns
#endif

// ============================================================================
// EOF //
