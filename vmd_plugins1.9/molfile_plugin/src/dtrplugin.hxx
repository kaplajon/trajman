//
// Version info for VMD plugin tree:
//   $Id: dtrplugin.hxx,v 1.2 2011/03/05 03:24:25 johns Exp $
//
// Version info for last sync with D. E. Shaw Research:
//  //depot/desrad/main/sw/libs/vmd_plugins,DIST/dtrplugin.cxx#3
//

/*
Copyright 2009, D. E. Shaw Research, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research, LLC nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef MOLFILE_DTRPLUGIN_HXX
#define MOLFILE_DTRPLUGIN_HXX

#if defined(_MSC_VER)
#ifndef DESRES_WIN32
#define DESRES_WIN32
#endif
#endif

#include <stdio.h>
#ifdef DESRES_WIN32
#include <io.h>
#include <direct.h>
#include <fcntl.h>
#include <windows.h>
typedef int int32_t;
typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
#if 1
typedef unsigned __int64 uint64_t;    // This also works with MVSC6
#else
typedef unsigned long long uint64_t;
#endif
typedef unsigned short uint16_t;
typedef unsigned int ssize_t;
typedef int mode_t;
#define mkdir(a,b) _mkdir(a)
#define rmdir(a)   _rmdir(a)
#define ftello(a)  ftell(a)
#else
#define O_BINARY 0
#include <inttypes.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#endif

#include <molfile_plugin.h>

#include <vector>
#include <string>


namespace desres { namespace molfile {
  
  const char * dtr_serialized_version();

  typedef struct key_record {
    uint32_t time_lo;       /* Time associated with this file (low bytes). */
    uint32_t time_hi;       /* Time associated with this file (high bytes). */
    uint32_t offset_lo;     /* Zero in the 1 frame/file case (low bytes) */
    uint32_t offset_hi;     /* Zero in the 1 frame/file case (high bytes) */
    uint32_t framesize_lo;  /* Number of bytes in frame (low bytes) */
    uint32_t framesize_hi;  /* Number of bytes in frame (high bytes) */

    double time() const;
    uint64_t offset() const;
    uint64_t size() const;

  } key_record_t;

  class DtrReader;

  class FrameSetReader {
  protected:
    std::string dtr;

  public:
    uint32_t natoms;
    bool with_velocity;

    FrameSetReader()
    : natoms(0), with_velocity(false)
    {}

    virtual ~FrameSetReader() {}

    bool has_velocities() const { return with_velocity; }

    const std::string &path() const { return dtr; }

    // initialize all members from frameset at given path.  Return success.
    // If changed is provided, set to true/false if timekeys were/were not
    // reloaded.
    virtual bool init(const std::string &path, int * changed = NULL) = 0;

    // number of frames
    virtual ssize_t size() const = 0;

    // read the next frame.  If ts is NULL, skip it.  If there are no more
    // frames to read, return MOLFILE_EOF.  Otherwise, fill in coordinates
    // in ts->coords, velocities in ts->velocities if ts->velocities is 
    // non-NULL, and return MOLFILE_SUCCESS if all goes well.
    virtual int next(molfile_timestep_t *ts) = 0;

    // Get the DtrReader component corresponding to frame n.  Change
    // n to the local index within the returned dtr.
    virtual const DtrReader * component(ssize_t &n) const = 0;

    // read a specified frame.
    virtual int frame(ssize_t n, molfile_timestep_t *ts) const = 0;

    // read up to count times beginning at index start into the provided space;
    // return the number of times actually read.
    virtual ssize_t times(ssize_t start, ssize_t count, double * times) const = 0;

    virtual std::ostream& dump(std::ostream &out) const = 0;
    virtual std::istream& load(std::istream &in) = 0;
  };

  class DtrReader : public FrameSetReader {
    uint32_t frames_per_file;
    mutable int m_ndir1;
    mutable int m_ndir2;
    ssize_t m_curframe;
    float * rmass;
    bool owns_rmass;

    bool eof() const { return m_curframe >= (ssize_t)keys.size(); }

  public:
    DtrReader() 
    : frames_per_file(0), m_ndir1(-1), m_ndir2(-1), m_curframe(0),
      rmass(NULL), owns_rmass(false) {}

    virtual ~DtrReader() {
      set_rmass(NULL);
    }

    std::vector<key_record_t> keys;  // one per frame

    // lazy-loaded DDparams
    int ndir1() const;
    int ndir2() const;

    /* rmass and owns_rmass are initially set to NULL and false.  In init(), 
     * if rmass is NULL and owns_rmass is false, we try to read the rmass
     * from the metadata frame (an expensive operation), in which case 
     * owns_rmass becomes true, whether or not rmass was actually read.
     * Otherwise, we leave those values alone.  In the destructor, we delete 
     * rmass if we own in.  The StkReader class can share the rmass between 
     * DtrReader instances in the following way: if the rmass it has to share
     * is NULL, it should set owns_rmass to true in the DtrReader instances,
     * so that the DtrReaders don't keep search for their own rmass.  If
     * the rmass it has to share is non-NULL, it should set owns_rmass to 
     * false so that the rmass pointer doesn't get double-freed.
     */
    float * get_rmass() const { return rmass; }
    void set_rmass(float * ptr) {
      if (rmass && owns_rmass) delete [] rmass;
      if (ptr) {
        rmass = ptr;
        owns_rmass = false;
      } else {
        rmass = NULL;
        owns_rmass = true;
      }
    }

    /* the keys we removed in order to fix overlapping frames.  These get
     * serialized/deserialized by dump() and load().
     */
    std::vector<key_record_t> overlapping_keys;

    ssize_t curframe() const { return m_curframe; }

    uint32_t framesperfile() const { return frames_per_file; }

    virtual bool init(const std::string &path, int * changed=NULL);
    virtual ssize_t size() const { return keys.size(); }
    virtual ssize_t times(ssize_t start, ssize_t count, double * times) const;

    virtual int next(molfile_timestep_t *ts);

    virtual const DtrReader * component(ssize_t &n) const {
      return this;
    }

    virtual int frame(ssize_t n, molfile_timestep_t *ts) const;


    // path for frame at index.  Empty string on not found.
    std::string framefile(ssize_t n) const;

    // read a frame from supplied bytes
    int frame_from_bytes( const void *buf, uint64_t len,
                          molfile_timestep_t *ts ) const;

    std::ostream& dump(std::ostream &out) const;
    std::istream& load(std::istream &in);
  };

  struct DtrWriter {
    std::string dtr;
    const uint32_t natoms;
    int frame_fd;        // for writing
    std::vector<key_record_t> keys;  // one per frame
    uint32_t frames_per_file;

    explicit DtrWriter(uint32_t natoms_) 
    : natoms(natoms_), frame_fd(0), frames_per_file(256)
    {}

    ~DtrWriter();

    // initialize for writing at path
    bool init(const std::string &path);

    // write another frame.  Return MOLFILE_SUCCESS or MOLFILE_ERROR.
    int next(const molfile_timestep_t *ts);
  };

  class StkReader : public FrameSetReader {
    std::vector<DtrReader*> framesets;
    size_t curframeset;

  public:
    StkReader() : curframeset(0) {}
    StkReader(DtrReader *reader);
    ~StkReader();
    virtual bool init(const std::string &path, int * changed=NULL);
    virtual ssize_t size() const;
    virtual ssize_t times(ssize_t start, ssize_t count, double * times) const;
    virtual int next(molfile_timestep_t *ts);
    virtual int frame(ssize_t n, molfile_timestep_t *ts) const;
    virtual const DtrReader * component(ssize_t &n) const;

    static bool recognizes(const std::string &path);

    std::ostream& dump(std::ostream &out) const;
    std::istream& load(std::istream &in);
  };
} }

#endif
