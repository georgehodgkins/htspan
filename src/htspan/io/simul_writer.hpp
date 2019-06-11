#include <ostream>
#include <string>
#include <sstream>
#include <list>
#include <fstream>
#include <stdexcept>

#ifndef _HTSPAN_SIMUL_HPP_
#define _HTSPAN_SIMUL_HPP_

namespace hts {

namespace frontend {

using namespace std;

// The rest of the code will access the log and result output via
// the global_log and global_out global objects defined at
// the bottom of the file (in frontend namespace)

//TODO: eliminate string/cstr juggling
class simul_writer {
	private:
	list<ofstream*> file_streams;
	bool u_cerr;
	bool u_cout;
	int verbosity; // verbosity level of writer (default: 1)
	// a value of -1 indicates that verbosity is disabled for the writer
	// a value of 0 indicates total silence
	int curr_vth; // tracker of verbosity level for the current write
	public:
	simul_writer ();
	~simul_writer ();
	void add_file (string);
	void use_cout (bool);
	void use_cerr (bool);
	void set_verbosity (int v) { verbosity = v; }
	simul_writer& operator<< (string);
	simul_writer& operator<< (int);
	simul_writer& operator<< (size_t);
	simul_writer& operator<< (double);
	simul_writer& operator<< (char);
	simul_writer& v (int);
	int v () const {return verbosity;}
};

// construct a writer initially set to be silent
simul_writer::simul_writer () {
	u_cerr = false;
	u_cout = false;
	verbosity = 1;
}

simul_writer::~simul_writer () {
	for (list<ofstream*>::iterator it = file_streams.begin(); it != file_streams.end(); ++it) {
		(*it)->close();
		delete *it;
	}
}

// construct and open an ofstream to the given path, and store a pointer in the list
void simul_writer::add_file (string fname) {
	ofstream *pt = new ofstream(fname.c_str(), ios::trunc);
	if (!pt->good()) {
		throw runtime_error ("Could not open file name passed to simul_writer.");
	}
	file_streams.push_back(pt);
}

// tell the writer whether to use cout (exclusive with cerr)
void simul_writer::use_cout (bool tf) {
	u_cout = tf;
	if (u_cout) {
		u_cerr = false;
	}
}

// tell the writer to use cerr (exclusive with cout)
void simul_writer::use_cerr (bool tf) {
	u_cerr = tf;
	if (u_cerr) {
		u_cout = false;
	}
}

// stream insertion operator alias, writes to all streams and cout or cerr if selected
simul_writer& simul_writer::operator<< (string param) {
	if (verbosity >= curr_vth) {
		for (list<ofstream*>::iterator it = file_streams.begin(); it != file_streams.end(); ++it) {
			(*it)->write(param.c_str(), param.size());//C++98 does not support inserting strings into streams
			(*it)->flush();
		}
		if (u_cerr) {
			cerr << param << flush;
		} else if (u_cout) {
			cout << param << flush;
		}
	}
	return *this;
}

// overloads for common types
// unfortunately, ios manipulators (endl, flush, setw) only work with proper ostreams
simul_writer& simul_writer::operator<< (int param) {
	ostringstream sstream;
	sstream << param;
	return operator<<(sstream.str());
}

simul_writer& simul_writer::operator<< (size_t param) {
	ostringstream sstream;
	sstream << param;
	return operator<<(sstream.str());
}

simul_writer& simul_writer::operator<< (double param) {
	ostringstream sstream;
	sstream << param;
	return operator<<(sstream.str());
}

simul_writer& simul_writer::operator<< (char param) {
	return operator<<(string(1, param));
}

/**
* This function allows code to easily specify verbosity levels for output,
* without external conditional wrappers. There is also a parameterless
* overload of this function which returns the verbosity level, so we can use
* external conditional wrappers where repeated calls to this function would be inefficient
*
* Usage:
* writer_object.v(verbosity_level) << "Message of this verbosity level"
*
* Note that if used it must be used consistently,
* as there is no good way to reset the current threshold automatically.
* 
* @param v_lev Verbosity level of write such that write will only be done if verbosity >= v_lev
* @return reference to the writer object, which can then be streamed to
*/
simul_writer& simul_writer::v (int v_lev) {
	if (verbosity == -1) {
		throw runtime_error("Error: attempting to use v() on a writer which has verbosity disabled.");
	}
	curr_vth = v_lev;
	if (v_lev == 0) {
		operator<<("Writer warning: using a verbosity level of 0 is not recommended.\n");
	}
	return *this;
}

// These are acessible globally, within the frontend namespace
// so any function can easily write to them
// initially defined inactive
simul_writer global_log;// logging output
simul_writer global_out;// result output

}// namespace frontend

}// namespace hts

#endif //_HTSPAN_SIMUL_HPP_