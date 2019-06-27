#ifndef _HTSPAN_RAZERS3_PAIR_HPP_
#define _HTSPAN_RAZERS3_PAIR_HPP_

namespace hts {

struct razers3_pair_f {
	seqan::FragmentStore store;

	razers3_pair_f (const char* ref_fname, const char* seq_fname) {

