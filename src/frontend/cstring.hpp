#ifndef _HTSPAN_CSTRING_HPP_
#define _HTSPAN_CSTRING_HPP_ 

namespace hts {

/**
* Case insensitive strcmp, returns bool instead of the traditional int
*/
int strcmpi (const char* s1, const char* s2) {
	if (s1 == NULL) {
		if (s2 == NULL) {
			// both s1 and s2 are NULL
			return 0;
		} else {
			// s1 is longer
			return 1;
		}
	} else {
		if (s2 == NULL) {
			// s2 is longer
			return -1;
		}
	}

	while (tolower(*s1) == tolower(*s2)) {
		++s1;
		++s2;
		if (*s1 == '\0') break;
		// if *s1 != '\0' and *s2 == '\0',
		// then the while condition will be false
	}

	return *s1 - *s2;
}

}

#endif // _HTSPAN_CSTRING_HPP_
