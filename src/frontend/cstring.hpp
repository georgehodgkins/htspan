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

/**
* Counts the occurences of a character in a string.
*/
size_t str_count_char (const char* s, char c) {
	size_t r = 0;
	for (; *s != '\0'; ++s) {
		if (*s == c) {
			++r;
		}
	}
	return r;
}

} // namespace hts

#endif // _HTSPAN_CSTRING_HPP_
