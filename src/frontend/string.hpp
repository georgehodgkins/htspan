#ifndef _HTSPAN_CSTRING_HPP_
#define _HTSPAN_CSTRING_HPP_ 

#include <cstring>
#include <string>

// A couple of c-string and std::string utilities

namespace hts {

/**
* Case insensitive strcmp
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
* Counts the occurrences of a character in a string.
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

/**
* Insert appropriate whitespace in a serialized JSON string,
* to make it more human-readable.
*
* Specifically, inserts a tab in and newline after an open curly brace,
* tab out and newline /before/ a closing curly brace, space after a colon,
* and an appropriately indented newline after a comma.
*
* Temporarily resizes the string to ~1.5x its current size, then
* shrinks it to fit when done inserting characters.
*/
void indent_serialized_json (string& s) {
	string indent;
	for (size_t i = 0; i < s.size(); ++i) {

		string ins; // string to be inserted
		bool before = false; // insert before or after position
		s.reserve(s.size() * 3 / 2); // add some extra space for the characters to be inserted

		if (s[i] == '{') {
			indent += '\t'; // one more tab
			ins = "\n" + indent;
		} else if (s[i] == '}') {
			indent.erase(indent.end()-1);
			ins = "\n" + indent;
			before = true;
		} else if (s[i] == ':') {
			ins = " ";
		} else if (s[i] == ',') {
			ins = "\n" + indent;
		}

		if (!ins.empty()) {
			if (before) {
				s.insert(i, ins);
			} else {
				s.insert(i + 1, ins);
			}
			i += ins.size(); // skip the sequence we just inserted
		}
	}

	s.reserve(s.size()); // request shrink-to-fit
}



} // namespace hts

#endif // _HTSPAN_CSTRING_HPP_
