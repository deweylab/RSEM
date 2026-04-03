#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>
#include <vector>

/**
 * Split a string by delimiters. Replaces boost::tokenizer for C++17 compatibility.
 * @param s The string to split
 * @param delimiters Characters that separate tokens
 * @param keep_empty If true, empty tokens between consecutive delimiters are kept
 */
inline std::vector<std::string> split(const std::string& s, const std::string& delimiters, bool keep_empty = false) {
  std::vector<std::string> result;
  std::string token;
  for (size_t i = 0; i < s.size(); ++i) {
    if (delimiters.find(s[i]) != std::string::npos) {
      if (keep_empty || !token.empty()) {
        result.push_back(token);
      }
      token.clear();
    } else {
      token += s[i];
    }
  }
  if (keep_empty || !token.empty()) {
    result.push_back(token);
  }
  return result;
}

#endif
