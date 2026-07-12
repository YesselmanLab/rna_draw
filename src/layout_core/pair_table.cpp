/**
 * @file pair_table.cpp
 * @brief Implements `make_pair_table` (`pair_table.hpp`).
 */

#include "rna_layout/pair_table.hpp"

#include <stdexcept>

namespace rna_layout {

std::vector<int> make_pair_table(const std::string& dot_bracket) {
  std::vector<int> pair_table(dot_bracket.size() + 1, 0);
  pair_table[0] = static_cast<int>(dot_bracket.size());

  std::vector<int> open_stack;
  open_stack.reserve(dot_bracket.size());

  for (std::size_t index = 0; index < dot_bracket.size(); ++index) {
    const char base = dot_bracket[index];
    const int position = static_cast<int>(index) + 1;  // 1-indexed, matching vrna_ptable

    if (base == '(') {
      open_stack.push_back(position);
    } else if (base == ')') {
      if (open_stack.empty()) {
        throw std::invalid_argument("unbalanced structure: an unmatched ')' at position " +
                                    std::to_string(position));
      }
      const int partner = open_stack.back();
      open_stack.pop_back();
      pair_table[position] = partner;
      pair_table[partner] = position;
    } else if (base != '.') {
      throw std::invalid_argument(
          "structure contains a character outside \"().\": not well-nested dot-bracket input");
    }
  }

  if (!open_stack.empty()) {
    throw std::invalid_argument("unbalanced structure: an unmatched '(' would never close");
  }
  return pair_table;
}

}  // namespace rna_layout
