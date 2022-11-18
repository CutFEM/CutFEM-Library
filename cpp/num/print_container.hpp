#ifndef NUM_PRINT_CONTAINER_HPP
#define NUM_PRINT_CONTAINER_HPP

#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

template <typename T>
std::ostream &operator<<(std::ostream &o, const std::vector<T> &v) {
   o << "[";
   if (v.empty()) {
      o << "]";
      return o;
   }
   // For every item except the last write "Item, "
   int i = 0;
   for (auto it = v.begin(); it != --v.end(); it++, ++i) {
      o << *it << ((i % 5) == 4 ? "\n\t" : "\t");
   }
   // Write out the last item
   o << v.back() << "]";
   return o;
}

template <size_t N, typename T>
std::ostream &operator<<(std::ostream &o, const std::array<T, N> &v) {
   o << "[";
   // For every item except the last write "Item, "
   int i = 0;
   for (const auto &a : v) {
      o << a << ((i % 5) == 4 ? "\n\t" : "\t");
      ++i;
   }
   // Write out the last item
   o << "]";
   return o;
}

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &o, const std::pair<T1, T2> &m) {
   o << "{ " << m.first << " , " << m.second << " }";
   return o;
}

template <typename KeyT, typename ValueT>
std::ostream &operator<<(std::ostream &o, const std::map<KeyT, ValueT> &m) {
   o << "{";
   if (m.empty()) {
      o << "}";
      return o;
   }
   // For every pair except the last write "Key: Value, "
   for (auto it = m.begin(); it != --m.end(); it++) {
      const auto &[key, value] = *it;
      o << key << ": " << value << "\n";
   }
   // Write out the last item
   const auto &[key, value] = *--m.end();
   o << key << ": " << value << "}";
   return o;
}

template <typename KeyT, typename ValueT>
std::ostream &operator<<(std::ostream &o,
                         const std::unordered_map<KeyT, ValueT> &m) {
   o << "{";
   if (m.empty()) {
      o << "}";
      return o;
   }
   for (auto it = m.begin(); it != m.end(); it++) {
      const auto &[key, value] = *it;
      o << key << ": " << value << ", ";
   }
   o << "}";
   return o;
}

template <class TupType, size_t... I>
std::ostream &print(std::ostream &o, const TupType &_tup,
                    std::index_sequence<I...>) {
   o << "(";
   (..., (o << (I == 0 ? "" : ", ") << std::get<I>(_tup)));
   o << ")\n";
   return o;
}
template <class... T>
std::ostream &operator<<(std::ostream &o, const std::tuple<T...> &_tup) {
   return print(o, _tup, std::make_index_sequence<sizeof...(T)>());
}

#endif

//
