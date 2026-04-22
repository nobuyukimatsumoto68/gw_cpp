#include <highfive/highfive.hpp>

int main() {
  auto file = HighFive::File("foo.h5", HighFive::File::Create);
  return 0;
}
